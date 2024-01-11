#!/usr/bin/python3.6
'''
Faster and better CRISPR guide RNA design with the Crackling method.
Jacob Bradford, Timothy Chappell, Dimitri Perrin
bioRxiv 2020.02.14.950261; doi: https://doi.org/10.1101/2020.02.14.950261


Purpose:    identify all offtarget sites in the whole genome

Input:      FASTA, or multi-FASTA, formatted file

Output:     one file with all the sites

To use:     python3.7 ExtractOfftargets.py output-file  (input-files... | input-dir>)

'''

import glob, multiprocessing, os, re, sys, tempfile, heapq, psutil, argparse
from crackling.Helpers import *
from crackling.Paginator import Paginator

# Defining the patterns used to detect sequences
pattern_forward_offsite = r"(?=([ACG][ACGT]{19}[ACGT][AG]G))"
pattern_reverse_offsite = r"(?=(C[CT][ACGT][ACGT]{19}[TGC]))"

def _explodeMultiFastaFile(fpInput, fpOutputTempDir):
    newFilesPaths = []

    with open(fpInput, 'r') as fRead:
        fWrite = None

        for line in fRead:
            line = line.strip()
            
            # just found a new fasta segment. open a new file
            if line[0] == '>':
                fpTemp = tempfile.NamedTemporaryFile(
                    mode = 'w+', 
                    delete = False,
                    dir = fpOutputTempDir
                )
                
                # close the current file if necessary
                if fWrite is not None:
                    fWrite.write('\n')
                    fWrite.close()
                    
                # open a new one
                fWrite = open(fpTemp.name, 'w+')
                newFilesPaths.append(fpTemp.name)
                
                fWrite.write(line)
                fWrite.write('\n')

            if line[0] != '>':
                fWrite.write(line.upper().strip())
        
        # close the last file if necessary
        if fWrite is not None:
            fWrite.close()
            
    return newFilesPaths

def _splitListOfFilesBySize(listOfFiles, maxBlockSize=None, sizeFactor=1.0):

    if maxBlockSize == None:
        maxBlockSize = psutil.virtual_memory().available # bytes

    pages = [[]]
    
    bufferedSizes = [os.path.getsize(x) * sizeFactor for x in listOfFiles]

    for fileNumber, bufferedSize in enumerate(bufferedSizes):
        predictedBlockSize = sum([
            bufferedSizes[i]
            for i in pages[-1]
        ], bufferedSize)
        
        if predictedBlockSize >= maxBlockSize and len(pages[-1]) != 0:
            pages.append([])
        pages[-1].append(fileNumber)
        
    return [
        [listOfFiles[i] for i in page]
        for page in pages
    ]

def _processingNode(fpInputs, fpOutputTempDir = None):
    # Create a temporary file
    fpTemp = tempfile.NamedTemporaryFile(
        mode = 'w+', 
        delete = False,
        dir = fpOutputTempDir
    )
    
    with open(fpTemp.name, 'w+') as outFile:

        for fpInput in fpInputs:
        
            # key: FASTA header, value: sequence
            seqsByHeader = {}

            with open(fpInput, 'r') as inFile:
                header = fpInput
                
                for line in inFile:
                    if line[0] == '>':
                        header = line[1:]
                        seqsByHeader[header] = []
                    else:
                        if header not in seqsByHeader:
                            # it could be a plain text file, without a header
                            seqsByHeader[header] = []
                        seqsByHeader[header].append(line.rstrip().upper())

            # For each FASTA sequence
            lineNumber = 0
            offtargets = [] 
            for header in seqsByHeader:
                seq = ''.join(seqsByHeader[header])
                
                for strand, pattern, seqModifier in [
                    ['positive', pattern_forward_offsite, lambda x : x],
                    ['negative', pattern_reverse_offsite, lambda x : rc(x)]
                ]:
                    match_chr = re.findall(pattern, seq)

                    for i in range(0,len(match_chr)):
                        offtargets.append(
                            seqModifier(match_chr[i][0:20])
                        )

                lineNumber += 1
                
            outFile.write(''.join(f'{offTarget}\n' for offTarget in offtargets))

def _sortSites(filesToSort, fpOutput, mpPool): 
    printer('Preparing for ISSL by sorting all intermediate files')
    printer('Then, writing to user-specified output file')
    
    # Create temp file directory
    sortedTempDir = tempfile.TemporaryDirectory()
    printer(f'Created temp directory {sortedTempDir.name} for sorting')

    # If the input is a directory then select all files
    if os.path.isdir(filesToSort):
        filesToSort = glob.glob(
            os.path.join(
                filesToSort, 
                '*'
            )
        )
        
    # Split sorting the files into pages to avoid saturating memory
    # Sizefactor: in Python, 70 bytes are needed for a 20 character string
    pages = _splitListOfFilesBySize(filesToSort, maxBlockSize=None, sizeFactor=4.0)
    
    if len(pages) > 1:
        printer(f'Files will be sorted in {len(pages)} pages to avoid saturating memory')

    for i, page in enumerate(pages):
        if len(pages) > 1:
            printer(f'Starting on page {i+1} of {len(pages)}, containing {len(page)} files')
        else:
            printer(f'Beginning to sort {len(page)} files')

        args = [
            (
                file,
                sortedTempDir.name 
            ) for file in page
        ]

        mpPool.starmap(
            _sortingNode,
            args
        )
        
    # Collect sorted files to merge
    sortedFiles = glob.glob(
        os.path.join(
            sortedTempDir.name,
            '*'
        )
    )
    
    printer(f'Beginning to merge {len(sortedFiles)} sorted files')
    
    # Open all the sorted files to merge
    sortedFilesPointers = [open(file, 'r') for file in sortedFiles]

    # Merge all the sorted files 
    with open(fpOutput, 'w') as f:
        f.writelines(heapq.merge(*sortedFilesPointers))
    
    # Close all of the sorted files
    for file in sortedFilesPointers:
        file.close()

def _sortingNode(fileToSort, sortedTempDir):
    # Create a temporary file
    sortedFile = tempfile.NamedTemporaryFile(
        mode = 'w+', 
        delete = False,
        dir = sortedTempDir
    )
    # Sort input file and store in new output dir
    with open(fileToSort, 'r') as input:
        lines = input.readlines()
        lines.sort()
        sortedFile.writelines(lines)
        sortedFile.close()

def _startProcessing(fpInputs, fpOutput, numProcessors, sizeFactor=1.0):
    printer('Extracting off-targets using multiprocessing approach')
    printer(f'Allowed processes: {numProcessors} (system has {os.cpu_count()})')
    
    mpPool = multiprocessing.Pool(numProcessors)
    
    fpTempDir = tempfile.TemporaryDirectory()
    printer(f'Created a temporary directory for intermediate files: {fpTempDir.name}')

    # Find all files in the input directory (if not a list of files)
    if len(fpInputs) == 1 and os.path.isdir(fpInputs[0]):
        fpInputs = glob.glob(
            os.path.join(
                fpInputs[0], 
                '*'
            )
        )
        
    # If there is only one file to process, split it up
    if len(fpInputs) == 1:
        fpExplodeTempDir = tempfile.TemporaryDirectory()
        
        printer('Only one input file to process')
        printer('Attempting to explode multi-FASTA file')
        printer(f'Writing each file to a temporary directory: {fpExplodeTempDir.name}')

        fpInputs = _explodeMultiFastaFile(
            fpInputs[0],
            fpExplodeTempDir.name
        )
        
        printer(f'Exploded into {len(fpInputs)} files')
    
    else:
        printer(f'Processing {len(fpInputs)} files')

    # Split processing the files into pages to avoid saturating memory
    pages = _splitListOfFilesBySize(fpInputs, sizeFactor=sizeFactor)
    
    if len(pages) > 1:
        printer(f'Files will be multiprocessed in {len(pages)} pages to avoid saturating memory')

    for i, page in enumerate(pages):
        if len(pages) > 1:
            printer(f'Starting on page {i+1} of {len(pages)}, containing {len(page)} files')
        else:
            printer(f'Beginning to process {len(page)} files')

        print(page)

        args = [
            (
                [fpInput],
                fpTempDir.name 
            ) for fpInput in page
        ]
    
    
        mpPool.starmap(
            _processingNode,
            args
        )

    printer('Processing completed')
    
    _sortSites(fpTempDir.name, fpOutput, mpPool)
    
    mpPool.close()

def extractOfftargets(output, inputs, memory=None, numProcessors=None, sizeFactor=None):
    if memory is None:
        memory = psutil.virtual_memory().available
        
    if numProcessors is None:
        numProcessors = os.cpu_count()
        
    if sizeFactor is None:
        '''
        sizeFactor is defaulted to 13.0
        
        Experiments using artificial input sequences with lengths between 10k
        and 300k show that approx. 13x the memory of the input sequence is 
        needed. i.e. the number of sites * 70 bytes (50 overhead + 20 seq bytes)
        
        Based on the 12 genomes used in (Bradford, 2022; CRISPR Journal), this
        experimental approximation matches real data.
        '''
        sizeFactor = 30.0
        
    _startProcessing(inputs, output, numProcessors, sizeFactor=sizeFactor)
    
def main():
    parser = argparse.ArgumentParser(description='Extract CRISPR-Cas9 target sites. Used primarily for creating an ISSL index of off-target sites.')
    parser.add_argument('output', help='The output file')
    parser.add_argument('inputs', help='The FASTA-formatted files to extract sites from', nargs='+')
    parser.add_argument('-m', '--memory', required=False, default=None,
        help='The maximum volume of memory to consume. Adjust `--sizeFactor` to scale the effective volume of memory available. Note, if one file exceeds the available memory, it will be processed by itself. If not specified, the current available volume of memory will be used.')
    parser.add_argument('-p', '--processors', required=False, default=None, type=int,
        help='The number of CPU cores to use in multiprocessing. Default is all cores.')
    parser.add_argument('--sizefactor', required=False, default=None, type=float,
        help='Input file sizes will be multiplied by this value. Larger values imply fewer files will be processed simultaneously. Note: Python has an overhead on strings, of 50 bytes; therefore, each extracted site (20 bp) needs 70 bytes in memory.')
    
    args = parser.parse_args()
    
    extractOfftargets(
        args.output,
        args.inputs,
        memory=args.memory,
        numProcessors=args.processors,
        sizeFactor=args.sizefactor
    )

    printer('Goodbye.')

if __name__ == '__main__':
    main()