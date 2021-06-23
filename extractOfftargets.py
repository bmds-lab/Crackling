'''
Faster and better CRISPR guide RNA design with the Crackling method.
Jacob Bradford, Timothy Chappell, Dimitri Perrin
bioRxiv 2020.02.14.950261; doi: https://doi.org/10.1101/2020.02.14.950261


Purpose:    identify all offtarget sites in the whole genome

Input:      FASTA, or multi-FASTA, formatted file

Output:     one file with all the sites

To use:     python3.7 ExtractOfftargets.py output-file  (input-files... | input-dir>)

'''

import glob, multiprocessing, os, re, shutil, string, sys, tempfile
from Helpers import *
from Paginator import Paginator

# Defining the patterns used to detect sequences
pattern_forward_offsite = r"(?=([ACG][ACGT]{19}[ACGT][AG]G))"
pattern_reverse_offsite = r"(?=(C[CT][ACGT][ACGT]{19}[TGC]))"

# The off-target sites need to be sorted so the ISSL index is space-optimised.
# In some cases, there are many files to sort, we can paginate the files.
# How many files should we sort per page?
# This is dictated by the number of arguments that we can pass to `sort`
# Default: 50
SORT_PAGE_SIZE = 500

# Set the number of processes to generate.
# This sets the --threads argument for `sort`, and
# the process count for multiprocessing.Map(..)
# Default: os.cpu_count()
PROCESSES_COUNT = os.cpu_count()

def processingNode(fpInputs, fpOutputTempDir = None):
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
                
            outFile.write(
                '\n'.join(offtargets)
            )

def explodeMultiFastaFile(fpInput, fpOutputTempDir):
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

def bashPaginatedSort(filesToSort):
    # sort in batches
    unsortedFiles = filesToSort

    while len(unsortedFiles) > 1:
        newUnsortedFiles = []
        
        onlyMerge = False
        if len(newUnsortedFiles) > 0:
            onlyMerge = True
        
        for pageNum, pageContents in Paginator(
            unsortedFiles,
            SORT_PAGE_SIZE
        ):
            strFiles = '"' + '" "'.join(pageContents) + '"'

            fpTemp = tempfile.NamedTemporaryFile(
                mode = 'w+', 
                delete = False
            )

            if onlyMerge:
                caller(f'sort --merge --parallel={PROCESSES_COUNT} {strFiles} > {fpTemp.name}', shell=True)
            else:
                caller(f'sort --parallel={PROCESSES_COUNT} {strFiles} > {fpTemp.name}', shell=True)
            
            newUnsortedFiles.append(fpTemp.name)

            for fp in pageContents:
                os.remove(fp)
            
        unsortedFiles = newUnsortedFiles
        onlyMerge = True
            
    # should only contain one list that is sorted, despite the name
    return unsortedFiles 

def startMultiprocessing(fpInputs, fpOutput):
    printer('Extracting off-targets using multiprocessing approach')
    
    printer(f'Allowed processes: {PROCESSES_COUNT}')
    
    fpTempDir = tempfile.TemporaryDirectory()
    printer(f'Created a temporary directory for intermediate files: {fpTempDir.name}')

    if len(fpInputs) == 1 and os.path.isdir(fpInputs[0]):
        fpInputs = glob.glob(
            os.path.join(
                fpInputs[0], 
                '*'
            )
        )

    if len(fpInputs) == 1:
        printer('Only one input file to process')
        
        fpExplodeTempDir = tempfile.TemporaryDirectory()
        printer('Attempting to explode multi-FASTA file')
        printer(f'Writing each file to a temporary directory: {fpExplodeTempDir.name}')
        
    
        fpInputs = explodeMultiFastaFile(
            fpInputs[0],
            fpExplodeTempDir.name
        )
        
        printer(f'Exploded into {len(fpInputs)} files')

    args = [
        (
            [fpInput],
            fpTempDir.name 
        ) for fpInput in fpInputs
    ]

    printer(f'Beginning to process {len(args)} files...')

    with multiprocessing.Pool(PROCESSES_COUNT) as mpPool:
        # https://docs.python.org/3/library/multiprocessing.html#multiprocessing.pool.Pool.starmap
        mpPool.starmap(
            processingNode,
            args
        )

    printer('Processing completed')
    
    printer('Preparing for ISSL by sorting all intermediate files')
    printer('Then, writing to user-specified output file')
    
    # sort in batches
    fpSorted = bashPaginatedSort(
        glob.glob(
            os.path.join(
                fpTempDir.name, 
                '*'
            )
        )
    )
    
    shutil.move(fpSorted[0], fpOutput)

if __name__ == '__main__':
    if (len(sys.argv) < 3):
        print('Error!')
        print('Expecting: ExtractOfftargets.py <output-file> [<input-file-1> <input-file-2> <input-file-n> | <input-dir>]')
        print('\n')
        exit()

    fpOutput = sys.argv[1]
    
    fpInputs = sys.argv[2:]
        
    startMultiprocessing(fpInputs, fpOutput)
    
    printer('Goodbye.')
