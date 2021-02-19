'''
Faster and better CRISPR guide RNA design with the Crackling method.
Jacob Bradford, Timothy Chappell, Dimitri Perrin
bioRxiv 2020.02.14.950261; doi: https://doi.org/10.1101/2020.02.14.950261


Purpose:    identify all offtarget sites in the whole genome

Input:      FASTA, or multi-FASTA, formatted file

Output:     one file with all the sites

To use:     python3.7 extractOfftargets.py <output-file> <input-file-1> <input-file-2> <input-file-n>

'''

import multiprocessing, os, re, string, sys, tempfile
from Helpers import *

# Defining the patterns used to detect sequences
pattern_forward_offsite = r"(?=([ACG][ACGT]{19}[ACGT][AG]G))"
pattern_reverse_offsite = r"(?=(C[CT][ACGT][ACGT]{19}[TGC]))"

def processingNode(fpInputs, fpOutput = None, fpOutputTempDir = None):
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
                previousHeader = fpInput
                
                for line in inFile:
                    line = line.strip()
                    if line[0] == '>':
                        seqsByHeader[line[1:]] = ""
                        previousHeader = line[1:]
                        continue
                    else:
                        if previousHeader not in seqsByHeader:
                            # it could be a plain text file, without a header
                            seqsByHeader[previousHeader] = ""
                            
                        seqsByHeader[previousHeader] += line.strip()
            
        
            # For each FASTA sequence
            lineNumber = 0
            for header in seqsByHeader:
            
                seq = seqsByHeader[header]
                
                for strand, pattern, seqModifier in [
                    ['positive', pattern_forward_offsite, lambda x : x],
                    ['negative', pattern_reverse_offsite, lambda x : rc(x)]
                ]:
                
                    match_chr = re.findall(pattern, line)

                    for i in range(0,len(match_chr)):
                        outFile.write(
                            seqModifier(match_chr[i][0:20])+'\n'
                        )

                lineNumber += 1

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

def startMultiprocessing(fpInputs, fpOutput):
    printer('Extracting off-targets using multiprocessing approach')
    
    printer(f'CPU count: {os.cpu_count()}')
    
    printer('Creating a temporary directory for intermediate files')

    fpTempDir = tempfile.TemporaryDirectory()

    if len(fpInputs) == 1:
        printer('Only one input file to process')
        printer('Attempting to separate file, if it is multi-FASTA formatted')
        
    
        fpExplodeTempDir = tempfile.TemporaryDirectory()
        fpInputs = explodeMultiFastaFile(
            fpInputs[0],
            fpExplodeTempDir.name
        )
        print(fpInputs)

    args = [
        (
            [fpInput],
            None, 
            fpTempDir.name 
        ) for fpInput in fpInputs
    ]

    printer('Beginning to process...')
    
    with multiprocessing.Pool() as mpPool:
        # https://docs.python.org/3/library/multiprocessing.html#multiprocessing.pool.Pool.starmap
        mpPool.starmap(
            processingNode,
            args
        )

    printer('Processing completed')
    
    printer('Preparing for ISSL by sorting all intermediate files')
    printer('Then, writing to user-specified output file')
    
    fpTempDirAllFiles = os.path.join(fpTempDir.name, '*')
    numberOfCores = os.cpu_count()
    caller(f'sort --parallel=64 {fpTempDirAllFiles} > {fpOutput}', shell=True)

if __name__ == '__main__':
    if (len(sys.argv) < 3):
        print('Error!')
        print('Expecting: prepareListOfftargetSites.py <output-file> <input-file-1> <input-file-2> <input-file-n>')
        print('\n')
        exit()

    fpOutput = sys.argv[1]
    
    fpInputs = sys.argv[2:]
    
    startMultiprocessing(fpInputs, fpOutput)
    
    printer('Goodbye.')
