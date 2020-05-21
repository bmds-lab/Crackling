'''

Faster and better CRISPR guide RNA design with the Crackling method.
Jacob Bradford, Timothy Chappell, Dimitri Perrin
bioRxiv 2020.02.14.950261; doi: https://doi.org/10.1101/2020.02.14.950261


Purpose:    identify all offtarget sites in the whole genome

Input:      one file with all the chromosome sequences, one chromosome per line

Output:     one file with all the sites

To use:     python prepareListOfftargetSites.py <input-file> <output-file>
            python prepareListOfftargetSites.py ~/genomes/mouse.txt ~/genomes/mouse_offtargets.txt


'''

from time import localtime, strftime
import re,string, sys

# Defining the patterns used to detect sequences
pattern_forward_offsite = r"(?=([ACG][ACGT]{19}[ACGT][AG]G))"
pattern_reverse_offsite = r"(?=(C[CT][ACGT][ACGT]{19}[TGC]))"

# Input and output
fileInput = None
fileOutput = None

# Function that returns the reverse-complement of a given sequence
def rc(dna):
    complements = string.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    rcseq = dna.translate(complements)[::-1]
    return rcseq

# Function that finds off-targets and writes them to the specified file
def findOfftargets():
    print(strftime("%H:%M:%S", localtime())+":\tRunning for: "+fileInput)
    print(strftime("%H:%M:%S", localtime())+":\tWriting to: "+fileOutput)

    with open(fileInput, 'r') as inFile, open(fileOutput, 'w+') as outFile:
        # For every chromosome
        lineNumber = 0
        for line in inFile:

            print strftime("%H:%M:%S", localtime())+":\tChromomose "+str(lineNumber)
            
            # we parse the line and look for forward sequences
            print(strftime("%H:%M:%S", localtime())+":\t\tForward-parsing the chromosome.")
            match_chr = re.findall(pattern_forward_offsite, line)
            print(strftime("%H:%M:%S", localtime())+":\t\tProcessing and saving the parsed sequences.")
            if match_chr:
                print("\t\t\tWe detected "+str(len(match_chr))+" possible off-target sites.")
                # we save each sequence
                for i in range(0,len(match_chr)):
                    outFile.write(match_chr[i][0:20]+"\n")
            else:
                print "\t\t\tWe did not detect any possible off-target sites."

            # we parse the line and look for reverse sequences
            print(strftime("%H:%M:%S", localtime())+":\t\tReverse-parsing the chromosome.")
            match_chr = re.findall(pattern_reverse_offsite, line)
            print(strftime("%H:%M:%S", localtime())+":\t\tProcessing the parsed sequences.")
            if match_chr:
                print("\t\t\tWe detected "+str(len(match_chr))+" possible off-target sites.")
                # we save each reverse-complement sequence
                for i in range(0,len(match_chr)):
                    # we reverse-complement the sequence and count the number of mismatches between the rc and the target
                    outFile.write(rc(match_chr[i])[0:20]+"\n")

            lineNumber += 1

    print("\n"+strftime("%H:%M:%S", localtime())+":\tDone.")

if __name__ == '__main__':
    if (len(sys.argv) != 3):
        print('Error!')
        print('Expecting: prepareListOfftargetSites.py <input-file> <output-file>')
        print('\n')
        exit()

    fileInput = sys.argv[1]
    fileOutput = sys.argv[2]
    
    findOfftargets()