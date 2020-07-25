'''
https://github.com/bmds-lab/Crackling

Author: Jake Bradford, Dimitri Perrin

This is an implementation of three CRISPR-Cas9 guide design methods:
    1. mm10db (20 >= AT% <= 65, secondary structure energy, polythymine) 
    2. CHOPCHOP (G20)
    3. sgRNAScorer 2.0 (efficacy model trained on Chari-2015 data)
    
It implements our consensus approach where at least two of the tools need to be
in agreeance for the guide to be accepted. See references list below.

Config:
    - See config.ini

References:
    - https://github.com/bmds-lab/mm10db
    
    - https://chopchop.cbu.uib.no/
    
    - https://crispr.med.harvard.edu/sgRNAScorer/
    
    - Bradford, J., & Perrin, D. (2019). Improving CRISPR guide design with consensus approaches. BMC genomics, 20(9), 931.
    
    - Bradford, J., & Perrin, D. (2019). A benchmark of computational CRISPR-Cas9 guide design methods. PLoS Computational Biology, 15(8), e1007274.
    
    - Chari, R., Yeo, N. C., Chavez, A., & Church, G. M. (2017). SgRNA Scorer 2.0: A Species-Independent Model to Predict CRISPR/Cas9 Activity. ACS Synthetic Biology, 6(5), 902-904. https://doi.org/10.1021/acssynbio.6b00343
    
    - Montague, T. G., Cruz, J. M., Gagnon, J. A., Church, G. M., & Valen, E. (2014). CHOPCHOP: A CRISPR/Cas9 and TALEN web tool for genome editing. Nucleic Acids Research, 42(W1), 401-407. https://doi.org/10.1093/nar/gku410
    
    - Whole-genome CRISPR target detection, implementing the method described in Sunagawa et al., "Mammalian Reverse Genetics without Crossing Reveals Nr3a as a Short-Sleeper Gene", Cell Reports 14(3)662-677, 2016

Input:
    - See config
    
    - A line-separated list of target sequences (input:exon-sequences in config)
    
    - A line-separated list of off-target sites (input:offtarget-sites in config)
    
Output:
    - A list of guides which were deemed efficient when accepted by at least two
    of: mm10db, CHOPCHOP-G20 and sgRNAScorer2.
    
    - A list of guides rejected due to being predicted as inefficient, or due to
    low specificity.

'''

import sys, argparse, os, re, joblib, timeit, ast, glob, string, random, json, time, glob
from sklearn.svm import SVC
from subprocess import call
from time import localtime, strftime, gmtime

from ConfigManager import ConfigManager
    
# flags indicate accept/reject, consensus count etc.
# order is important here, so we'll keep one list of flags (integers)
# and another for the corresponding flag labels
NUM_FLAGS = 10
FLAG_IDX = 0
FLAGS = []

# keep all the data relating to a target (eg: offtarget score, secondary structure formation, etc)
# the order is NOT important here, so we'll just construct a nested dictionary
# key: target sequence, value: a dictionary which is built on the fly
targetsData = {}

# binary encoding
encoding = {
    'A' : '0001',
    'C' : '0010',
    'T' : '0100',
    'G' : '1000',
    'K' : '1100',
    'M' : '0011',
    'R' : '1001',
    'Y' : '0110',
    'S' : '1010',
    'W' : '0101',
    'B' : '1110',
    'V' : '1011',
    'H' : '0111',
    'D' : '1101',
    'N' : '1111'
}

CODE_UNTESTED = "-"
CODE_ACCEPTED = 1
CODE_REJECT = 0
CODE_UNTESTED = '?'

# Function that returns the reverse-complement of a given sequence
def rc(dna):
    complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    rcseq = dna.translate(complements)[::-1]
    return rcseq

# Function that replaces U with T in the sequence (to go back from RNA to DNA)
def transToDNA(rna):
    switch_UT = str.maketrans('U', 'T')
    dna = rna.translate(switch_UT)
    return dna

# Function that calculates the AT% of a given sequence
def AT_percentage(seq):
    total = 0.0
    length = float(len(seq))
    for c in seq:
        if c in "AT":
            total += 1
    return 100.0*total/length

    
def printer(stringFormat):
    print('{}:\t{}\n'.format(
        strftime("%H:%M:%S", localtime()),
        stringFormat
    ))
    
#####################################
##  Begin analysis of single file  ##
##  or directory of files          ##
#####################################

def Crackling(configMngr):
    totalSizeBytes = configMngr.getDatasetSizeBytes()
    completedSizeBytes = 0

    sys.stdout = configMngr.getLogMethod()

    seqFileCounter = 0   
    lastRunTimeSec = 0
    lastScaffoldSizeBytes = 0
    totalRunTimeSec = 0
    
    printer('Crackling is starting...')
    
    for seqFilePath in configMngr.getIterFilesToProcess():
        ###################################
        ##        Begin this run         ##
        ###################################
        printer('*' * 100)

        FLAG_IDX = 0
        FLAGS = []

        sys.stdout.log.flush()

        start_time = time.time()

        PROCESS_CODE = configMngr.getConfigName()
        
        printer('{} of {} bytes processed ({}%)'.format(
            completedSizeBytes,
            totalSizeBytes,
            round((float(completedSizeBytes) / float(totalSizeBytes) * 100.0), 3)
        ))

        lastScaffoldSizeBytes = os.path.getsize(seqFilePath)

        #if os.path.isfile(configMngr['output']['file']):
        #    printer('{} already processed, skipping...'.format(PROCESS_CODE))
        #    seqFileCounter += 1
        #    completedSizeBytes += lastScaffoldSizeBytes
        #    print('-------------------------------------\n'*5)
        #    continue

        completedSizeBytes += lastScaffoldSizeBytes
        
        seqFileCounter += 1



        ###################################
        ##   Processing the input file   ##
        ###################################
        FLAGS.append('seqCount')

        printer('Identifying possible target sites in: {}'.format(
            seqFilePath
        ))

        possibleTargets = {}

        pattern_forward = r"(?=([ATCG]{21}GG))"
        pattern_reverse = r"(?=(CC[ACGT]{21}))"

        with open(seqFilePath, 'r') as inFile:
            for line_seq in inFile:
                if line_seq[0] == '>':
                    continue
                    
                # once for forward, once for reverse
                for pattern, seqModifier in [
                    [pattern_forward, lambda x : x], 
                    [pattern_reverse, lambda x : rc(x)]
                ]:
                    match_seq = re.findall(pattern, line_seq)
                    if match_seq:
                        for i in range(0, len(match_seq)):
                            target23 = seqModifier(match_seq[i])
                            if target23 in possibleTargets:
                                possibleTargets[target23][FLAG_IDX] += 1
                            else:
                                possibleTargets[target23] = [0] + [CODE_UNTESTED] * (NUM_FLAGS - 1)
                                targetsData[target23] = {}


        printer('Identified {} possible target sites.'.format(
            len(possibleTargets)
        ))

        ############################################
        ##   Removing targets that contain TTTT   ##
        ############################################
        FLAGS.append('passedTTTT')
        FLAG_IDX += 1

        if (configMngr['consensus'].getboolean('mm10db')):
            printer('Removing all targets that contain TTTT.')

            failedCount = 0
            for target23 in possibleTargets:
                if "TTTT" in target23:
                    possibleTargets[target23][FLAG_IDX] = CODE_REJECT # reject due to this reason
                    failedCount += 1
                else:
                    possibleTargets[target23][FLAG_IDX] = CODE_ACCEPTED # accept due to this reason

            printer('\t{} of {} failed here.'.format(
                failedCount,
                len(possibleTargets)
            ))

                
        #########################################
        ##    AT% ideally is between 20-65%    ##
        ######################################### 
        FLAGS.append('passedAT2065')
        FLAG_IDX += 1

        if (configMngr['consensus'].getboolean('mm10db')):
            printer('Calculating AT 20-65.')

            failedCount = 0
            for target23 in possibleTargets:
                AT = AT_percentage(target23[0:20])
                if AT < 20 or AT > 65:
                    possibleTargets[target23][FLAG_IDX] = CODE_REJECT # reject due to this reason
                    failedCount += 1
                else:
                    possibleTargets[target23][FLAG_IDX] = CODE_ACCEPTED # accept due to this reason
                
                targetsData[target23]['AT'] = AT

            printer('\t{} of {} failed here.'.format(
                failedCount,
                len(possibleTargets)
            ))
           
        #########################################
        ##                 G20                 ##
        #########################################
        FLAGS.append('passedG20')
        FLAG_IDX += 1

        if (configMngr['consensus'].getboolean('CHOPCHOP')):
            printer('Check for G20.')

            failedCount = 0
            for target23 in possibleTargets:
                if target23[19] != 'G':
                    possibleTargets[target23][FLAG_IDX] = CODE_REJECT # reject due to this reason
                    failedCount += 1
                else:
                    possibleTargets[target23][FLAG_IDX] = CODE_ACCEPTED # accept due to this reason

            printer('\t{} of {} failed here.'.format(
                failedCount,
                len(possibleTargets)
            ))      
                
                
        ##########################################
        ##   Calculating secondary structures   ##
        ##########################################
        # Do not check SS if any other of mm10db's checks have already failed   

        FLAGS.append('passedSecondaryStructure')
        FLAG_IDX += 1

        if (configMngr['consensus'].getboolean('mm10db')):

            printer('Secondary structure analysis.')

            guide = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU"
            pattern_RNAstructure = r".{28}\({4}\.{4}\){4}\.{3}\){4}.{21}\({4}\.{4}\){4}\({7}\.{3}\){7}\.{3}\s\((.+)\)"
            pattern_RNAenergy = r"\s\((.+)\)"

            call(
                ["rm -f \"{}\"".format(configMngr['rnafold']['output'])], 
                shell=True
            )

            printer('\tConstructing the RNAfold input file.')


            testedCount = 0
            with open(configMngr['rnafold']['input'], 'w+') as fRnaInput:
                for target23 in possibleTargets:
                
                    
                    # we don't want guides on + starting with T, or on - ending with A
                    # and only things that have passed everything else so far
                    if  not ( 
                            (target23[-2:] == 'GG' and target23[0] == 'T') or 
                            (target23[:2] == 'CC' and target23[-1] == 'A')
                        ) and \
                        possibleTargets[target23][FLAGS.index('passedAT2065')] == 1 and \
                        possibleTargets[target23][FLAGS.index('passedTTTT')] == 1:
                        fRnaInput.write(
                            "G{}{}\n".format(
                                target23[1:20], 
                                guide
                            )
                        )
                        
                        testedCount += 1
                        
            printer('\tFile ready. {} to calculate.'.format(
                testedCount
            ))

            printer('\n\n=========== (Start RNAfold) =========')

            call(
                "{} --noPS -j{} -i \"{}\" >> \"{}\"".format(
                    configMngr['rnafold']['binary'],
                    configMngr['rnafold']['threads'],
                    configMngr['rnafold']['input'],
                    configMngr['rnafold']['output']
                ), 
                shell=True
            )

            printer('\n\n============ (End RNAfold) ==========')

            total_number_structures = len(possibleTargets)

            printer('\tStarting to process the RNAfold results.')

            RNA_structures = None
            with open(configMngr['rnafold']['output'], 'r') as fRnaOutput:
                RNA_structures = fRnaOutput.readlines()

            i=0
            failedCount = 0
            for target23 in possibleTargets:
                targetsData[target23]['L1'] = None
                targetsData[target23]['structure'] = None
                targetsData[target23]['energy'] = None

                # we don't want guides on + starting with T, or on - ending with A
                # and only things that have passed everything else so far
                if  not ( 
                        (target23[-2:] == 'GG' and target23[0] == 'T') or 
                        (target23[:2] == 'CC' and target23[-1] == 'A')
                    ) and possibleTargets[target23][FLAGS.index('passedAT2065')] == 1 and possibleTargets[target23][FLAGS.index('passedTTTT')] == 1:
                
                    
                    L1 = RNA_structures[2*i].rstrip()
                    L2 = RNA_structures[2*i+1].rstrip()
                    
                    structure = L2.split(" ")[0]
                    energy = L2.split(" ")[1][1:-1]
                    
                    targetsData[target23]['L1'] = L1
                    targetsData[target23]['structure'] = structure
                    targetsData[target23]['energy'] = energy
                    
                    target = L1[:20]
                    if transToDNA(target) != target23[0:20] and transToDNA("C"+target[1:]) != target23[0:20] and transToDNA("A"+target[1:]) != target23[0:20]:
                        print("Error? "+target23+"\t"+target)
                        quit()

                    match_structure = re.search(pattern_RNAstructure, L2)
                    if match_structure:
                        energy = ast.literal_eval(match_structure.group(1))
                        if energy < float(configMngr['rnafold']['low_energy_threshold']):
                            possibleTargets[transToDNA(target23)][FLAG_IDX] = CODE_REJECT # reject due to this reason
                            failedCount += 1
                        else:
                            possibleTargets[target23][FLAG_IDX] = CODE_ACCEPTED # accept due to this reason
                    else:
                        match_energy = re.search(pattern_RNAenergy, L2)
                        if match_energy:
                            energy = ast.literal_eval(match_energy.group(1))
                            if energy <= float(configMngr['rnafold']['high_energy_threshold']):
                                possibleTargets[transToDNA(target23)][FLAG_IDX] = CODE_REJECT # reject due to this reason
                                failedCount += 1
                            else:
                                possibleTargets[target23][FLAG_IDX] = CODE_ACCEPTED # accept due to this reason
                    i+=1

            printer('\t{} of {} failed here.'.format(
                failedCount,
                testedCount
            ))
        
        #########################################
        ##         Calc mm10db result          ##
        ######################################### 
        FLAGS.append('acceptedByMm10db')
        FLAG_IDX += 1

        if (configMngr['consensus'].getboolean('mm10db')):
            printer('Calculating mm10db result.')

            acceptedCount = 0
            failedCount = 0
            
            for target23 in possibleTargets:
                if (
                    # reject if: failed any tests or guide is on + starting with T or guide is on - ending with A
                    possibleTargets[target23][FLAGS.index('passedAT2065')] == 0               or     # accepted when tested on AT2065
                    possibleTargets[target23][FLAGS.index('passedTTTT')] == 0                 or     # accepted when tested on TTTT
                    possibleTargets[target23][FLAGS.index('passedSecondaryStructure')] == 0   or
                    ((target23[-2:] == 'GG' and target23[0] == 'T') or (target23[:2] == 'CC' and target23[-1] == 'A'))        # accepted when tested on SS
                ):
                    # mm10db rejected the guide
                    possibleTargets[target23][FLAG_IDX] = CODE_REJECT # reject due to this reason
                    failedCount += 1
                else:
                    possibleTargets[target23][FLAG_IDX] = CODE_ACCEPTED # accept due to this reason
                    acceptedCount += 1

            printer('\t{} accepted.'.format(
                acceptedCount
            ))

            printer('\t{} failed.'.format(
                failedCount
            ))
                
            del acceptedCount
                 
        #########################################
        ##         sgRNAScorer 2.0 model       ##
        ######################################### 
        FLAGS.append('acceptedBySgRnaScorer')
        FLAG_IDX += 1

        if (configMngr['consensus'].getboolean('sgRNAScorer2')):
            printer('Calculating sgRNAScorer2 scores.')

            clfLinear = joblib.load(configMngr['sgrnascorer2']['model'])

            failedCount = 0
            testedCount = 0
            for target23 in possibleTargets:
                # seeing that this is the last efficacy test...
                # only test if the guide hasn't passed the consensus approach yet AND has the potential to pass
                # eg: for a guide to be deemed efficient, it must pass 2 of 3 tests
                # if it has passed none so far then don't bother using the sgRNAscorer2 model
                # if it is only one test away from passing (threshold - 1) then use the sgRNAScorer2 model
                # if it has already passed the threshold then don't bother
                
                # how many tools were potentially ran?
                consensusToolsCount = configMngr.getNumberToolsInConsensus()
                
                currentConsensus = ((int)(possibleTargets[target23][FLAGS.index('acceptedByMm10db')] == 1) +    # mm10db accepted
                    (int)(possibleTargets[target23][FLAGS.index('passedG20')] == 1))                            # chopchop-g20 accepted
                if (currentConsensus == (int(configMngr['consensus']['n']) - 1)):
                    sequence = target23.upper()
                    entryList = []
                    testedCount += 1
                    
                    x = 0
                    while x < 20:
                        y = 0
                        while y < 4:
                            entryList.append(int(encoding[sequence[x]][y]))
                            y += 1
                        x += 1
                        
                    # predict based on the entry
                    prediction = clfLinear.predict([entryList])
                    score = clfLinear.decision_function([entryList])[0]

                    targetsData[target23]['sgrnascorer2score'] = score

                    if float(score) < float(float(configMngr['sgrnascorer2']['score-threshold'])):
                        possibleTargets[target23][FLAG_IDX] = CODE_REJECT # reject due to this reason
                        failedCount += 1
                    else:
                        possibleTargets[target23][FLAG_IDX] = CODE_ACCEPTED # accept due to this reason
                        
                else:
                    targetsData[target23]['sgrnascorer2score'] = None
              
            printer('\t{} of {} failed here.'.format(
                failedCount,
                testedCount
            ))

            
        #########################################
        ##      Begin efficacy consensus       ##
        ######################################### 
        FLAGS.append('consensusCount')
        FLAG_IDX += 1

        printer('Beginning efficacy consensus.')

        for target23 in possibleTargets:
            possibleTargets[target23][FLAG_IDX] = (
                (int)(possibleTargets[target23][FLAGS.index('acceptedByMm10db')] == 1)          +     # mm10db accepted
                (int)(possibleTargets[target23][FLAGS.index('acceptedBySgRnaScorer')] == 1)     +     # sgrnascorer2 accepted
                (int)(possibleTargets[target23][FLAGS.index('passedG20')] == 1)                       # chopchop-g20 accepted
            )

           
                
                
        ###############################################
        ##         Using Bowtie for positioning      ##
        ###############################################
        FLAGS.append('passedBowtie')
        FLAG_IDX += 1

        printer('Bowtie analysis.')

        printer('\tConstructing the Bowtie input file.')

        tempTargetDict_offset = {}
        testedCount = 0
        with open(configMngr['bowtie2']['input'], 'w') as fWriteBowtie:
            for target23 in possibleTargets:
                
                targetsData[target23]['chr'] = None
                targetsData[target23]['start'] = None
                targetsData[target23]['end'] = None
            
                if possibleTargets[target23][FLAGS.index('consensusCount')] >= int(configMngr['consensus']['n']):
                    testedCount += 1
                    similarTargets = [
                        target23[0:20] + "AGG", 
                        target23[0:20] + "CGG", 
                        target23[0:20] + "GGG", 
                        target23[0:20] + "TGG", 
                        target23[0:20] + "AAG", 
                        target23[0:20] + "CAG", 
                        target23[0:20] + "GAG", 
                        target23[0:20] + "TAG"
                    ]
                    
                    for seq in similarTargets:
                        fWriteBowtie.write(seq + "\n")
                        tempTargetDict_offset[seq] = target23
                
        printer('\tFile ready. Calling Bowtie.')

        printer('\n\n=========== (Start Bowtie2) =========')

        call("{} -x {} -p {} --reorder --no-hd -t -r -U \"{}\" -S \"{}\"".format(
           configMngr['bowtie2']['binary'],
           configMngr['input']['bowtie2-index'],
           configMngr['bowtie2']['threads'],
           configMngr['bowtie2']['input'],
           configMngr['bowtie2']['output'])
        ,shell=True)       

        printer('\n\n============ (End Bowtie2) ==========\n\n')

        printer('\tStarting to process the Bowtie results.')

        inFile = open(configMngr['bowtie2']['output'],'r')
        bowtieLines = inFile.readlines()
        inFile.close()

        i=0
        failedCount = 0
        while i<len(bowtieLines):
            nb_occurences = 0
            # we extract the read and use the dictionnary to find the corresponding target
            line = bowtieLines[i].rstrip().split("\t")
            chr = line[2]
            pos = ast.literal_eval(line[3])
            read = line[9]
            seq = ""
            
            if read in tempTargetDict_offset:
                seq = tempTargetDict_offset[read]
            elif rc(read) in tempTargetDict_offset:
                seq = tempTargetDict_offset[rc(read)]
            else:
                print("Problem? "+read)

            if seq[:-2] == 'GG':
                targetsData[seq]['chr'] = chr
                targetsData[seq]['start'] = pos
                targetsData[seq]['end'] = pos + 22
            elif rc(seq)[:2] == 'CC':
                targetsData[seq]['chr'] = chr
                targetsData[seq]['start'] = pos
                targetsData[seq]['end'] = pos + 22
            else:
                print("Error? "+seq)
                quit()  
                
            # we count how many of the eight reads for this target have a perfect alignment
            for j in range(i,i+8):
            
                # http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#sam-output
                # XM:i:<N>    The number of mismatches in the alignment. Only present if SAM record is for an aligned read.
                # XS:i:<N>    Alignment score for the best-scoring alignment found other than the alignment reported.
                
                if "XM:i:0" in bowtieLines[j]:
                    nb_occurences += 1
                    
                    # we also check whether this perfect alignment also happens elsewhere
                    if "XS:i:0"  in bowtieLines[j]:
                        nb_occurences += 1

            # if that number is at least two, the target is removed
            if nb_occurences > 1:
                possibleTargets[seq][FLAG_IDX] = CODE_REJECT # reject due to this reason
                failedCount += 1
            else:
                possibleTargets[seq][FLAG_IDX] = CODE_ACCEPTED # accept due to this reason
                
            # we continue with the next target
            i+=8

        # we can remove the dictionnary
        del tempTargetDict_offset

        printer('\t{} of {} failed here.'.format(
            failedCount,
            testedCount
        ))


        #########################################
        ##      Begin off-target scoring       ##
        #########################################   
        FLAGS.append('passedOffTargetScore')
        FLAG_IDX += 1

        printer('Beginning off-target scoring.')

        i=0
        targetsToRemove=[]

        # prepare the list of candidate guides to score
        testedCount = 0
        with open(configMngr['offtargetscore']['input'], 'w') as fTargetsToScore:
            for target23 in possibleTargets:
                targetsData[target23]['offtargetscore'] = None
                
                if possibleTargets[target23][FLAGS.index('consensusCount')] >= int(configMngr['consensus']['n']) and \
                    possibleTargets[target23][FLAGS.index('passedBowtie')] == 1:
                    target = target23[0:20]
                    fTargetsToScore.write(target+"\n")
                    testedCount += 1
         
        printer('\t{} to calculate scores.'.format(
            testedCount
        ))

        printer('\n\n============= (Start Offtarget) ============')

        # call the scoring method
        call(
            ["{} \"{}\" \"{}\" \"{}\" \"{}\" > \"{}\"".format(
                configMngr['offtargetscore']['binary'],
                configMngr['input']['offtarget-sites'],
                configMngr['offtargetscore']['input'],
                str(configMngr['offtargetscore']['max-distance']),
                str(configMngr['offtargetscore']['score-threshold']),
                configMngr['offtargetscore']['output'],
            )],
            shell = True
        )

        printer('\n\n============== (End Offtarget) =============\n\n')

        targetsScored = {}
        with open(configMngr['offtargetscore']['output'], 'r') as fTargetsScored:
            for targetScored in [x.split('\t') for x in fTargetsScored.readlines()]:
                if len(targetScored) == 2:
                    targetsScored[targetScored[0]] = float(targetScored[1].strip())

        failedCount = 0
        for target23 in possibleTargets:
            if target23[0:20] in targetsScored:
                score = targetsScored[target23[0:20]]
                targetsData[target23]['offtargetscore'] = score
                
                if score < float(configMngr['offtargetscore']['score-threshold']):
                    possibleTargets[target23][FLAG_IDX] = CODE_REJECT # reject due to this reason
                    failedCount += 1
                else:
                    possibleTargets[target23][FLAG_IDX] = CODE_ACCEPTED # accept due to this reason

        printer('\t{} of {} failed here.'.format(
            failedCount,
            testedCount
        ))
          
        #########################################
        ##           Begin output              ##
        #########################################   
        printer('Writing results to file.')

        # Write guides to file. Include scores etc.
        with open(configMngr['output']['file'], 'w+') as fOpen:
            
            fOpen.write('{}\n'.format(
                        configMngr['output']['delimiter'].join(
                            ['seq'] +
                            FLAGS +
                            [
                                'start',
                                'end',
                                'chr',
                                'offtargetscore',
                                'sgrnascorer2score',
                                'L1',
                                'structure',
                                'energy',
                            ]
                        )
                    )
                )
            
            for target23 in possibleTargets:
                output = [
                    target23
                ] + \
                    possibleTargets[target23] + \
                [
                    targetsData[target23].get('start', CODE_UNTESTED),
                    targetsData[target23].get('end', CODE_UNTESTED),
                    targetsData[target23].get('chr', CODE_UNTESTED),
                    targetsData[target23].get('offtargetscore', CODE_UNTESTED),
                    targetsData[target23].get('sgrnascorer2score', CODE_UNTESTED),
                    targetsData[target23].get('L1', CODE_UNTESTED),
                    targetsData[target23].get('structure', CODE_UNTESTED),
                    targetsData[target23].get('energy', CODE_UNTESTED),
                ]
                
                fOpen.write('{}\n'.format(
                        configMngr['output']['delimiter'].join(map(str, output))
                    )
                )


        #########################################
        ##              Clean up               ##
        #########################################
        printer('Clearning auxiliary files')
        for f in [
            configMngr['rnafold']['input'],
            configMngr['rnafold']['output'],
            configMngr['offtargetscore']['input'],
            configMngr['offtargetscore']['output'],
            configMngr['bowtie2']['input'],
            configMngr['bowtie2']['output'],
        ]:
            try:
                os.remove(f)
            except:
                pass
       


        #########################################
        ##               Done                  ##
        #########################################
        printer('Done.')

        printer('{} guides evaluated.'.format(
            len(possibleTargets)
        ))

        printer('Ran in {} (dd hh:mm:ss) or {} seconds'.format(
            strftime("%d %H:%M:%S", gmtime((time.time() - start_time))), 
            (time.time() - start_time)
        ))
        
        lastRunTimeSec = time.time() - start_time
        totalRunTimeSec += lastRunTimeSec
        
if __name__ == '__main__':
    # load in config
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', help='File to process', default=None, required=True)
    args = parser.parse_args()

    configMngr = ConfigManager(args.c, lambda x : print(f'configMngr says: {x}'))

    if not configMngr.isConfigured():
        print('Something went wrong with reading the configMngruration.')
        exit()
    else:
        Crackling(configMngr)
        
    print('Goodbye.')