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

import sys, argparse, os, re, joblib, timeit, ast, glob, string, random, json, time, glob, csv
from sklearn.svm import SVC
from subprocess import call
from time import localtime, strftime, gmtime

from ConfigManager import ConfigManager

# binary encoding
encoding = {
    'A' : '0001',    'C' : '0010',    'T' : '0100',    'G' : '1000',
    'K' : '1100',    'M' : '0011',    'R' : '1001',    'Y' : '0110',
    'S' : '1010',    'W' : '0101',    'B' : '1110',    'V' : '1011',
    'H' : '0111',    'D' : '1101',    'N' : '1111'
}

CODE_ACCEPTED = 1
CODE_REJECTED = 0
CODE_UNTESTED = "?"
CODE_AMBIGUOUS = "-"
CODE_ERROR = "!"

DEFAULT_GUIDE_PROPERTIES = {
    'seq'                       : "",
    'header'                    : "",
    'seqCount'                  : 1,
    'start'                     : CODE_UNTESTED,
    'end'                       : CODE_UNTESTED,
    'strand'                    : CODE_UNTESTED,
    'passedTTTT'                : CODE_UNTESTED,
    'passedATPercent'           : CODE_UNTESTED,
    'passedG20'                 : CODE_UNTESTED,
    'passedSecondaryStructure'  : CODE_UNTESTED,
    'ssL1'                      : CODE_UNTESTED,
    'ssStructure'               : CODE_UNTESTED,
    'ssEnergy'                  : CODE_UNTESTED,
    'acceptedByMm10db'          : CODE_UNTESTED,
    'acceptedBySgRnaScorer'     : CODE_UNTESTED,
    'consensusCount'            : CODE_UNTESTED,
    'passedBowtie'              : CODE_UNTESTED,
    'passedOffTargetScore'      : CODE_UNTESTED,
    'sgrnascorer2score'         : CODE_UNTESTED,
    'AT'                        : CODE_UNTESTED,
    'bowtieChr'                 : CODE_UNTESTED,
    'bowtieStart'               : CODE_UNTESTED,
    'bowtieEnd'                 : CODE_UNTESTED,
    'offtargetscore'            : CODE_UNTESTED,
    'passedAvoidLeadingT'       : CODE_UNTESTED,
    #'passedReversePrimer'      : CODE_UNTESTED,
}

DEFAULT_GUIDE_PROPERTIES_ORDER = [
    'seq',
    'sgrnascorer2score',
    'header',
    'start',
    'end',
    'strand',
    'seqCount',
    'passedG20',
    'passedTTTT',
    'passedATPercent',
    'passedSecondaryStructure',
    'ssL1',
    'ssStructure',
    'ssEnergy',
    'acceptedByMm10db',
    'acceptedBySgRnaScorer',
    'consensusCount',
    'passedBowtie',
    'passedOffTargetScore',
    'AT',
    'bowtieChr',
    'bowtieStart',
    'bowtieEnd',
    'offtargetscore',
    'passedAvoidLeadingT',
    #'passedReversePrimer',
]

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

# Function that aligns two given sequences and returns their similarity
def NeedlemanWunsch(seq1,seq2):
    d = -5 # Gap penalty
    A = seq1 # First sequence to be compared
    B = seq2 # Second ""
    I = range(len(seq1)) # To help iterate (Pythonic)
    J = range(len(seq2)) # ""
    F = [[0 for i in seq1] for j in seq2] # Fill a 2D array with zeroes
    # Similarity matrix from Wikipedia:
    S = \
    {'A': {'A': 10, 'G': -1, 'C': -3, 'T': -4},
        'G': {'A': -1, 'G':  7, 'C': -5, 'T': -3},
        'C': {'A': -3, 'G': -5, 'C':  9, 'T':  0},
        'T': {'A': -4, 'G': -3, 'C':  0, 'T':  8}}

    # Initialization
    for i in I:
        F[i][0] = d * i
    for j in J:
        F[0][j] = d * j

    # Scoring
    for i in I[1:]:
        for j in J[1:]:
            Match = F[i-1][j-1] + S[A[i]][B[j]]
            Delete = F[i-1][j] + d
            Insert = F[i][j-1] + d
            F[i][j] = max(Match, Insert, Delete)

    # Traceback
    AlignmentA = ""
    AlignmentB = ""
    i = len(seq1) - 1
    j = len(seq2) - 1

    while (i > 0 and j > 0):
        Score = F[i][j]
        ScoreDiag = F[i - 1][j - 1]
        ScoreUp = F[i][j - 1]
        ScoreLeft = F[i - 1][j]
        if (Score == ScoreDiag + S[A[i]][B[j]]):
            AlignmentA = A[i] + AlignmentA
            AlignmentB = B[j] + AlignmentB
            i -= 1
            j -= 1
        elif (Score == ScoreLeft + d):
            AlignmentA = A[i] + AlignmentA
            AlignmentB = "-" + AlignmentB
            i -= 1
        elif (Score == ScoreUp + d):
            AlignmentA = "-" + AlignmentA
            AlignmentB = B[j] + AlignmentB
            j -= 1
        else:
            print("algorithm error?")
    while (i > 0):
        AlignmentA = A[i] + AlignmentA
        AlignmentB = "-" + AlignmentB
        i -= 1
    while (j > 0):
        AlignmentA = "-" + AlignmentA
        AlignmentB = B[j] + AlignmentB
        j -= 1

    # Similarity
    lenA = len(AlignmentA)
    lenB = len(AlignmentB)
    sim1 = ""
    sim2 = ""
    len0 = 0
    k = 0
    total = 0.0
    similarity = 0.0

    if (lenA > lenB):
        sim1 = AlignmentA
        sim2 = AlignmentB
        len0 = lenA
    else:
        sim1 = AlignmentB
        sim2 = AlignmentA
        len0 = lenB

    while (k < len0):
        if (sim1[k] == sim2[k]):
            total += 1
        k += 1

    similarity = total / len0 * 100

    return similarity
    
def printer(stringFormat):
    print('>>> {}:\t{}\n'.format(
        strftime("%Y-%m-%d %H:%M:%S", localtime()),
        stringFormat
    ))

def caller(*args, **kwargs):
    printer(f"| Calling: {args}")
    call(*args, **kwargs)
    printer(f"| Finished")
 
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

        sys.stdout.log.flush()

        start_time = time.time()

        PROCESS_CODE = configMngr.getConfigName()
        
        printer('{} of {} bytes processed ({}%)'.format(
            completedSizeBytes,
            totalSizeBytes,
            round((float(completedSizeBytes) / float(totalSizeBytes) * 100.0), 3)
        ))

        lastScaffoldSizeBytes = os.path.getsize(seqFilePath)

        completedSizeBytes += lastScaffoldSizeBytes
        
        seqFileCounter += 1

        ###################################
        ##   Processing the input file   ##
        ###################################
        printer('Identifying possible target sites in: {}'.format(
            seqFilePath
        ))
        
        # key: FASTA header, value: sequence
        seqsInFile = {}

        with open(seqFilePath, 'r') as inFile:
            previousHeader = seqFilePath
            for line in inFile:
                line = line.strip()
                if line[0] == '>':
                    seqsInFile[line[1:]] = ""
                    previousHeader = line[1:]
                    continue
                else:
                    if previousHeader not in seqsInFile:
                        # it could be a plain text file, without a header
                        seqsInFile[previousHeader] = ""
                        
                    seqsInFile[previousHeader] += line.strip()

        printer(f"The file contained {len(seqsInFile)} sequence(s). The length of the sequence(s): {', '.join([str(len(seqsInFile[x])) for x in seqsInFile])}")

        possibleTargets = {}

        pattern_forward = r"(?=([ATCG]{21}GG))"
        pattern_reverse = r"(?=(CC[ACGT]{21}))"

        for seqHeader in seqsInFile:
            seq = seqsInFile[seqHeader]
            
            # once for forward, once for reverse
            for pattern, strand, seqModifier in [
                [pattern_forward, '+', lambda x : x], 
                [pattern_reverse, '-', lambda x : rc(x)]
            ]:
                p = re.compile(pattern)
                for m in p.finditer(seq):
                    target23 = seqModifier(seq[m.start() : m.start() + 23])
                    if target23 not in possibleTargets:
                        possibleTargets[target23] = DEFAULT_GUIDE_PROPERTIES.copy()
                        possibleTargets[target23]['seq'] = target23
                        possibleTargets[target23]['header'] = seqHeader
                        possibleTargets[target23]['start'] = m.start()
                        possibleTargets[target23]['end'] = m.start() + 23
                        possibleTargets[target23]['strand'] = strand
                    else:
                        # we've already seen this guide, make the positioning ambiguous
                        possibleTargets[target23]['seqCount'] += 1
                        possibleTargets[target23]['header'] = CODE_AMBIGUOUS
                        possibleTargets[target23]['start'] = CODE_AMBIGUOUS
                        possibleTargets[target23]['end'] = CODE_AMBIGUOUS
                        possibleTargets[target23]['strand'] = CODE_AMBIGUOUS


        printer('Identified {} possible target sites.'.format(
            len(possibleTargets)
        ))

        ############################################
        ##     Removing targets with leading T    ##
        ############################################
        if (configMngr['consensus'].getboolean('mm10db')):
            printer('mm10db - remove all targets with a leading T (+) or trailing A (-).')

            failedCount = 0
            for target23 in possibleTargets:
                if (target23[-2:] == 'GG' and target23[0] == 'T') or \
                    (target23[:2] == 'CC' and target23[-1] == 'A'):
                    possibleTargets[target23]['passedAvoidLeadingT'] = CODE_REJECTED # reject due to this reason
                    failedCount += 1
                else:
                    possibleTargets[target23]['passedAvoidLeadingT'] = CODE_ACCEPTED # accept due to this reason

            printer('\t{} of {} failed here.'.format(
                failedCount,
                len(possibleTargets)
            ))
            
        #########################################
        ##    AT% ideally is between 20-65%    ##
        ######################################### 
        if (configMngr['consensus'].getboolean('mm10db')):
            printer('mm10db - remove based on AT percent.')

            failedCount = 0
            for target23 in possibleTargets:
                AT = AT_percentage(target23[0:20])
                #if AT < 45:
                if AT < 20 or AT > 65:
                    possibleTargets[target23]['passedATPercent'] = CODE_REJECTED # reject due to this reason
                    failedCount += 1
                else:
                    possibleTargets[target23]['passedATPercent'] = CODE_ACCEPTED # accept due to this reason
                
                possibleTargets[target23]['AT'] = AT

            printer('\t{} of {} failed here.'.format(
                failedCount,
                len(possibleTargets)
            ))
            
            #printer(f"\tmm10db equiv: {sum([(possibleTargets[target23]['passedAvoidLeadingT'] and possibleTargets[target23]['passedATPercent']) for target23 in possibleTargets])}")

        ############################################
        ##   Removing targets that contain TTTT   ##
        ############################################
        if (configMngr['consensus'].getboolean('mm10db')):
            printer('mm10db - remove all targets that contain TTTT.')

            failedCount = 0
            for target23 in possibleTargets:
                if "TTTT" in target23:
                    possibleTargets[target23]['passedTTTT'] = CODE_REJECTED # reject due to this reason
                    failedCount += 1
                else:
                    possibleTargets[target23]['passedTTTT'] = CODE_ACCEPTED # accept due to this reason

            printer('\t{} of {} failed here.'.format(
                failedCount,
                len(possibleTargets)
            ))
            
            #printer(f"\tmm10db equiv: {sum([(possibleTargets[target23]['passedAvoidLeadingT'] and possibleTargets[target23]['passedATPercent'] and possibleTargets[target23]['passedTTTT']) for target23 in possibleTargets])}")
            
        # #########################################
        # ##     Distance to reverse primer      ##
        # #########################################
        # if (configMngr['consensus'].getboolean('mm10db')):
        #     printer('Check if guide is too close to reverse primer.')
        # 
        #     failedCount = 0
        #     for target23 in possibleTargets:
        #         if NeedlemanWunsch(target23[0:20], "AAAAGCACCGACTCGGTGCC") > 60:
        #             possibleTargets[target23]['passedReversePrimer'] = CODE_REJECTED # reject due to this reason
        #             failedCount += 1
        #         else:
        #             possibleTargets[target23]['passedReversePrimer'] = CODE_ACCEPTED # accept due to this reason
        # 
        #     printer('\t{} of {} failed here.'.format(
        #         failedCount,
        #         len(possibleTargets)
        #     ))      
                
                
        ##########################################
        ##   Calculating secondary structures   ##
        ##########################################
        # Do not check SS if any other of mm10db's checks have already failed   
        if (configMngr['consensus'].getboolean('mm10db')):
            printer('mm10db - check secondary structure.')

            guide = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU"
            pattern_RNAstructure = r".{28}\({4}\.{4}\){4}\.{3}\){4}.{21}\({4}\.{4}\){4}\({7}\.{3}\){7}\.{3}\s\((.+)\)"
            pattern_RNAenergy = r"\s\((.+)\)"

            caller(
                ["rm -f \"{}\"".format(configMngr['rnafold']['output'])], 
                shell=True
            )

            printer('\tConstructing the RNAfold input file.')


            testedCount = 0
            with open(configMngr['rnafold']['input'], 'w+') as fRnaInput:
                for target23 in possibleTargets:
                    # we only guides that have passed for mm10db's test so far
                    if  possibleTargets[target23]['passedAvoidLeadingT'] == 1 and \
                        possibleTargets[target23]['passedATPercent'] == 1 and \
                        possibleTargets[target23]['passedTTTT'] == 1:
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

            caller(
                "{} --noPS -j{} -i \"{}\" > \"{}\"".format(
                    configMngr['rnafold']['binary'],
                    configMngr['rnafold']['threads'],
                    configMngr['rnafold']['input'],
                    configMngr['rnafold']['output']
                ), 
                shell=True
            )

            printer('\tStarting to process the RNAfold results.')

            RNAstructures = {}
            with open(configMngr['rnafold']['output'], 'r') as fRnaOutput:
                i = 0
                L1, L2, target = None, None, None
                for line in fRnaOutput:
                    if i % 2 == 0:
                        # 0th, 2nd, 4th, etc.
                        L1 = line.rstrip()
                        target = L1[0:20]
                    else:
                        # 1st, 3rd, 5th, etc.
                        L2 = line.rstrip()
                        RNAstructures[transToDNA(target[1:20])] = [
                            L1, L2, target
                        ]

                    i += 1

            failedCount = 0
            errorCount = 0
            notFoundCount = 0 
            for target23 in possibleTargets:
                # we only guides that have passed for mm10db's test so far
                if  possibleTargets[target23]['passedAvoidLeadingT'] == 1 and \
                    possibleTargets[target23]['passedATPercent'] == 1 and \
                    possibleTargets[target23]['passedTTTT'] == 1:
                    
                    key = target23[1:20]
                    if key not in RNAstructures:
                        print(f"Could not find: {target23[0:20]}")
                        notFoundCount += 1
                        continue
                    else:
                        L1 = RNAstructures[key][0]
                        L2 = RNAstructures[key][1]
                        target = RNAstructures[key][2]

                    structure = L2.split(" ")[0]
                    energy = L2.split(" ")[1][1:-1]
                    
                    possibleTargets[target23]['ssL1'] = L1
                    possibleTargets[target23]['ssStructure'] = structure
                    possibleTargets[target23]['ssEnergy'] = energy
                    
                    if transToDNA(target) != target23[0:20] and transToDNA("C"+target[1:]) != target23[0:20] and transToDNA("A"+target[1:]) != target23[0:20]:
                        possibleTargets[target23]['passedSecondaryStructure'] = CODE_ERROR
                        #print(
                        #    f"Error?",
                        #    target23,
                        #    target
                        #)
                        errorCount += 1
                        continue

                    match_structure = re.search(pattern_RNAstructure, L2)
                    if match_structure:
                        energy = ast.literal_eval(match_structure.group(1))
                        if energy < float(configMngr['rnafold']['low_energy_threshold']):
                            possibleTargets[transToDNA(target23)]['passedSecondaryStructure'] = CODE_REJECTED # reject due to this reason
                            failedCount += 1
                        else:
                            possibleTargets[target23]['passedSecondaryStructure'] = CODE_ACCEPTED # accept due to this reason
                    else:
                        match_energy = re.search(pattern_RNAenergy, L2)
                        if match_energy:
                            energy = ast.literal_eval(match_energy.group(1))
                            if energy <= float(configMngr['rnafold']['high_energy_threshold']):
                                possibleTargets[transToDNA(target23)]['passedSecondaryStructure'] = CODE_REJECTED # reject due to this reason
                                failedCount += 1
                            else:
                                possibleTargets[target23]['passedSecondaryStructure'] = CODE_ACCEPTED # accept due to this reason


            printer('\t{} of {} failed here.'.format(
                failedCount,
                testedCount
            ))
            
            if errorCount > 0:
                printer('\t{} of {} erred here.'.format(
                    errorCount,
                    testedCount
                ))
            
            if notFoundCount > 0:
                printer('\t{} of {} not found in RNAfold output.'.format(
                    notFoundCount,
                    testedCount
                ))

            #printer(f"\tmm10db equiv: {sum([(possibleTargets[target23]['passedAvoidLeadingT'] and possibleTargets[target23]['passedATPercent'] and possibleTargets[target23]['passedTTTT'] and possibleTargets[target23]['passedSecondaryStructure']) for target23 in possibleTargets])}")
            
        #########################################
        ##         Calc mm10db result          ##
        ######################################### 
        if (configMngr['consensus'].getboolean('mm10db')):
            printer('Calculating mm10db final result.')

            acceptedCount = 0
            failedCount = 0
            
            for target23 in possibleTargets:
                if (
                    # reject if: failed any tests or guide is on + starting with T or guide is on - ending with A
                    possibleTargets[target23]['passedATPercent'] != CODE_ACCEPTED            or     # rejected when tested on ATPercent
                    possibleTargets[target23]['passedTTTT'] != CODE_ACCEPTED                 or     # rejected when tested on TTTT
                    possibleTargets[target23]['passedSecondaryStructure'] != CODE_ACCEPTED   or     # rejected due to secondary structure
                    possibleTargets[target23]['passedAvoidLeadingT'] != CODE_ACCEPTED               # rejected as it started with a T (+) or ended with an A (-)
                ):
                    # mm10db rejected the guide
                    possibleTargets[target23]['acceptedByMm10db'] = CODE_REJECTED # reject due to this reason
                    failedCount += 1
                else:
                    possibleTargets[target23]['acceptedByMm10db'] = CODE_ACCEPTED # accept due to this reason
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
        if (configMngr['consensus'].getboolean('sgRNAScorer2')):
            printer('sgRNAScorer2 - score using model.')

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
                
                currentConsensus = ((int)(possibleTargets[target23]['acceptedByMm10db'] == CODE_ACCEPTED) +    # mm10db accepted
                    (int)(possibleTargets[target23]['passedG20'] == CODE_ACCEPTED))                            # chopchop-g20 accepted

                # 5991
                #if possibleTargets[target23]['seqCount'] != 1:# or possibleTargets[target23]['start'] > :
                    #continue
                
                #if (currentConsensus == (int(configMngr['consensus']['n']) - 1)):
                if True:
                    sequence = target23.upper()
                    entryList = []
                    testedCount += 1

                    for x in range(0, 20):
                        for y in range(0, 4):
                            entryList.append(int(encoding[sequence[x]][y]))

                    # predict based on the entry
                    prediction = clfLinear.predict([entryList])
                    score = clfLinear.decision_function([entryList])[0] 

                    possibleTargets[target23]['sgrnascorer2score'] = score

                    if float(score) < float(float(configMngr['sgrnascorer2']['score-threshold'])):
                        possibleTargets[target23]['acceptedBySgRnaScorer'] = CODE_REJECTED # reject due to this reason
                        failedCount += 1
                    else:
                        possibleTargets[target23]['acceptedBySgRnaScorer'] = CODE_ACCEPTED # accept due to this reason
                            
            printer('\t{} of {} failed here.'.format(
                failedCount,
                testedCount
            ))
           
        #########################################
        ##                 G20                 ##
        #########################################
        if (configMngr['consensus'].getboolean('CHOPCHOP')):
            printer('CHOPCHOP - remove those without G in position 20.')

            failedCount = 0
            for target23 in possibleTargets:
                if target23[19] != 'G':
                    possibleTargets[target23]['passedG20'] = CODE_REJECTED # reject due to this reason
                    failedCount += 1
                else:
                    possibleTargets[target23]['passedG20'] = CODE_ACCEPTED # accept due to this reason

            printer('\t{} of {} failed here.'.format(
                failedCount,
                len(possibleTargets)
            ))     

        #########################################
        ##      Begin efficacy consensus       ##
        ######################################### 
        printer('Beginning efficacy consensus.')
        
        failedCount = 0
        for target23 in possibleTargets:
            possibleTargets[target23]['consensusCount'] = (
                (int)(possibleTargets[target23]['acceptedByMm10db'] == CODE_ACCEPTED)          +     # mm10db accepted
                (int)(possibleTargets[target23]['acceptedBySgRnaScorer'] == CODE_ACCEPTED)     +     # sgrnascorer2 accepted
                (int)(possibleTargets[target23]['passedG20'] == CODE_ACCEPTED)                       # chopchop-g20 accepted
            )
            if possibleTargets[target23]['consensusCount'] < int(configMngr['consensus']['n']):
                failedCount += 1
                
        printer('\t{} of {} failed here.'.format(
            failedCount,
            len(possibleTargets)
        ))     

        if (configMngr['offtargetscore'].getboolean('enabled')):
            #### BEGIN: REMOVE OFF-TARGET SCORING
            
            ###############################################
            ##         Using Bowtie for positioning      ##
            ###############################################
            printer('Bowtie analysis.')
            
            printer('\tConstructing the Bowtie input file.')
            
            tempTargetDict_offset = {}
            testedCount = 0
            with open(configMngr['bowtie2']['input'], 'w') as fWriteBowtie:
                for target23 in possibleTargets:
                
                    if possibleTargets[target23]['consensusCount'] >= int(configMngr['consensus']['n']):
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
            
            caller("{} -x {} -p {} --reorder --no-hd -t -r -U \"{}\" -S \"{}\"".format(
            configMngr['bowtie2']['binary'],
            configMngr['input']['bowtie2-index'],
            configMngr['bowtie2']['threads'],
            configMngr['bowtie2']['input'],
            configMngr['bowtie2']['output'])
            ,shell=True)       
            
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
                    possibleTargets[seq]['bowtieChr'] = chr
                    possibleTargets[seq]['bowtieStart'] = pos
                    possibleTargets[seq]['bowtieEnd'] = pos + 22
                elif rc(seq)[:2] == 'CC':
                    possibleTargets[seq]['bowtieChr'] = chr
                    possibleTargets[seq]['bowtieStart'] = pos
                    possibleTargets[seq]['bowtieEnd'] = pos + 22
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
                    possibleTargets[seq]['passedBowtie'] = CODE_REJECTED # reject due to this reason
                    failedCount += 1
                else:
                    possibleTargets[seq]['passedBowtie'] = CODE_ACCEPTED # accept due to this reason
                    
                # we continue with the next target
                i+=8
            
            # we can remove the dictionary
            del tempTargetDict_offset
            
            printer('\t{} of {} failed here.'.format(
                failedCount,
                testedCount
            ))
            
            
            #########################################
            ##      Begin off-target scoring       ##
            #########################################   
            printer('Beginning off-target scoring.')
            
            i=0
            targetsToRemove=[]
            
            # prepare the list of candidate guides to score
            testedCount = 0
            
            with open(configMngr['offtargetscore']['input'], 'w') as fTargetsToScore:
                for target23 in possibleTargets:
                    
                    if possibleTargets[target23]['consensusCount'] >= int(configMngr['consensus']['n']) and \
                        possibleTargets[target23]['passedBowtie'] == 1:
                        target = target23[0:20]
                        fTargetsToScore.write(target+"\n")
                        testedCount += 1
            
            printer('\t{} to calculate scores.'.format(
                testedCount
            ))
            
            # call the scoring method
            caller(
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
            
            targetsScored = {}
            with open(configMngr['offtargetscore']['output'], 'r') as fTargetsScored:
                for targetScored in [x.split('\t') for x in fTargetsScored.readlines()]:
                    if len(targetScored) == 2:
                        targetsScored[targetScored[0]] = float(targetScored[1].strip())
            
            failedCount = 0
            for target23 in possibleTargets:
                if target23[0:20] in targetsScored:
                    score = targetsScored[target23[0:20]]
                    possibleTargets[target23]['offtargetscore'] = score
                    
                    if score < float(configMngr['offtargetscore']['score-threshold']):
                        possibleTargets[target23]['passedOffTargetScore'] = CODE_REJECTED # reject due to this reason
                        failedCount += 1
                    else:
                        possibleTargets[target23]['passedOffTargetScore'] = CODE_ACCEPTED # accept due to this reason
            
            printer('\t{} of {} failed here.'.format(
                failedCount,
                testedCount
            ))
            
            #### END: REMOVE OFF-TARGET SCORING
          
        #########################################
        ##           Begin output              ##
        #########################################   
        printer('Writing results to file.')

        # Write guides to file. Include scores etc.
        with open(configMngr['output']['file'], 'a+') as fOpen:
            csvWriter = csv.writer(fOpen, delimiter=configMngr['output']['delimiter'],
                            quotechar='"', quoting=csv.QUOTE_MINIMAL)
                         
            csvWriter.writerow(DEFAULT_GUIDE_PROPERTIES_ORDER)
            
            for target23 in possibleTargets:
                output = [possibleTargets[target23][x] for x in DEFAULT_GUIDE_PROPERTIES_ORDER]
                
                csvWriter.writerow(output)

        #########################################
        ##              Clean up               ##
        #########################################
        printer('Cleaning auxiliary files')
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
    
    printer('Total run time (dd hh:mm:ss) {} or {} seconds'.format(
        strftime("%d %H:%M:%S", gmtime(totalRunTimeSec)), 
        totalRunTimeSec
    ))   
    
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