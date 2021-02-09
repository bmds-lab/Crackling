'''
https://github.com/bmds-lab/Crackling

Author: Jake Bradford, Dimitri Perrin

Config:
    - See config.ini
'''

import argparse, ast, csv, joblib, os, re, sys, time
from sklearn.svm import SVC

from ConfigManager import ConfigManager
from Paginator import Paginator
from Constants import *
from Helpers import * 

def Crackling(configMngr):
    totalSizeBytes = configMngr.getDatasetSizeBytes()
    completedSizeBytes = 0

    sys.stdout = configMngr.getLogMethod()

    lastRunTimeSec = 0
    lastScaffoldSizeBytes = 0
    totalRunTimeSec = 0

    candidateGuides = {}

    start_time = time.time()

    ####################################
    ###     Run-time Optimisation     ##
    ####################################
    def filterCandidateGuides(dictCandidateGuides, module):
        nonlocal configMngr
        
        module = module.lower()
        
        optimisation = configMngr['general']['optimisation']
        
        consensusN = int(configMngr['consensus']['n'])
        
        for target23 in dictCandidateGuides:
            doAssess = True
        
            if optimisation == 'ultralow':
                doAssess = True
        
            if optimisation == 'low':
                # Never assess guides that appear twice
                if (candidateGuides[target23]['seqCount'] > 1):
                    doAssess = False
                
            if optimisation == 'medium':
                # Never assess guides that appear twice
                if (candidateGuides[target23]['seqCount'] > 1):
                    doAssess = False
            
                # For mm10db:
                if (module == MODULE_MM10DB):
                    # if any of the mm10db tests have failed, then fail them all
                    if (CODE_REJECTED in [
                        candidateGuides[target23]['passedAvoidLeadingT'],
                        candidateGuides[target23]['passedATPercent'],
                        candidateGuides[target23]['passedTTTT'],
                        candidateGuides[target23]['passedSecondaryStructure'],
                        candidateGuides[target23]['acceptedByMm10db'],
                    ]):
                        doAssess = False

                # For CHOPCHOP:
                # Always assess, unless the guide is seen multiple times
                        
                # For sgRNAScorer2:
                # Always assess, unless the guide is seen multiple times
                
                # For specificity:
                if (module == MODULE_SPECIFICITY):
                    # don't assess if they failed consensus
                    if (int(candidateGuides[target23]['consensusCount']) < consensusN):
                        doAssess = False
                    
                    # don't assess if they failed Bowtie
                    if (candidateGuides[target23]['passedBowtie'] == CODE_REJECTED):
                        doAssess = False
                
            if optimisation == 'high':
                # Never assess guides that appear twice
                if (candidateGuides[target23]['seqCount'] > 1):
                    doAssess = False
                    
                # For mm10db:
                if (module == MODULE_MM10DB):
                    # if any of the mm10db tests have failed, then fail them all
                    if (CODE_REJECTED in [
                        candidateGuides[target23]['passedAvoidLeadingT'],
                        candidateGuides[target23]['passedATPercent'],
                        candidateGuides[target23]['passedTTTT'],
                        candidateGuides[target23]['passedSecondaryStructure'],
                        candidateGuides[target23]['acceptedByMm10db'],
                    ]):
                        doAssess = False

                # For CHOPCHOP:
                if (module == MODULE_CHOPCHOP):
                    if consensusN == 1 and candidateGuides[target23]['acceptedByMm10db'] == CODE_ACCEPTED:
                        doAssess = False
                        
                # For sgRNAScorer2:
                if (module == MODULE_SGRNASCORER2):
                    currentConsensus = ((int)(candidateGuides[target23]['acceptedByMm10db'] == CODE_ACCEPTED) +    # mm10db accepted
                    (int)(candidateGuides[target23]['passedG20'] == CODE_ACCEPTED))                            # chopchop-g20 accepted

                    # if the guide is further than one test away from passing
                    # the consensus approach, then skip it.
                    if (currentConsensus < (consensusN - 1)):
                        doAssess = False
                        
                # For specificity:
                if (module == MODULE_SPECIFICITY):
                    # don't assess if they failed consensus
                    if (int(candidateGuides[target23]['consensusCount']) < consensusN):
                        doAssess = False

            if doAssess:
                yield target23
    
    ###################################
    ##   Processing the input file   ##
    ###################################
    printer('Analysing files...') 
    
    # key: FASTA header, value: sequence
    seqsByHeader = {}
    
    for seqFilePath in configMngr.getIterFilesToProcess():
        printer('Identifying possible target sites in: {}'.format(
            seqFilePath
        ))

        sys.stdout.log.flush()

        printer('{} of {} bytes processed ({}%)'.format(
            completedSizeBytes,
            totalSizeBytes,
            round((float(completedSizeBytes) / float(totalSizeBytes) * 100.0), 3)
        ))

        lastScaffoldSizeBytes = os.path.getsize(seqFilePath)

        completedSizeBytes += lastScaffoldSizeBytes

        with open(seqFilePath, 'r') as inFile:
            previousHeader = seqFilePath
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

        pattern_forward = r"(?=([ATCG]{21}GG))"
        pattern_reverse = r"(?=(CC[ACGT]{21}))"

        for seqHeader in seqsByHeader:
            seq = seqsByHeader[seqHeader]
            
            # once for forward, once for reverse
            for pattern, strand, seqModifier in [
                [pattern_forward, '+', lambda x : x], 
                [pattern_reverse, '-', lambda x : rc(x)]
            ]:
                p = re.compile(pattern)
                for m in p.finditer(seq):
                    target23 = seqModifier(seq[m.start() : m.start() + 23])
                    if target23 not in candidateGuides:
                        candidateGuides[target23] = DEFAULT_GUIDE_PROPERTIES.copy()
                        candidateGuides[target23]['seq'] = target23
                        candidateGuides[target23]['header'] = seqHeader
                        candidateGuides[target23]['start'] = m.start()
                        candidateGuides[target23]['end'] = m.start() + 23
                        candidateGuides[target23]['strand'] = strand
                    else:
                        # we've already seen this guide, make the positioning ambiguous
                        candidateGuides[target23]['seqCount'] += 1
                        candidateGuides[target23]['header'] = CODE_AMBIGUOUS
                        candidateGuides[target23]['start'] = CODE_AMBIGUOUS
                        candidateGuides[target23]['end'] = CODE_AMBIGUOUS
                        candidateGuides[target23]['strand'] = CODE_AMBIGUOUS

    del seqsByHeader

    printer('Identified {} possible target sites.'.format(
        len(candidateGuides)
    ))

    ############################################
    ##     Removing targets with leading T    ##
    ############################################
    if (configMngr['consensus'].getboolean('mm10db')):
        printer('mm10db - remove all targets with a leading T (+) or trailing A (-).')

        failedCount = 0
        testedCount = 0
        for target23 in filterCandidateGuides(candidateGuides, MODULE_MM10DB):
            if (target23[-2:] == 'GG' and target23[0] == 'T') or \
                (target23[:2] == 'CC' and target23[-1] == 'A'):
                candidateGuides[target23]['passedAvoidLeadingT'] = CODE_REJECTED
                failedCount += 1
            else:
                candidateGuides[target23]['passedAvoidLeadingT'] = CODE_ACCEPTED
            
            testedCount += 1
            
        printer('\t{} of {} failed here.'.format(
            failedCount,
            testedCount
        ))

    #########################################
    ##    AT% ideally is between 20-65%    ##
    ######################################### 
    if (configMngr['consensus'].getboolean('mm10db')):
        printer('mm10db - remove based on AT percent.')

        failedCount = 0
        testedCount = 0
        for target23 in filterCandidateGuides(candidateGuides, MODULE_MM10DB):
            AT = AT_percentage(target23[0:20])

            if AT < 20 or AT > 65:
                candidateGuides[target23]['passedATPercent'] = CODE_REJECTED
                failedCount += 1
            else:
                candidateGuides[target23]['passedATPercent'] = CODE_ACCEPTED
            
            candidateGuides[target23]['AT'] = AT
            
            testedCount += 1
            
        printer('\t{} of {} failed here.'.format(
            failedCount,
            testedCount
        ))

    ############################################
    ##   Removing targets that contain TTTT   ##
    ############################################
    if (configMngr['consensus'].getboolean('mm10db')):
        printer('mm10db - remove all targets that contain TTTT.')

        failedCount = 0
        testedCount = 0
        for target23 in filterCandidateGuides(candidateGuides, MODULE_MM10DB):
            if "TTTT" in target23:
                candidateGuides[target23]['passedTTTT'] = CODE_REJECTED
                failedCount += 1
            else:
                candidateGuides[target23]['passedTTTT'] = CODE_ACCEPTED
           
            testedCount += 1
            
        printer('\t{} of {} failed here.'.format(
            failedCount,
            testedCount
        ))

    ##########################################
    ##   Calculating secondary structures   ##
    ##########################################
    if (configMngr['consensus'].getboolean('mm10db')):
        printer('mm10db - check secondary structure.')
        
        # RNAFold is memory intensive for very large datasets.
        # We will paginate in order not to overflow memory.

        guide = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU"
        pattern_RNAstructure = r".{28}\({4}\.{4}\){4}\.{3}\){4}.{21}\({4}\.{4}\){4}\({7}\.{3}\){7}\.{3}\s\((.+)\)"
        pattern_RNAenergy = r"\s\((.+)\)"

        testedCount = 0
        failedCount = 0
        errorCount = 0
        notFoundCount = 0 

        pgLength = int(configMngr['rnafold']['page-length'])

        for pgIdx, pageCandidateGuides in Paginator(candidateGuides, pgLength):
            
            if pgLength > 0:
                numPages = int(len(candidateGuides) / pgLength) + 1
                printer(f'\tProcessing page {(pgIdx+1)} of {numPages}')
            
            if os.path.exists(configMngr['rnafold']['output']):
                os.remove(configMngr['rnafold']['output'])

            printer('\t\tConstructing the RNAfold input file.')

            with open(configMngr['rnafold']['input'], 'w+') as fRnaInput:
                for target23 in candidateGuides:
                    # we only guides that have passed for mm10db's test so far
                    if  candidateGuides[target23]['passedAvoidLeadingT'] == 1 and \
                        candidateGuides[target23]['passedATPercent'] == 1 and \
                        candidateGuides[target23]['passedTTTT'] == 1:
                        fRnaInput.write(
                            "G{}{}\n".format(
                                target23[1:20], 
                                guide
                            )
                        )

            caller(
                "{} --noPS -j{} -i \"{}\" > \"{}\"".format(
                    configMngr['rnafold']['binary'],
                    configMngr['rnafold']['threads'],
                    configMngr['rnafold']['input'],
                    configMngr['rnafold']['output']
                ), 
                shell=True
            )

            printer('\t\tStarting to process the RNAfold results.')

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


            for target23 in candidateGuides:
                # we only guides that have passed for mm10db's test so far
                if  candidateGuides[target23]['passedAvoidLeadingT'] == 1 and \
                    candidateGuides[target23]['passedATPercent'] == 1 and \
                    candidateGuides[target23]['passedTTTT'] == 1:
                    
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
                    
                    candidateGuides[target23]['ssL1'] = L1
                    candidateGuides[target23]['ssStructure'] = structure
                    candidateGuides[target23]['ssEnergy'] = energy
                    
                    if transToDNA(target) != target23[0:20] and transToDNA("C"+target[1:]) != target23[0:20] and transToDNA("A"+target[1:]) != target23[0:20]:
                        candidateGuides[target23]['passedSecondaryStructure'] = CODE_ERROR
                        errorCount += 1
                        continue

                    match_structure = re.search(pattern_RNAstructure, L2)
                    if match_structure:
                        energy = ast.literal_eval(match_structure.group(1))
                        if energy < float(configMngr['rnafold']['low_energy_threshold']):
                            candidateGuides[transToDNA(target23)]['passedSecondaryStructure'] = CODE_REJECTED
                            failedCount += 1
                        else:
                            candidateGuides[target23]['passedSecondaryStructure'] = CODE_ACCEPTED
                    else:
                        match_energy = re.search(pattern_RNAenergy, L2)
                        if match_energy:
                            energy = ast.literal_eval(match_energy.group(1))
                            if energy <= float(configMngr['rnafold']['high_energy_threshold']):
                                candidateGuides[transToDNA(target23)]['passedSecondaryStructure'] = CODE_REJECTED
                                failedCount += 1
                            else:
                                candidateGuides[target23]['passedSecondaryStructure'] = CODE_ACCEPTED


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

    #########################################
    ##         Calc mm10db result          ##
    ######################################### 
    if (configMngr['consensus'].getboolean('mm10db')):
        printer('Calculating mm10db final result.')

        acceptedCount = 0
        failedCount = 0
        
        for target23 in candidateGuides:
            if (CODE_REJECTED in [
                candidateGuides[target23]['passedATPercent'],
                candidateGuides[target23]['passedTTTT'],
                candidateGuides[target23]['passedSecondaryStructure'],
                candidateGuides[target23]['passedAvoidLeadingT'],
            ]):
                # mm10db rejected the guide
                candidateGuides[target23]['acceptedByMm10db'] = CODE_REJECTED
                failedCount += 1
            else:
                candidateGuides[target23]['acceptedByMm10db'] = CODE_ACCEPTED
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
        
        # binary encoding
        encoding = {
            'A' : '0001',    'C' : '0010',    'T' : '0100',    'G' : '1000',
            'K' : '1100',    'M' : '0011',    'R' : '1001',    'Y' : '0110',
            'S' : '1010',    'W' : '0101',    'B' : '1110',    'V' : '1011',
            'H' : '0111',    'D' : '1101',    'N' : '1111'
        }

        clfLinear = joblib.load(configMngr['sgrnascorer2']['model'])

        failedCount = 0
        testedCount = 0
        for target23 in filterCandidateGuides(candidateGuides, MODULE_SGRNASCORER2):
            sequence = target23.upper()
            entryList = []
            testedCount += 1

            for x in range(0, 20):
                for y in range(0, 4):
                    entryList.append(int(encoding[sequence[x]][y]))

            # predict based on the entry
            prediction = clfLinear.predict([entryList])
            score = clfLinear.decision_function([entryList])[0] 

            candidateGuides[target23]['sgrnascorer2score'] = score

            if float(score) < float(float(configMngr['sgrnascorer2']['score-threshold'])):
                candidateGuides[target23]['acceptedBySgRnaScorer'] = CODE_REJECTED
                failedCount += 1
            else:
                candidateGuides[target23]['acceptedBySgRnaScorer'] = CODE_ACCEPTED
                        
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
        testedCount = 0
        for target23 in filterCandidateGuides(candidateGuides, MODULE_CHOPCHOP):
            if target23[19] != 'G':
                candidateGuides[target23]['passedG20'] = CODE_REJECTED
                failedCount += 1
            else:
                candidateGuides[target23]['passedG20'] = CODE_ACCEPTED
            testedCount += 1

        printer('\t{} of {} failed here.'.format(
            failedCount,
            testedCount
        ))     

    #########################################
    ##      Begin efficacy consensus       ##
    ######################################### 
    printer('Evaluating efficiency via consensus approach.')
    
    failedCount = 0
    testedCount = 0
    for target23 in candidateGuides:
        candidateGuides[target23]['consensusCount'] = (
            (int)(candidateGuides[target23]['acceptedByMm10db'] == CODE_ACCEPTED)          +     # mm10db accepted
            (int)(candidateGuides[target23]['acceptedBySgRnaScorer'] == CODE_ACCEPTED)     +     # sgrnascorer2 accepted
            (int)(candidateGuides[target23]['passedG20'] == CODE_ACCEPTED)                       # chopchop-g20 accepted
        )
        if candidateGuides[target23]['consensusCount'] < int(configMngr['consensus']['n']):
            failedCount += 1
            
        testedCount += 1
            
    printer('\t{} of {} failed here.'.format(
        failedCount,
        testedCount
    ))   

    if (configMngr['offtargetscore'].getboolean('enabled')):
        ###############################################
        ##         Using Bowtie for positioning      ##
        ###############################################
        printer('Bowtie analysis.')

        testedCount = 0
        failedCount = 0

        pgLength = int(configMngr['bowtie2']['page-length'])

        for pgIdx, pageCandidateGuides in Paginator(
            filterCandidateGuides(candidateGuides, MODULE_SPECIFICITY), 
            pgLength
        ):

            if pgLength > 0:
                numPages = int(len(candidateGuides) / pgLength) + 1
                printer(f'\tProcessing page {(pgIdx+1)} of {numPages}')
            
            if os.path.exists(configMngr['bowtie2']['output']):
                os.remove(configMngr['bowtie2']['output'])

            printer('\tConstructing the Bowtie input file.')
            
            tempTargetDict_offset = {}
            with open(configMngr['bowtie2']['input'], 'w') as fWriteBowtie:
                for target23 in pageCandidateGuides:
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
            , shell=True)       
            
            printer('\tStarting to process the Bowtie results.')
            
            inFile = open(configMngr['bowtie2']['output'], 'r')
            bowtieLines = inFile.readlines()
            inFile.close()

            i=0
            while i<len(bowtieLines):
                nb_occurences = 0
                # we extract the read and use the dictionary to find the corresponding target
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
                    candidateGuides[seq]['bowtieChr'] = chr
                    candidateGuides[seq]['bowtieStart'] = pos
                    candidateGuides[seq]['bowtieEnd'] = pos + 22
                elif rc(seq)[:2] == 'CC':
                    candidateGuides[seq]['bowtieChr'] = chr
                    candidateGuides[seq]['bowtieStart'] = pos
                    candidateGuides[seq]['bowtieEnd'] = pos + 22
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
                    candidateGuides[seq]['passedBowtie'] = CODE_REJECTED
                    failedCount += 1
                else:
                    candidateGuides[seq]['passedBowtie'] = CODE_ACCEPTED
                    
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

        testedCount = 0
        
        # prepare the list of candidate guides to score
        with open(configMngr['offtargetscore']['input'], 'w') as fTargetsToScore:
            for target23 in filterCandidateGuides(candidateGuides, MODULE_SPECIFICITY):
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
        for target23 in filterCandidateGuides(candidateGuides, MODULE_SPECIFICITY):
            if target23[0:20] in targetsScored:
                score = targetsScored[target23[0:20]]
                candidateGuides[target23]['offtargetscore'] = score
                
                if score < float(configMngr['offtargetscore']['score-threshold']):
                    candidateGuides[target23]['passedOffTargetScore'] = CODE_REJECTED
                    failedCount += 1
                else:
                    candidateGuides[target23]['passedOffTargetScore'] = CODE_ACCEPTED
        
        printer('\t{} of {} failed here.'.format(
            failedCount,
            testedCount
        ))

    #########################################
    ##           Begin output              ##
    #########################################   
    printer('Writing results to file.')

    # Write guides to file. Include scores etc.
    with open(configMngr['output']['file'], 'a+') as fOpen:
        csvWriter = csv.writer(fOpen, delimiter=configMngr['output']['delimiter'],
                        quotechar='"', quoting=csv.QUOTE_MINIMAL)
                     
        csvWriter.writerow(DEFAULT_GUIDE_PROPERTIES_ORDER)
        
        for target23 in candidateGuides:
            output = [candidateGuides[target23][x] for x in DEFAULT_GUIDE_PROPERTIES_ORDER]
            
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
        len(candidateGuides)
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
    parser.add_argument('-c', help='Configuration file', default=None, required=True)
    args = parser.parse_args()

    configMngr = ConfigManager(args.c, lambda x : print(f'configMngr says: {x}'))

    if not configMngr.isConfigured():
        print('Something went wrong with reading the configuration.')
        exit()
    else:
        printer('Crackling is starting...')
        
        Crackling(configMngr)
        
    print('Goodbye.')