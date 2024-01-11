'''
https://github.com/bmds-lab/Crackling

Author: Jake Bradford, Dimitri Perrin

Config:
    - See config.ini
'''

import ast, csv, joblib, os, re, sys, time, tempfile, shutil

from crackling.Paginator import Paginator
from crackling.Batchinator import Batchinator
from crackling.Constants import *
from crackling.Helpers import *

def Crackling(configMngr):
    totalSizeBytes = configMngr.getDatasetSizeBytes()
    completedSizeBytes = 0

    _stdout = sys.stdout
    _stderr = sys.stderr

    sys.stdout = configMngr.getLogMethod()
    sys.stderr = configMngr.getErrLogMethod()

    lastRunTimeSec = 0
    seqFileSize = 0
    totalRunTimeSec = 0

    startTime = time.time()

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
                if (candidateGuides[target23]['isUnique'] == CODE_REJECTED):
                    doAssess = False

            if optimisation == 'medium':
                # Never assess guides that appear twice
                if (candidateGuides[target23]['isUnique'] == CODE_REJECTED):
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
                if (candidateGuides[target23]['isUnique'] == CODE_REJECTED):
                    doAssess = False

                # For efficiency
                if module in [MODULE_CHOPCHOP, MODULE_MM10DB, MODULE_SGRNASCORER2]:

                    # `consensusCount` cannot be used yet as it may not have been
                    # calculated. instead, calculate it on-the-go.
                    countAlreadyAccepted = sum([
                        candidateGuides[target23]['acceptedByMm10db'] == CODE_ACCEPTED,
                        candidateGuides[target23]['passedG20'] == CODE_ACCEPTED,
                        candidateGuides[target23]['acceptedBySgRnaScorer'] == CODE_ACCEPTED,
                    ])

                    countAlreadyAssessed = sum([
                        candidateGuides[target23]['acceptedByMm10db'] in [CODE_ACCEPTED, CODE_REJECTED],
                        candidateGuides[target23]['passedG20'] in [CODE_ACCEPTED, CODE_REJECTED],
                        candidateGuides[target23]['acceptedBySgRnaScorer'] in [CODE_ACCEPTED, CODE_REJECTED],
                    ])

                    countToolsInConsensus = sum([
                        configMngr['consensus'].getboolean('mm10db'),
                        configMngr['consensus'].getboolean('chopchop'),
                        configMngr['consensus'].getboolean('sgRNAScorer2'),
                    ])

                    # Do not assess if passed consensus already
                    if countAlreadyAccepted >= consensusN:
                        doAssess = False

                    # Do not assess if there are not enough remaining tests to pass consensus
                    #   i.e. if the number of remaining tests is less than what is needed to pass then do not assess
                    if countToolsInConsensus - countAlreadyAssessed < consensusN - countAlreadyAccepted:
                        doAssess = False

                    # For mm10db:
                    if module == MODULE_MM10DB:
                        # if any of the mm10db tests have failed, then fail them all
                        if (CODE_REJECTED in [
                            candidateGuides[target23]['passedAvoidLeadingT'],
                            candidateGuides[target23]['passedATPercent'],
                            candidateGuides[target23]['passedTTTT'],
                            candidateGuides[target23]['passedSecondaryStructure'],
                            candidateGuides[target23]['acceptedByMm10db'],
                        ]):
                            doAssess = False

                # For specificity:
                if (module == MODULE_SPECIFICITY):
                    # don't assess if they failed consensus
                    if (int(candidateGuides[target23]['consensusCount']) < consensusN):
                        doAssess = False

                    # don't assess if they failed Bowtie
                    if (candidateGuides[target23]['passedBowtie'] == CODE_REJECTED):
                        doAssess = False

            if doAssess:
                yield target23

    def processSequence(sequence):
        # Patterns for guide matching
        pattern_forward = r'(?=([ATCG]{21}GG))'
        pattern_reverse = r'(?=(CC[ACGT]{21}))'

        # New sequence deteced, process sequence
        # once for forward, once for reverse
        for pattern, strand, seqModifier in [
            [pattern_forward, '+', lambda x : x],
            [pattern_reverse, '-', lambda x : rc(x)]
        ]:
            p = re.compile(pattern)
            for m in p.finditer(sequence):
                target23 = seqModifier(seq[m.start() : m.start() + 23])
                yield [target23, seqHeader,  m.start(),  m.start() + 23, strand]

    ###################################
    ##   Processing the input file   ##
    ###################################

    printer('Analysing files...')

    # Sets to keep track of Guides and sequences seen before
    candidateGuides = set()
    duplicateGuides = set()
    recordedSequences = set()

    guideBatchinator = Batchinator(int(configMngr['input']['batch-size']))

    printer(f'Batchinator is writing to: {guideBatchinator.workingDir.name}')

    for seqFilePath in configMngr.getIterFilesToProcess():

        numIdentifiedGuides = 0
        numDuplicateGuides = 0

        printer(f'Identifying possible target sites in: {seqFilePath}')
        seqFileSize = os.path.getsize(seqFilePath)

        completedSizeBytes += seqFileSize

        # We first remove all the line breaks within a given sequence (FASTA format)
        with open(seqFilePath, 'r') as inFile, tempfile.NamedTemporaryFile(mode='w',delete=False) as parsedFile:
            for line in inFile:
                line = line.strip()
                if line[0] == '>':
                    # this is the header line for a new sequence, so we break the previous line and write the header as a new line
                    parsedFile.write('\n'+line+'\n')
                else:
                    # this is (part of) the sequence; we write it without line break
                    parsedFile.write(line.strip())


        with open(parsedFile.name, 'r') as inFile:
            seqHeader = ''
            seq = ''
            for line in inFile:
                # Remove garbage from line
                line = line.strip()
                # Some lines (e.g., first line in file) can be just a line break, move to next line
                if line=='':
                    continue
                # Header line, start of a new sequence
                elif line[0]=='>':
                    # If we haven't seen the sequence OR we have found a sequence without header
                    if (seqHeader not in recordedSequences) or (seqHeader=='' and seq!=''):
                        # Record header
                        recordedSequences.add(seqHeader)
                        # Process the sequence
                        for guide in processSequence(seq):
                            numIdentifiedGuides += 1
                            # Check if guide has been seen before
                            if guide[0] not in candidateGuides:
                                # Record guide
                                candidateGuides.add(guide[0])
                                # Record candidate guide to temp file
                                guideBatchinator.recordEntry(guide)
                            else:
                                # Record duplicate guide
                                duplicateGuides.add(guide[0])
                                numDuplicateGuides += 1
                    # Update sequence and sequence header
                    seqHeader = line[1:]
                    seq = ''
                # Sequence line, section of existing sequence
                else:
                    # Append section to total sequence
                    seq += line.strip()

            # Process the last sequence
            for guide in processSequence(seq):
                numIdentifiedGuides += 1
                # Check if guide has been seen before
                if guide[0] not in candidateGuides:
                    # Record guide
                    candidateGuides.add(guide[0])
                    # Record candidate guide to temp file
                    guideBatchinator.recordEntry(guide)
                else:
                    # Record duplicate guide
                    duplicateGuides.add(guide[0])
                    numDuplicateGuides += 1

        duplicatePercent = round(numDuplicateGuides / numIdentifiedGuides * 100.0, 3)
        printer(f'\tIdentified {numIdentifiedGuides:,} possible target sites in this file.')
        printer(f'\tOf these, {len(duplicateGuides):,} are not unique. These sites occur a total of {numDuplicateGuides} times.')
        printer(f'\tRemoving {numDuplicateGuides:,} of {numIdentifiedGuides:,} ({duplicatePercent}%) guides.')
        printer(f'\t{len(candidateGuides):,} distinct guides have been discovered so far.')

        completedPercent = round(completedSizeBytes / totalSizeBytes * 100.0, 3)
        printer(f'\tExtracted from {completedPercent}% of input')

    # Write header line for output file
    with open(configMngr['output']['file'], 'a+') as fOpen:
        csvWriter = csv.writer(fOpen, delimiter=configMngr['output']['delimiter'],
                        quotechar='"',dialect='unix', quoting=csv.QUOTE_MINIMAL)

        csvWriter.writerow(DEFAULT_GUIDE_PROPERTIES_ORDER)

    # Clean up unused variables
    os.unlink(parsedFile.name)
    del candidateGuides
    del recordedSequences


    batchFileId = 0
    for batchFile in guideBatchinator:
        batchStartTime = time.time()

        printer(f'Processing batch file {(batchFileId+1):,} of {len(guideBatchinator)}')

        # Create new candidate guide dictionary
        candidateGuides = {}
        # Load guides from temp file
        with open(batchFile, 'r') as inputFp:
            # Create csv reader to parse temp file
            csvReader = csv.reader(inputFp, delimiter=configMngr['output']['delimiter'],
                quotechar='"',dialect='unix', quoting=csv.QUOTE_MINIMAL)
            # Rebuild dictonary from temp file
            for row in csvReader:
                candidateGuides[row[0]] = DEFAULT_GUIDE_PROPERTIES.copy()
                candidateGuides[row[0]]['seq'] = row[0]
                if row[0] in duplicateGuides:
                    candidateGuides[row[0]]['header'] = CODE_AMBIGUOUS
                    candidateGuides[row[0]]['start'] = CODE_AMBIGUOUS
                    candidateGuides[row[0]]['end'] = CODE_AMBIGUOUS
                    candidateGuides[row[0]]['strand'] = CODE_AMBIGUOUS
                    candidateGuides[row[0]]['isUnique'] = CODE_REJECTED
                else:
                    candidateGuides[row[0]]['header'] = row[1]
                    candidateGuides[row[0]]['start'] = row[2]
                    candidateGuides[row[0]]['end'] = row[3]
                    candidateGuides[row[0]]['strand'] = row[4]

        printer(f'\tLoaded {len(candidateGuides):,} guides')

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

            printer(f'\t{failedCount:,} of {testedCount:,} failed here.')

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

            printer(f'\t{failedCount:,} of {testedCount:,} failed here.')

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

            printer(f'\t{failedCount:,} of {testedCount:,} failed here.')

        ############################################
        ##   Removing targets that contain TTTT   ##
        ############################################
        if (configMngr['consensus'].getboolean('mm10db')):
            printer('mm10db - remove all targets that contain TTTT.')

            failedCount = 0
            testedCount = 0
            for target23 in filterCandidateGuides(candidateGuides, MODULE_MM10DB):
                if 'TTTT' in target23:
                    candidateGuides[target23]['passedTTTT'] = CODE_REJECTED
                    failedCount += 1
                else:
                    candidateGuides[target23]['passedTTTT'] = CODE_ACCEPTED
                testedCount += 1

            printer(f'\t{failedCount:,} of {testedCount:,} failed here.')

        ##########################################
        ##   Calculating secondary structures   ##
        ##########################################
        if (configMngr['consensus'].getboolean('mm10db')):
            printer('mm10db - check secondary structure.')

            # RNAFold is memory intensive for very large datasets.
            # We will paginate in order not to overflow memory.

            guide = 'GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU'
            pattern_RNAstructure = r'.{28}\({4}\.{4}\){4}\.{3}\){4}.{21}\({4}\.{4}\){4}\({7}\.{3}\){7}\.{3}\s\((.+)\)'
            pattern_RNAenergy = r'\s\((.+)\)'

            testedCount = 0
            failedCount = 0
            errorCount = 0
            notFoundCount = 0

            pgLength = int(configMngr['rnafold']['page-length'])

            for pgIdx, pageCandidateGuides in Paginator(
                filterCandidateGuides(candidateGuides, MODULE_MM10DB),
                pgLength
            ):
                if pgLength > 0:
                    printer(f'\tProcessing page {(pgIdx+1)} ({pgLength:,} per page).')

                if os.path.exists(configMngr['rnafold']['output']):
                    os.remove(configMngr['rnafold']['output'])

                printer('\t\tConstructing the RNAfold input file.')

                guidesInPage = 0
                with open(configMngr['rnafold']['input'], 'w+') as fRnaInput:
                    for target23 in pageCandidateGuides:
                        fRnaInput.write(f'G{target23[1:20]}{guide}\n')
                        guidesInPage += 1

                printer(f'\t\t{guidesInPage:,} guides in this page.')

                runner('{} --noPS -j{} -i {} -o'.format(
                        configMngr['rnafold']['binary'],
                        configMngr['rnafold']['threads'],
                        configMngr['rnafold']['input']
                    ),
                    shell=True,
                    check=True
                )

                os.replace('RNAfold_output.fold' ,configMngr['rnafold']['output'])

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


                for target23 in pageCandidateGuides:
                    key = target23[1:20]
                    if key not in RNAstructures:
                        print(f'Could not find: {target23[0:20]}')
                        notFoundCount += 1
                        continue
                    else:
                        L1 = RNAstructures[key][0]
                        L2 = RNAstructures[key][1]
                        target = RNAstructures[key][2]

                    structure = L2.split(' ')[0]
                    energy = L2.split(' ')[1][1:-1]

                    candidateGuides[target23]['ssL1'] = L1
                    candidateGuides[target23]['ssStructure'] = structure
                    candidateGuides[target23]['ssEnergy'] = energy

                    if transToDNA(target) != target23[0:20] and transToDNA('C'+target[1:]) != target23[0:20] and transToDNA('A'+target[1:]) != target23[0:20]:
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
                    testedCount += 1


            printer(f'\t{failedCount:,} of {testedCount:,} failed here.')

            if errorCount > 0:
                printer(f'\t{errorCount} of {testedCount} erred here.')

            if notFoundCount > 0:
                printer(f'\t{notFoundCount} of {testedCount} not found in RNAfold output.')

        #########################################
        ##         Calc mm10db result          ##
        #########################################
        if (configMngr['consensus'].getboolean('mm10db')):
            printer('Calculating mm10db final result.')

            acceptedCount = 0
            failedCount = 0

            for target23 in candidateGuides:
                if not all([
                    candidateGuides[target23]['passedATPercent'] == CODE_ACCEPTED,
                    candidateGuides[target23]['passedTTTT'] == CODE_ACCEPTED,
                    candidateGuides[target23]['passedSecondaryStructure'] == CODE_ACCEPTED,
                    candidateGuides[target23]['passedAvoidLeadingT'] == CODE_ACCEPTED,
                ]):
                    # mm10db rejected the guide
                    candidateGuides[target23]['acceptedByMm10db'] = CODE_REJECTED
                    failedCount += 1
                else:
                    candidateGuides[target23]['acceptedByMm10db'] = CODE_ACCEPTED
                    acceptedCount += 1

            printer(f'\t{acceptedCount} accepted.')

            printer(f'\t{failedCount} failed.')

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

            printer(f'\t{failedCount:,} of {testedCount:,} failed here.')

        #########################################
        ##      Begin efficacy consensus       ##
        #########################################
        printer('Evaluating efficiency via consensus approach.')

        failedCount = 0
        testedCount = 0
        for target23 in candidateGuides:
            candidateGuides[target23]['consensusCount'] = sum([
                candidateGuides[target23]['acceptedByMm10db'] == CODE_ACCEPTED,
                candidateGuides[target23]['acceptedBySgRnaScorer'] == CODE_ACCEPTED,
                candidateGuides[target23]['passedG20'] == CODE_ACCEPTED,
            ])

            if candidateGuides[target23]['consensusCount'] < int(configMngr['consensus']['n']):
                failedCount += 1

            testedCount += 1

        printer(f'\t{failedCount:,} of {testedCount:,} failed here.')

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
                    printer(f'\tProcessing page {(pgIdx+1)} ({pgLength:,} per page).')

                if os.path.exists(configMngr['bowtie2']['output']):
                    os.remove(configMngr['bowtie2']['output'])

                printer('\tConstructing the Bowtie input file.')

                tempTargetDict_offset = {}
                guidesInPage = 0
                with open(configMngr['bowtie2']['input'], 'w') as fWriteBowtie:
                    for target23 in pageCandidateGuides:
                        similarTargets = [
                            target23[0:20] + 'AGG',
                            target23[0:20] + 'CGG',
                            target23[0:20] + 'GGG',
                            target23[0:20] + 'TGG',
                            target23[0:20] + 'AAG',
                            target23[0:20] + 'CAG',
                            target23[0:20] + 'GAG',
                            target23[0:20] + 'TAG'
                        ]

                        for seq in similarTargets:
                            fWriteBowtie.write(seq + '\n')
                            tempTargetDict_offset[seq] = target23

                        testedCount += 1
                        guidesInPage += 1

                printer(f'\t\t{guidesInPage:,} guides in this page.')

                runner('{} -x {} -p {} --reorder --no-hd -t -r -U {} -S {}'.format(
                        configMngr['bowtie2']['binary'],
                        configMngr['input']['bowtie2-index'],
                        configMngr['bowtie2']['threads'],
                        configMngr['bowtie2']['input'],
                        configMngr['bowtie2']['output']
                    ),
                    shell=True,
                    check=True
                )

                printer('\tStarting to process the Bowtie results.')

                inFile = open(configMngr['bowtie2']['output'], 'r')
                bowtieLines = inFile.readlines()
                inFile.close()

                i=0
                while i<len(bowtieLines):
                    nb_occurences = 0
                    # we extract the read and use the dictionary to find the corresponding target
                    line = bowtieLines[i].rstrip().split('\t')
                    chr = line[2]
                    pos = ast.literal_eval(line[3])
                    read = line[9]
                    seq = ''

                    if read in tempTargetDict_offset:
                        seq = tempTargetDict_offset[read]
                    elif rc(read) in tempTargetDict_offset:
                        seq = tempTargetDict_offset[rc(read)]
                    else:
                        print('Problem? '+read)

                    if seq[:-2] == 'GG':
                        candidateGuides[seq]['bowtieChr'] = chr
                        candidateGuides[seq]['bowtieStart'] = pos
                        candidateGuides[seq]['bowtieEnd'] = pos + 22
                    elif rc(seq)[:2] == 'CC':
                        candidateGuides[seq]['bowtieChr'] = chr
                        candidateGuides[seq]['bowtieStart'] = pos
                        candidateGuides[seq]['bowtieEnd'] = pos + 22
                    else:
                        print('Error? '+seq)
                        quit()

                    # we count how many of the eight reads for this target have a perfect alignment
                    for j in range(i,i+8):

                        # http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#sam-output
                        # XM:i:<N>    The number of mismatches in the alignment. Only present if SAM record is for an aligned read.
                        # XS:i:<N>    Alignment score for the best-scoring alignment found other than the alignment reported.

                        if 'XM:i:0' in bowtieLines[j]:
                            nb_occurences += 1

                            # we also check whether this perfect alignment also happens elsewhere
                            if 'XS:i:0'  in bowtieLines[j]:
                                nb_occurences += 1

                    # if that number is at least two, the target is removed
                    if nb_occurences > 1:

                        # increment the counter if this guide has not already been rejected by bowtie
                        if candidateGuides[seq]['passedBowtie'] != CODE_REJECTED:
                            failedCount += 1

                        candidateGuides[seq]['passedBowtie'] = CODE_REJECTED
                    else:
                        candidateGuides[seq]['passedBowtie'] = CODE_ACCEPTED

                    # we continue with the next target
                    i+=8

                # we can remove the dictionary
                del tempTargetDict_offset

            printer(f'\t{failedCount:,} of {testedCount:,} failed here.')

            #########################################
            ##      Begin off-target scoring       ##
            #########################################
            printer('Beginning off-target scoring.')

            testedCount = 0
            failedCount = 0

            pgLength = int(configMngr['offtargetscore']['page-length'])

            for pgIdx, pageCandidateGuides in Paginator(
                filterCandidateGuides(candidateGuides, MODULE_SPECIFICITY),
                pgLength
            ):

                if pgLength > 0:
                    printer(f'\tProcessing page {(pgIdx+1)} ({pgLength:,} per page).')

                # prepare the list of candidate guides to score
                guidesInPage = 0
                with open(configMngr['offtargetscore']['input'], 'w') as fTargetsToScore:
                    for target23 in pageCandidateGuides:
                        target = target23[0:20]
                        fTargetsToScore.write(target+'\n')
                        testedCount += 1
                        guidesInPage += 1

                if guidesInPage != pgLength:
                    printer(f'\t\t{guidesInPage:,} guides in this page.')

                # Convert line endings (Windows)
                if os.name == 'nt':
                    runner('dos2unix {}'.format(
                            configMngr['offtargetscore']['input']
                        ),
                        shell=True,
                        check=True
                    )

                # call the scoring method
                runner('{} {} {} {} {} {} > {}'.format(
                        configMngr['offtargetscore']['binary'],
                        configMngr['input']['offtarget-sites'],
                        configMngr['offtargetscore']['input'],
                        str(configMngr['offtargetscore']['max-distance']),
                        str(configMngr['offtargetscore']['score-threshold']),
                        str(configMngr['offtargetscore']['method']),
                        configMngr['offtargetscore']['output'],
                    ),
                    shell=True,
                    check=True
                )

                targetsScored = {}
                with open(configMngr['offtargetscore']['output'], 'r') as fTargetsScored:
                    for targetScored in [x.split('\t') for x in fTargetsScored.readlines()]:
                        if len(targetScored) == 3:
                            targetsScored[targetScored[0]] = {'MIT': -1.0, 'CFD': -1.0}
                            targetsScored[targetScored[0]]['MIT'] = float(targetScored[1].strip())
                            targetsScored[targetScored[0]]['CFD'] = float(targetScored[2].strip())

                failedCount = 0
                for target23 in pageCandidateGuides:
                    if target23[0:20] in targetsScored:
                        score = targetsScored[target23[0:20]]
                        candidateGuides[target23]['mitOfftargetscore'] = score['MIT']
                        candidateGuides[target23]['cfdOfftargetscore'] = score['CFD']
                        scoreThreshold = float(configMngr['offtargetscore']['score-threshold'])
                        scoreMethod = str(configMngr['offtargetscore']['method']).strip().lower()

                        # MIT
                        if scoreMethod == 'mit':
                            if score['MIT'] < scoreThreshold:
                                candidateGuides[target23]['passedOffTargetScore'] = CODE_REJECTED
                                failedCount += 1
                            else:
                                candidateGuides[target23]['passedOffTargetScore'] = CODE_ACCEPTED

                        # CFD
                        elif scoreMethod == 'cfd':
                            if score['CFD'] < scoreThreshold:
                                candidateGuides[target23]['passedOffTargetScore'] = CODE_REJECTED
                                failedCount += 1
                            else:
                                candidateGuides[target23]['passedOffTargetScore'] = CODE_ACCEPTED

                        # AND
                        elif scoreMethod == 'and':
                            if (score['MIT'] < scoreThreshold) and (score['CFD'] < scoreThreshold):
                                candidateGuides[target23]['passedOffTargetScore'] = CODE_REJECTED
                                failedCount += 1
                            else:
                                candidateGuides[target23]['passedOffTargetScore'] = CODE_ACCEPTED

                        # OR
                        elif scoreMethod == 'or':
                            if (score['MIT'] < scoreThreshold) or (score['CFD'] < scoreThreshold):
                                candidateGuides[target23]['passedOffTargetScore'] = CODE_REJECTED
                                failedCount += 1
                            else:
                                candidateGuides[target23]['passedOffTargetScore'] = CODE_ACCEPTED

                        # AVERAGE
                        elif scoreMethod == 'avg':
                            if ((score['MIT'] + score['CFD'])/2) < scoreThreshold:
                                candidateGuides[target23]['passedOffTargetScore'] = CODE_REJECTED
                                failedCount += 1
                            else:
                                candidateGuides[target23]['passedOffTargetScore'] = CODE_ACCEPTED

                printer(f'\t{failedCount:,} of {testedCount:,} failed here.')

        #########################################
        ##           Begin output              ##
        #########################################
        printer('Writing results to file.')

        # Write guides to file. Include scores etc.
        with open(configMngr['output']['file'], 'a+') as fOpen:
            csvWriter = csv.writer(fOpen, delimiter=configMngr['output']['delimiter'],
                            quotechar='"',dialect='unix', quoting=csv.QUOTE_MINIMAL)

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

        printer(f'{len(candidateGuides)} guides evaluated.')

        batchEndTime = time.time()
        timeDelta = batchEndTime - batchStartTime

        printer(f'This batch ran in {elapsedTimeString(batchStartTime, batchEndTime)} (dd hh:mm:ss) or {timeDelta:.3} seconds')

        batchFileId += 1

    endTime = time.time()
    timeDelta = endTime - startTime

    printer(f'Total run time {elapsedTimeString(startTime, endTime)} (dd hh:mm:ss) or {timeDelta:.3} seconds')

    sys.stdout.log.close()
    sys.stderr.log.close()
    sys.stdout = _stdout
    sys.stderr = _stderr

if __name__ == '__main__':
    print('This file is not callable')
    exit(1)
