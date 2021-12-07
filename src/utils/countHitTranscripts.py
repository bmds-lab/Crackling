"""
This utility counts the number of transcripts which a guide will target.
Author: Jake Bradford

Here, an example is provided. Exons are presented by `|||||`

    >      (Guide A)           (Guide B)
    >         *                   *
    >  ----|||*|||--------||||||||*||----|||||||-----  (Transcript 1)
    >         *                   *
    >  ----|||*|||--------||||||||*||----|||||||-----  (Transcript 2)
    >         *                   *
    >  -------*-----------||||||||*||----|||||||-----  (Transcript 3)
    >         *                   *
    >       (2/3)               (3/3)


    i.e. Of the three transcripts, guide A targets two, and guide B targets three.

To achieve this, two inputs are required:
    1. GFF annotation
    2. Crackling-designed guides (provided as Crackling's output file)
    
GFF annotation is formatted as:
    NC_008394.4	RefSeq	gene	169798	172144	.	-	.	ID=gene21;Dbxref=GeneID:9266052;Name=Os01g0102950;gbkey=Gene;gene=Os01g0102950;gene_biotype=protein_coding
    NC_008394.4	RefSeq	mRNA	169798	172144	.	-	.	ID=rna21;Parent=gene21;Dbxref=Genbank:NM_001185188.1,GeneID:9266052;Name=NM_001185188.1;Note=supported by AK108287;gbkey=mRNA;gene=Os01g0102950;transcript_id=NM_001185188.1
    NC_008394.4	RefSeq	exon	172004	172144	.	-	.	ID=id217;Parent=rna21;Dbxref=Genbank:NM_001185188.1,GeneID:9266052;Note=supported by AK108287;gbkey=mRNA;gene=Os01g0102950;transcript_id=NM_001185188.1
    NC_008394.4	RefSeq	exon	171770	171921	.	-	.	ID=id218;Parent=rna21;Dbxref=Genbank:NM_001185188.1,GeneID:9266052;Note=supported by AK108287;gbkey=mRNA;gene=Os01g0102950;transcript_id=NM_001185188.1
    NC_008394.4	RefSeq	exon	171578	171671	.	-	.	ID=id219;Parent=rna21;Dbxref=Genbank:NM_001185188.1,GeneID:9266052;Note=supported by AK108287;gbkey=mRNA;gene=Os01g0102950;transcript_id=NM_001185188.1

The following fields from Crackling output are needed:
    - seq,
    - bowtieChr
    - bowtieStart
    - bowtieEnd

Usage:
    1. Import as a module and call main(..)
    2. Call from the command line
        countHitTranscripts.py --annotation file.gff --crackling results.csv --output results.hits.csv

"""
import os, tempfile, csv, pickle

def loadAnnotation(annotationFile, forceReload=False):
    """
    Load an GFF3 annotation file into memory. Once loaded, the result is Pickle'd to a file.
    
    Args:
        annotationFile: a path to the GFF3 file
        forceReload: True - recompute even if the Pickle file exists. False - load from Pickle file.
    
    """
    annot = {}
    geneData = {}
    geneToMrnaMap = {}
    seqToGeneMap = {}
    geneToSeqMap = {}
    mrnaToGeneMap = {}

    # Check if a Pickle'd version exists
    annotationFilePickled = f"{annotationFile}.p"
    successUnpickling = False
    if os.path.exists(annotationFilePickled) and not forceReload:
        try:
            with open(annotationFilePickled, 'rb') as fp:
                data = pickle.load(fp)
                annot, geneData, geneToMrnaMap, seqToGeneMap, geneToSeqMap, mrnaToGeneMap = data
            successUnpickling = True
        except:
            pass

    if not successUnpickling or forceReload:
        # Load the data
        with open(annotationFile, 'r') as fp:
            i = 0
            for line in fp:
                line = [x.strip() for x in line.split('\t')]
                if len(line) != 9:
                    continue

                seqId, source, type, start, end, score, strand, phase = line[0:8]
                seqId = seqId.replace('.', '_')
                
                attributes = {
                    a.split('=')[0] : a.split('=')[1]
                    for a in line[8].split(';')
                }
                
                # Check the minimum attributes exist
                if any([x not in attributes for x in ['ID', 'Parent']]):
                    continue
                
                # Ignore features that we are not interested in
                if type not in ['gene', 'mRNA', 'exon']:
                    continue
                
                # New sequence (generally a chromosome)
                if seqId not in annot:
                    annot[seqId] = {}

                if type == 'gene':
                    if attributes['ID'] not in geneData:
                        geneData[attributes['ID']] = {
                            'seqId' : seqId,
                            'start' : start,
                            'end' : end,
                            'strand' : strand
                        }
                        
                    if seqId not in seqToGeneMap:
                        seqToGeneMap[seqId] = []
                    seqToGeneMap[seqId].append(attributes['ID'])
                    
                    if attributes['ID'] not in geneToSeqMap:
                        geneToSeqMap[attributes['ID']] = []
                    geneToSeqMap[attributes['ID']].append(seqId)

                if type == 'mRNA':
                    if attributes['ID'] not in annot[seqId]:
                        annot[seqId][attributes['ID']] = []
                        
                    if attributes['Parent'] not in geneToMrnaMap:
                        geneToMrnaMap[attributes['Parent']] = []
                    geneToMrnaMap[attributes['Parent']].append(attributes['ID'])
                    
                    if attributes['ID'] not in mrnaToGeneMap:
                        mrnaToGeneMap[attributes['ID']] = attributes['Parent']
                    
                if type == 'exon':
                    if attributes['Parent'] not in annot[seqId]:
                        annot[seqId][attributes['Parent']] = []
                        
                    annot[seqId][attributes['Parent']].append(
                        (int(start), int(end))
                    )
                    
                i += 1

        # Pickle the data
        data = [annot, geneData, geneToMrnaMap, seqToGeneMap, geneToSeqMap, mrnaToGeneMap] 
        with open(annotationFilePickled, 'wb') as fp:
            pickle.dump(data, fp)
            print(f'Pickled to: {annotationFilePickled}')

    return annot, geneData, geneToMrnaMap, seqToGeneMap, geneToSeqMap, mrnaToGeneMap

def countTranscripts(
    annot, 
    geneData, 
    geneToMrnaMap, 
    seqToGeneMap, 
    geneToSeqMap, 
    mrnaToGeneMap, 
    querySeqId, 
    queryStart, 
    queryEnd
):
    """
    Given annotation data and coordinates of a guide sequence, count the number 
    transcripts that the guide would hit.
    
    Args:
        Use loadAnnotation(..) to obtain most data
        querySeqId - provided by Crackling (field `bowtieChr`)
        queryStart - provided by Crackling (field `bowtieStart`)
        queryEnd - provided by Crackling (field `bowtieEnd`)
    """
    inMrna = []
    if querySeqId in annot:
        for mRNA in annot[querySeqId]:
            skipToNextMRNA = False
            for exon in annot[querySeqId][mRNA]:
                eStart = exon[0]
                eEnd = exon[1]
                if queryStart >= eStart and queryStart <= eEnd:
                    inMrna.append(mRNA)
                    skipToNextMRNA = True
                    break
            if skipToNextMRNA:
                continue

    if len(inMrna) == 0:
        return [0, 0]
   
    # how many transcripts for this gene
    gene = set([mrnaToGeneMap[x] for x in inMrna if x in mrnaToGeneMap])
    if len(gene) > 1:
        raise RuntimeError('Mapped to multiple genes - logical error?')
    else:
        gene = mrnaToGeneMap[inMrna[0]]
        
    return [len(inMrna), len(geneToMrnaMap[gene])]



def process(GFF_FP, CRACKLING_FP):
    annot, geneData, geneToMrnaMap, seqToGeneMap, geneToSeqMap, mrnaToGeneMap = loadAnnotation(GFF_FP, forceReload=True)
    # Read the Crackling file
    CracklingResults = []

    with open(CRACKLING_FP, 'r') as fp:
        fpCsv = csv.reader(fp, delimiter=',', quotechar='"')
        lineNum = 0
        
        idxSeq = None
        idxBowtieChr = None
        idxBowtieStart = None
        idxBowtieEnd = None
        
        for line in fpCsv:
            if lineNum == 0:
                #try:
                idxSeq = line.index('seq')
                idxBowtieChr = line.index('bowtieChr')
                idxBowtieStart = line.index('bowtieStart')
                idxBowtieEnd = line.index('bowtieEnd') 

                line.append('hits')
                
            else:
       
                if line[idxBowtieChr] != '?':
                    seq = line[idxSeq],
                    chr = line[idxBowtieChr]
                    start = int(line[idxBowtieStart])
                    end = int(line[idxBowtieEnd])
                    
                    try:
                        count = countTranscripts(annot, geneData, geneToMrnaMap, seqToGeneMap, geneToSeqMap, mrnaToGeneMap, chr, start, end)
                    except Exception as e:
                        count = ['?', '?']
                        pass
                        
                else:
                    count = ['?', '?']
                    
                line.append(f'{count[0]}/{count[1]}')
                
            CracklingResults.append(line)
            lineNum += 1
            
    return CracklingResults


def useSampleData():
    """
        
        Position indicators, multiply by ten
        0. 0123456789
        1.           0123456789
        2.                     0123456789
        3.                               0123456789
        4.                                         0123456789
        
        Chromosome 1:
                 *             *           *             *                           
           ----||*|||-------|||*|||------||*|||----------*--- (Gene 1 - Transcript 1)
           ----||*|||----------*---------||*|||----------*--- (Gene 1 - Transcript 2)
           ------*----------|||*|||------||*|||----------*--- (Gene 1 - Transcript 3)
           ------*-----------------------||*|||----------*--- (Gene 1 - Transcript 4)
                 *             *           *             *                           

    """

    cracklingData = '''seq,bowtieChr,bowtieStart,bowtieEnd
AAAA,Chr1,60,83
AAAT,Chr1,200,223
AATA,Chr1,320,343
ATAA,Chr1,460,483
'''

    annotationData = '''Chr1	JakeSeq	gene	5	540	.	-	.	ID=gene1
Chr1	JakeSeq	mRNA	10	530	.	-	.	ID=rna1;Parent=gene1
Chr1	JakeSeq	exon	40	100	.	-	.	ID=exon1;Parent=rna1
Chr1	JakeSeq	exon	170	220	.	-	.	ID=exon2;Parent=rna1
Chr1	JakeSeq	exon	300	360	.	-	.	ID=exon3;Parent=rna1
Chr1	JakeSeq	mRNA	50	533	.	-	.	ID=rna2;Parent=gene1
Chr1	JakeSeq	exon	40	100	.	-	.	ID=exon4;Parent=rna2
Chr1	JakeSeq	exon	300	360	.	-	.	ID=exon5;Parent=rna2
Chr1	JakeSeq	mRNA	47	530	.	-	.	ID=rna3;Parent=gene1
Chr1	JakeSeq	exon	170	220	.	-	.	ID=exon6;Parent=rna3
Chr1	JakeSeq	exon	300	360	.	-	.	ID=exon7;Parent=rna3
Chr1	JakeSeq	mRNA	47	530	.	-	.	ID=rna4;Parent=gene1
Chr1	JakeSeq	exon	300	360	.	-	.	ID=exon8;Parent=rna4
'''

    expectedResult = '''AAAA 2/4
AAAT 2/4
AATA 4/4
ATAA 0/0
'''

    print('Writing test data to file.')
    print('The expected results from the test are:')
    print(expectedResult)

    cracklingFp = tempfile.NamedTemporaryFile(mode='w', delete=False)
    annotationFp = tempfile.NamedTemporaryFile(mode='w', delete=False)
    
    #with open(cracklingFp, 'w') as cFp, open(annotationFp, 'w') as aFp:
    cracklingFp.write(cracklingData)
    annotationFp.write(annotationData)
    
    return annotationFp.name, cracklingFp.name

def main():
    import argparse
    
    parser = argparse.ArgumentParser()
    group = parser.add_argument_group('group')
    group.add_argument('-a', '--annotation', help='The GFF3 annotation file', default=None)
    group.add_argument('-c', '--crackling', help='The Crackling output file', default=None)
    group.add_argument('-o', '--output', help='The output file', default=None)
    parser.add_argument('-s', '--sample', help='Run sample', action='store_true', required=False)

    args = parser.parse_args()

    if args.sample:
        for r in process(*useSampleData()):
            print(r)
    else:
    
        results = process(
            args.annotation,
            args.crackling,
        )
        
        with open(args.output, 'w') as fp:
            csvW = csv.writer(
                fp, 
                delimiter=',',
                quotechar='"',
                dialect='unix', 
                quoting=csv.QUOTE_MINIMAL
            )
            for r in results:
                csvW.writerow(r)

if __name__ == '__main__':
    main()
            