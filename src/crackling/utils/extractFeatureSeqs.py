"""
This utility extracts sequences from FASTA files as according a GFF3 annotation.

Author: Jake Bradford

This can be imported into a Python script. Use `extractSequences(..)`.

Otherwise, run from CLI:
    py -3 extractSequences -g annotation.gff -f mouse_*.fa -o output_directory

    There exist more CLI arguments. To see help documentation, run:
        py -3 extractSequences --help
"""

from Bio import SeqIO
import os

class GffRecord:
    def __init__( 
        self, seqId, source, type, start, end, score, strand, phase, attributes
    ):
        """
        A convenient object to store a GFF record in.
        See the GFF3 specification here: https://m.ensembl.org/info/website/upload/gff3.html
        
        The attributes of the GFF record must contain an `ID` field.
        """
        
        self.seqId = seqId
        self.source = source
        self.type = type
        self.start = int(start)
        self.end = int(end)
        self.score = score
        self.strand = strand
        self.phase = phase
        self.attributes = self.strAttributesToDict(attributes)
    
    def __str__(self):
        attrId = self.attributes.get('ID', None)
        return f'<GffRecord seqId="{self.seqId}" type="{self.type}" attrId="{attrId}">'
    
    def strAttributesToDict(self, attrs):
        dictAttrs = {
            attr.split('=')[0] : attr.split('=')[1]
            for attr in attrs.split(';')
        }
        if 'ID' not in dictAttrs:
            raise RunTime('GFF record contains no attribute `ID`')
        return dictAttrs

class GffParser:
    def __init__(self, filePath, sampleSize=None):
        """
        Loads in a GFF3-formatted annotation.
        
        Arguments:
            filePath: the path to the GFF3 file
            sampleSize (optional): None indicates the entire file will be loaded.
                If you wish to load only a portion of the file, specify the
                number of lines that should be processed (ignoring GFF comments).
        """
        self.records = []
        
        with open(filePath, 'r') as fpGff:
            for line in fpGff:
                line = line.strip()
                
                if len(line) > 0 and line[0] != '#':
                
                    if sampleSize is not None:
                        if sampleSize == 0:
                            break
                        else:
                            sampleSize -= 1
                            
                    self.records.append(
                        GffRecord(*line.split('\t'))
                    )
                    
    def __str__(self):
        return f"<GffParser size={len(self.records)}>"
    
    def getRecordsByTypes(self, types):
        for record in self.records: 
            if record.type in types:
                yield record
        
def _writeFastaRecord(fpHandle, title, seq, wrap=0):
    """
    This is an internal method for writing a FASTA record to a file.
    
    Arguments:
        fpHandle (file object in write mode)
        
        title (string): to be written as `>{title}`
        
        seq (string): to be written on a single line or wrapped (according to 
            `wrap`)
            
        wrap (int): the maximum line length for the sequence. 
            Setting to zero indicates no wrap.
    """
    fpHandle.write(f">{title}\n")
    if wrap == 0:
        fpHandle.write(f"{seq}\n")
    else:
        for i in range(0, len(seq), wrap):
            fpHandle.write(f"{seq[i : i + wrap]}\n")

def extractSequences(
    listFastaPaths,
    gffPath, 
    outputDirectory,
    featureLevels=['exon'],
    outputOneFilePerSequence=False,
    outputOneFileName='sequences.fa',
    outputSequenceWrap=60,
    upstreamBuffer=0,
    downstreamBuffer=0,
    verbose=False
):
    """
    This function performs the extraction
    
    Arguments:
        listFastaPaths (list): paths to input FASTA files
        
        gffPath (string): a path to the GFF3-formatted annotation
        
        outputDirectory (string): a directory to write the output files to.
            If not does not exist, it will be created.
            
        featureLevels (list): which features should be extracted
        
        outputOneFilePerSequence (bool): 
            True one file is created. 
            False one file per feature is created.
            Default: False
            
        outputOneFileName (string):
            If output is written to a single file then provide the name of this
            file here.
            
        outputSequenceWrap (int):
            The number of characters per line in the output file.
            Setting to zero indicates the sequence should be on one line.
            Default: 60
        
        upstreamBuffer (int):
            The number of positions to expand the sequence at the 5' end.
        
        downstreamBuffer (int):
            The number of positions to expand the sequence at the 3' end.
    """
    
    # Make sure the output directory exists
    if not os.path.exists(outputDirectory):
        if verbose:
            print(f'Creating directory: {outputDirectory}')
        os.mkdir(outputDirectory)
    else:
        if verbose:
            print(f'Output directory already exists: {outputDirectory}')
            
    # A dictionary containing all input sequences
    # Key: FASTA header
    # Value: sequence
    seqs = {}
    
    # Load the input sequences from the FASTA files
    if verbose:
        print('Beginning to load FASTA sequences into memory')
        
    for fastaPath in listFastaPaths:
        with open(fastaPath, 'r') as fpFasta:
            for record in SeqIO.parse(fpFasta, "fasta"):
                if verbose:
                    print(f"Loading: {fastaPath}")
            
                # The ID returned by Bio.SeqIO.FastaIO is 
                #  `title.split(None, 1)[0]` where `title` is the FASTA header.
                # To access the entire FASTA header, use `record.description`.
                seqs[record.id] = record.seq
    
    if verbose:
        print(f"Loaded {len(seqs)} sequences")
        
    # Load the GFF annotation
    if verbose:
        print(f"Beginning to load GFF annotation into memory: {gffPath}")
    gff = GffParser(gffPath)

    if outputOneFilePerSequence:
        # Generate one file per output sequence
        if verbose:
            print('Generating one file per output sequence')
            
        # The ID of the feature will be used as the filename
        for feature in gff.getRecordsByTypes(featureLevels):

            # Check if this annotation record has an accompanying sequence
            if feature.seqId not in seqs:
                continue
            
            # Open a new file for this annotation record
            path = os.path.join(
                outputDirectory,
                f"{feature.attributes['ID']}.fa"
            )
            
            with open(path, 'w') as fpWrite:
                if verbose:
                    print(f"Writing to: {path}")
            
                # Write the title
                title = feature.attributes['ID']

                # Write the sequence (either on one line or wrap it)
                seq = seqs[feature.seqId][
                    feature.start - 1 - upstreamBuffer:
                    feature.end + downstreamBuffer
                ]
                
                _writeFastaRecord(fpWrite, title, seq, wrap=outputSequenceWrap)
    else:
        # Generate one file for all output sequences
        if verbose:
            print('Generating one file for all output sequences')
            
        path = os.path.join(
            outputDirectory,
            outputOneFileName
        )
        
        with open(path, 'w') as fpWrite:
            if verbose:
                print(f"Writing to: {path}")
                
            for feature in gff.getRecordsByTypes(featureLevels):
            
                # Check if this annotation record has an accompanying sequence
                if feature.seqId not in seqs:
                    continue

                # Write the title
                title = feature.attributes['ID']

                # Write the sequence (either on one line or wrap it)
                seq = seqs[feature.seqId][
                    feature.start - 1 - upstreamBuffer:
                    feature.end + downstreamBuffer
                ]
                
                _writeFastaRecord(fpWrite, title, seq, wrap=outputSequenceWrap)

    if verbose:
        print('Done')

def main():
    import argparse
    from glob import glob

    parser = argparse.ArgumentParser(
        description='Extract sequences from FASTA files according to GFF annotation records.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('-f', '--fasta',
        action='append', nargs='+', required=True,
        help='Paths to the FASTA files to process. You may use wildcards.'
    )
    
    parser.add_argument('-g', '--gff',
        required=True,
        help='The path to the GFF file.'
    )
    
    parser.add_argument('-o', '--output-dir',
        required=True,
        help='The output directory to write in.'
    )
    
    parser.add_argument('-l', '--feature-levels',
        default='exon', required=False,
        help='A comma-separated list of features (i.e. the third column) to use from the annotation.'
    )
    
    parser.add_argument('--one-file-per-seq',
        action='store_true', required=False,
        help='Write each output sequence to its own file.'
    )
    
    parser.add_argument('--output-filename',
        required=False, default='sequences.fa',
        help='If writing all sequences to one output file, specify the name of the file here.'
    )
    
    parser.add_argument('--seq-line-length',
        required=False, type=int, default=60,
        help='The line length for sequences written to file.'
    )
    
    parser.add_argument('--upstream-buffer',
        required=False, type=int, default=0,
        help="The number of positions to expand the sequence at the 5' end."
    )
    
    parser.add_argument('--downstream-buffer',
        required=False, type=int, default=0,
        help="The number of positions to expand the sequence at the 3' end."
    )
    parser.add_argument('-v', '--verbose',
        action='store_true', required=False,
        help='Print what I am doing.'
    )

    args = parser.parse_args()

    if args.verbose:
        print("Arguments:")
        for key, value in parser.parse_args()._get_kwargs():
            print(f"\t{key} = {value}")

    extractSequences(
        args.fasta[0],
        args.gff,
        args.output_dir,
        featureLevels = args.feature_levels.split(','),
        outputOneFilePerSequence = args.one_file_per_seq,
        outputOneFileName = args.output_filename,
        outputSequenceWrap = args.seq_line_length,
        upstreamBuffer = args.upstream_buffer,
        downstreamBuffer = args.downstream_buffer,
        verbose = args.verbose
    )

if __name__ == '__main__':
    main()
