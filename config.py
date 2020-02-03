CONFIG = {

    # Run name
    # If not set, then the start time is used.
    'name' : 'sample',

    # The consensus approach
    'consensus' : {
    
        # The number of methods, at least, that must agree
        'n' : 2
    },

    # Pipeline input
    'input' : {

        # Either a directory which contains one file per sequence, or
        # a file which contains line-separated sequences
        'exon-sequences' : r'/sample/scaffolds/',
        
        # A line-separated values file containing a list of all off-target sites
        'offtarget-sites' : r'/sample/offtargetSites.txt',
        
        # GFF formatted annotation of the genome
        'gff-annotation' : r'/sample.gff',
        
        # Bowtie2 index
        'bowtie2-index' : r'/sample/sample.fa',
    },
    
    # Pipeline output
    'output' : {
    
        # Directory to write output to
        'dir' : r'./sample-output/',
        
        # File to write guides to file
        'fileName' : r'guides.txt',

        # Output delimiter to separate values
        'delimiter' : '\t',
    },
    
    # Off-target scoring method config
    'offtargetscore' : {
        'binary' : r'./search_ots_score',
        'threads' : 128,
        'score-threshold' : 75, # if score < score-threshold then reject
        
    },
    
    # sgRNAScorer2 config
    'sgrnascorer2' : {
        'model' : r'model-py3.txt',
        'score-threshold' : 0, # if score < score-threshold then reject
    },
    
    # Bowtie2 config for specificity evaluation
    'bowtie2' : {
        'binary' : r'bowtie2',
        'threads' : 128,
    },
    
    # Secondary structure calculation by RNAfold
    'rnafold' : {
        'binary' : r'/usr/local/bin/RNAfold',
        'threads' : 128,
        'low_energy_threshold' : -30,
        'high_energy_threshold' : -18
    }
}