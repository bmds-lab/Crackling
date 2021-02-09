CODE_ACCEPTED = 1
CODE_REJECTED = 0
CODE_UNTESTED = "?"
CODE_AMBIGUOUS = "-"
CODE_ERROR = "!"

MODULE_MM10DB = 'mm10db'
MODULE_SGRNASCORER2 = 'sgrnascorer2'
MODULE_CHOPCHOP = 'chopchop'
MODULE_CONSENSUS = 'consensus'
MODULE_SPECIFICITY = 'specificity'

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
