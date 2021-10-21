from time import localtime, strftime, gmtime
from subprocess import call, run
import datetime

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
    print('>>> {}:\t{}\n'.format(
        datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S:%f"),
        stringFormat
    ))



def caller(*args, **kwargs):
    printer(f"| Calling: {args}")
    call(*args, **kwargs)
    printer(f"| Finished")

def runner(*args, **kwargs):
    printer(f"| Calling: {args}")
    run(*args, **kwargs)
    printer(f"| Finished")
