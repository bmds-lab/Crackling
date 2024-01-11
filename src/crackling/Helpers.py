from subprocess import run
from datetime import datetime

__all__ = ['rc','transToDNA','AT_percentage','printer','runner','elapsedTimeString']

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


# Function that formats provided text with time stamp
def printer(stringFormat):
    print('>>> {}:\t{}\n'.format(
        datetime.now().strftime("%Y-%m-%d %H:%M:%S:%f"),
        stringFormat
    ))


# Function that runs given external call and records start and finish times using printer
def runner(*args, **kwargs):
    printer(f"| Calling: {args}")
    run(*args, **kwargs)
    printer(f"| Finished")


# Function that generates a string representing elapsed time
def elapsedTimeString(start_time, end_time):
    elapsed_seconds = end_time - start_time

    days, remainder = divmod(elapsed_seconds, 86400)
    hours, remainder = divmod(remainder, 3600)
    minutes, seconds = divmod(remainder, 60)

    return f"{int(days)} {int(hours):02d}:{int(minutes):02d}:{int(seconds):02d}"
