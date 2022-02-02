import re
import os
import joblib

from crackling.Helpers import printer


COMPLIMENTS = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')

def reverse_complement(sequence):
    """
        Returns the reverse-complement of a given sequence
    """
    return sequence.translate(COMPLIMENTS)[::-1]


def process_sequence(sequence, sequence_header):
    # Patterns for guide matching
    pattern_forward = r'(?=([ATCG]{21}GG))'
    pattern_reverse = r'(?=(CC[ACGT]{21}))'

    # New sequence deteced, process sequence
    # once for forward, once for reverse
    for pattern, strand, sequence_modifier in [
        [pattern_forward, '+', lambda x : x], 
        [pattern_reverse, '-', lambda x : reverse_complement(x)]
    ]:
        p = re.compile(pattern)
        for m in p.finditer(sequence):
            target23 = sequence_modifier(sequence[m.start() : m.start() + 23])
            yield [target23, sequence_header,  m.start(),  m.start() + 23, strand]


def find_guides(sequence_header, sequence):
    guides = {}

    for guide in process_sequence(sequence, sequence_header):
        if guide[0] in guides.keys():
            guides[guide[0]] += 1
        else:
            guides[guide[0]] = 1

    return (sequence_header, guides)


def load_fasta_sequence_file(filename):
    sequence_header = None
    sequence = []

    with open(filename, 'r') as file:
        for line in file:
            # Clean line
            line = line.strip()

            # Sequence Header
            if (line[0] == '>'):
                if (sequence_header is not None):
                    yield sequence_header, ''.join(sequence)

                sequence_header = line[1:]
                sequence = []

            # Regular line
            else:
                sequence.append(line)
        
        # Return last line
        yield sequence_header, ''.join(sequence)


def find_candidates_in_file(guide_batchinator, target_file, candidate_guides, duplicate_guides, recorded_sequences):
    assert isinstance(candidate_guides, set)
    assert isinstance(duplicate_guides, set)
    assert isinstance(recorded_sequences, set)

    duplicate_guide_count = 0
    identified_guide_count = 0

    target_file_size = os.path.getsize(target_file)

    printer(f'Identifying possible target sites in: {target_file}')
    results = joblib.Parallel(n_jobs=-1)(joblib.delayed(find_guides)(sequence_header, sequence) for sequence_header, sequence in load_fasta_sequence_file(target_file))

    printer(f'Combining results from {len(results)} sequence headers')
    # Combine Results
    for (sequence_header, guides) in results:
        recorded_sequences.add(sequence_header)

        for guide, guide_count in guides.items():
            identified_guide_count += 1
            if (guide_count > 1) or (guide in candidate_guides):
                duplicate_guides.add(guide)
                duplicate_guide_count += 1
            else:
                guide_batchinator.recordEntry(guide)

            candidate_guides.add(guide)

    return candidate_guides, duplicate_guides, recorded_sequences, target_file_size, duplicate_guide_count, identified_guide_count
