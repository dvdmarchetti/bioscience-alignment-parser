import hashlib
import json
import os

def chunk(iterable, size):
    """Yield successive n-sized chunks from iterable."""
    for i in range(0, len(iterable), size):
        yield iterable[i:i + size]


def group_ranges(values):
    """
    Chunk the input list into groups of consecutive elements.

    Each group is in form of [start_index, end_index)
    """
    units = []
    prev = values[0]

    for value in values:
        if value == prev + 1:
            units[-1].append(value)
            # None
        else:
            units.append([value])
        prev = value

    return [[u[0], u[-1] + 1] if len(u) > 1 else [u[0], u[0] + 1] for u in units]

def removeLn(file):
    import os, sys
    """
    remove last line file
    """
    readFile = open(file)
    lines = readFile.readlines()
    readFile.close()
    w = open(file,'w')
    w.writelines([item for item in lines[:-1]])
    w.close()

def save(alignment, analyzer, path=None, reference_id='NC_045512.2'):
    if not analyzer.unmatches:
        raise AssertionError('You must run an analysis to save it results.')

    if not path:
        path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

    unmatching_groups = {}
    for sequence, groups in analyzer.unmatches.items():
        for start, end in groups:
            key = hash((start, end))
            if key in unmatching_groups:
                unmatching_groups[key]['sequences'].append(sequence)
            else:
                unmatching_groups[key] = {
                    'from': start,
                    'to': end,
                    'sequences': [sequence]
                }

    sequences = tuple(analyzer.unmatches.keys())
    payload = {
        'reference': reference_id,
        'analyzed_sequences': sequences,
        'unmatches': unmatching_groups
    }

    seq_hash = hashlib.sha1(bytes('-'.join(sequences), encoding='utf8'))
    filename = '{}_{}.json'.format(alignment.reference, seq_hash.hexdigest())
    with open(os.path.join(path, filename), 'w') as output:
        output.write(json.dumps(payload, indent=True))
