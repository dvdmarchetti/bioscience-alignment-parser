import hashlib
import json
import os
from datetime import datetime

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
    import sys
    """
    remove last line file
    """
    readFile = open(file)
    lines = readFile.readlines()
    readFile.close()
    w = open(file,'w')
    w.writelines([item for item in lines[:-1]])
    w.close()

def save(alignment, analyzer, reference_id=None, tool=None, path=None):
    """
    save results of all mismatches in a file.json named: reference_hash(Sequences+timestamp).json
    and return the name"""
    if not analyzer.unmatches:
        raise AssertionError('You must run an analysis to save it results.')

    if not path:
        path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

    now = datetime.now()
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
        'tool': tool,
        'analyzed_sequences': sequences,
        'unmatches': unmatching_groups
    }

    #datetime.timestamp(now).toString() in order to permit 2 outputs of same sequences in multiple tools
    # seq_hash = hashlib.sha1(bytes(str(datetime.timestamp(now)).join(sequences), encoding='utf8'))
    # filename = '{}_{}.json'.format(alignment.reference, seq_hash.hexdigest())
    filename = '{}_{}.json'.format(alignment.reference, str(datetime.timestamp(now)))
    with open(os.path.join(path, filename), 'w') as output:
        output.write(json.dumps(payload, indent=True))

    return filename

def jsonComp(file1, file2):
    """return differences file json output"""
    with open(file1) as f1:
        left = json.load(f1)

    with open(file2) as f2:
        right = json.load(f2)

    # different_items = []#{k: left[k] for k in left if k in right and left[k] != right[k]}
    diff_left = []
    diff_right = []

    value = []
    for key in left.keys():
        if (left.get(key, None) != right.get(key, None)):
            value.append(key)

    if left.get('unmatches') is None or right.get('unmatches') is None:
        print('Invalid output file provided')
        return [[], []]

    leftKeys = set(left.get('unmatches').keys())
    rightKeys = set(right.get('unmatches').keys())

    # Symmetric diff: take elements that are either in xkeys or ykeys, but not in both
    for key in leftKeys ^ rightKeys:
        if left.get('unmatches').get(key) is not None:
            diff_left.append(left.get('unmatches').get(key))
        if right.get('unmatches').get(key) is not None:
            diff_right.append(right.get('unmatches').get(key))


    diff = []
    cont = 0
    while cont < len(diff_left) and cont < len(diff_right):
        if cont > len(diff_left):
            diff.append([None, diff_right[cont]])
        else:
            if cont > len(diff_right):
                diff.append([diff_left[cont], None])
            else:
                diff.append([diff_left[cont], diff_right[cont]])
        cont = cont + 1

    return [diff]

def saveCompareFile(filename = "differences.txt", country = "", diff = [], path='output'):
    if not path:
        path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

    #append new
    file = open(os.path.join(path, filename),"a+")
    file.write(country + '\n')
    for i in range(0, len(diff)):
        file.write(str(diff[i]) + '\n')

    file.write('\n')
    file.close()
