import hashlib
import json
import os
import time

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


def save(alignment, analyzer, reference_id=None, tool=None, path=None):
    """
    save results of all mismatches in a file.json named: reference_hash(Sequences+timestamp).json
    and return the name"""
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
                    'reference': alignment.sequences[alignment.reference][start:end],
                    'alt': alignment.sequences[sequence][start:end],
                    'sequences': [ sequence ]
                }

    now = time.localtime()
    sequences = tuple(analyzer.unmatches.keys())
    payload = {
        'tool': tool,
        'timestamp': time.strftime('%Y/%m/%d %H:%M:%S UTC%z', now),
        'reference': reference_id,
        'analyzed_sequences': sequences,
        'unmatches': unmatching_groups
    }

    #datetime.timestamp(now).toString() in order to permit 2 outputs of same sequences in multiple tools
    # seq_hash = hashlib.sha1(bytes(timestamp.join(sequences), encoding='utf8'))
    # filename = '{}_{}.json'.format(alignment.reference, seq_hash.hexdigest())
    filename = '{}-{}_{}.json'.format(tool, alignment.reference, time.strftime('%Y-%m-%d_%H-%M', now))
    with open(os.path.join(path, filename), 'w') as output:
        output.write(json.dumps(payload))

    return filename


def compare_outputs(clustal=None, muscle=None, path=None):
    """return differences file json output"""
    if not clustal or not muscle:
        print('You must provide valid alignment output for both clustal and muscle')
        return

    with open(clustal) as f1:
        left = json.load(f1)

    with open(muscle) as f2:
        right = json.load(f2)

    merge = left
    diff_left = []
    diff_right = []

    value = []
    for key in left.keys():
        if (left.get(key, None) != right.get(key, None)):
            value.append(key)

    if left.get('unmatches') is None or right.get('unmatches') is None:
        print('Invalid output file provided')
        return [], []

    leftKeys = set(left.get('unmatches').keys())
    rightKeys = set(right.get('unmatches').keys())

    # Symmetric diff: take elements that are either in xkeys or ykeys, but not in both
    for key in leftKeys ^ rightKeys:
        if left.get('unmatches').get(key) is not None:
            diff_left.append(left.get('unmatches').get(key))
        if right.get('unmatches').get(key) is not None:
            diff_right.append(right.get('unmatches').get(key))

    tools = [left.get('tool'), right.get('tool')]
    for key in leftKeys & rightKeys:
        item = left.get('unmatches').get(key)
        item['tools'] = tools
        merge['unmatches'][key] = item

    for key in rightKeys - leftKeys:
        item = right.get('unmatches').get(key)
        item['tools'] = [right.get('tool')]
        merge['unmatches'][key] = item

    for key in leftKeys - rightKeys:
        item = left.get('unmatches').get(key)
        item['tools'] = [left.get('tool')]
        merge['unmatches'][key] = item

    now = time.localtime()
    merge['tools'] = tools
    merge.pop('tool')
    merge['timestamp'] = time.strftime('%Y/%m/%d %H:%M:%S UTC%z', now)

    filename = '{}-{}_{}.json'.format('_'.join(tools), merge['reference'], time.strftime('%Y-%m-%d_%H-%M', now))
    with open(os.path.join(path, filename), 'w') as output:
        output.write(json.dumps(merge, indent=True))

    return diff_left, diff_right


def saveCompareFile(filename="differences.txt", country="", diff=[], path='output'):
    if not path:
        path = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

    #append new
    file = open(os.path.join(path, filename),"a+")
    file.write(country + '\n')
    for i in range(0, len(diff)):
        file.write(str(diff[i]) + '\n')

    file.write('\n')
    file.close()
