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
        'analyzed_sequences': sequences,
        'unmatches': unmatching_groups
    }
#datetime.timestamp(now).toString() in order to permit 2 outputs of same sequences in multiple tools
    seq_hash = hashlib.sha1(bytes(str(datetime.timestamp(now)).join(sequences), encoding='utf8'))
    filename = '{}_{}.json'.format(alignment.reference, seq_hash.hexdigest())
    with open(os.path.join(path, filename), 'w') as output:
        output.write(json.dumps(payload, indent=True))

    return filename

def jsonComp(file1, file2):
    """return differences file json output"""
    with open(file1) as f1:
        data1 = json.load(f1)

    with open(file2) as f2:
        data2 = json.load(f2)

    different_items = {}#{k: data1[k] for k in data1 if k in data2 and data1[k] != data2[k]}

    value = []
    for key in data1.keys():
        if (data1.get(key, None) != data2.get(key, None)):
            value.append(key)

#salta sta parte e parti dagli allineamenti
    """if 'reference' in value:
        return "insert wrong reference alignment WTF?"

    #if i get error for used different sequences in the alignment
    if 'analyzed_sequences' in value:
        list  = data1.get('analyzed_sequences', None) + data2.get('analyzed_sequences', None) #sum for then compare
        extras = [value for value in list if (value in data1.get('analyzed_sequences', None)) ^ (value in data2.get('analyzed_sequences', None))]
        return "different sequences used: {}".format(extras) """

#real differences in alignments
    lFrom,lTo = [], []
    if 'unmatches' in value: #TODO
        for x in data1.get('unmatches').keys():
            #print(data1.get('unmatches').get(x))
            lFrom.append(data1.get('unmatches').get(x).get('from'))
            lTo.append(data1.get('unmatches').get(x).get('to'))

        print(lFrom)
        for y in data2.get('unmatches').keys():
            #if not data2.get('unmatches').get(y).get('from') in lFrom:
                #different_items[y] = [data2.get('unmatches').get(y)]
            if data2.get('unmatches').get(y).get('from') in lFrom:
                pos = lFrom.index(data2.get('unmatches').get(y).get('from'))
                #print(pos)
                if data2.get('unmatches').get(y).get('to') != lTo[pos]:
                    """TO FIX"""
                    different_items[y] = [data1.get('unmatches').get(x), data2.get('unmatches').get(y)]
                else:
                    #remove matching_groups
                    del lFrom[pos]
                    del lTo[pos]
            else:
                different_items[y] = [data2.get('unmatches').get(y)]

        #insert remaining


        """print(len(unmatches))
        cont = 0;
        for y in data2.get('unmatches').keys():
            #print(data2.get('unmatches').get(y))
            if ((data2.get('unmatches').get(y).get('from') != unmatches[cont].get('from')) or
             (data2.get('unmatches').get(y).get('to') != unmatches[cont].get('to'))):
                print(cont)
                different_items[y] = [unmatches[cont], data2.get('unmatches').get(y)]
                #print(type(unmatches[cont]))
            cont = cont + 1
        #for v in value:
            #print(type(v))
            #print(type(data1.get(v, None)))
            #print("data2")
            #print(type(data2.get(v, None)))
            """

    return different_items
