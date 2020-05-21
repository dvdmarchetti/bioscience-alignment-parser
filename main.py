import hashlib
import json
import os


import utils
from parsers import ClustalParser, MuscleParser
from matcher import AligmentDifferenceFinder


def save(alignment, analyzer, path=None):
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


reference_id = open("input/reference.fasta", "r").readline().split(' ')[0][1:] # NC_045512.2

# Muscle Parser
muscle_parser = MuscleParser(nseq=8)  #parser.py #(quante sequenze+file) --> sequenze lette
muscle_alignment = muscle_parser.parse('analysis/muscle-I20200512-170208-0225-69454386-p1m.clw', reference=reference_id)

print('Read {} bases'.format(len(muscle_alignment)))
analyzer = AligmentDifferenceFinder() # matcher.py
groups = analyzer.analyze(muscle_alignment) # vettori indici mismatch rispetto al reference per sequenza
save(muscle_alignment, analyzer, path='output')


# Clustal Parser
parser = ClustalParser(nseq=3) #parser.py #(quante sequenze+file) --> sequenze lette
alignment = parser.parse('analysis/iran-ref.txt', reference=reference_id)
# alignment = parser.parse('analysis/israel-ref.txt', reference=reference_id)

print('Read {} bases'.format(len(muscle_alignment)))
analyzer = AligmentDifferenceFinder()
groups = analyzer.analyze(alignment)
save(alignment, analyzer, path='output')

# Clustal Global Parser
parser = ClustalParser(nseq=14)
alignment = parser.parse('analysis/global.txt', reference=reference_id)

print('Read {} bases'.format(len(alignment)))
analyzer = AligmentDifferenceFinder()
groups = analyzer.analyze(alignment)
save(alignment, analyzer, path='output')
