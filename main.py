# Differenze tra reference e altri allineamenti
# conta le basi distinte
# python ad-hoc (secondo anno)
# formato: inizio, lunghezza (compatta mismatch consecutivi), tipo di match tra ref e sequenza allineata, dove per sequenza singola
import re

import utils
from parsers import ClustalParser
from matcher import AligmentDifferenceFinder


# reference_id = open("input/reference.fasta", "r").readline().split(' ')[0][1:] # NC_045512.2
# analysis = open("analysis/iran-ref.txt", "r")

aligned_sequences = 3

parser = ClustalParser(nseq=aligned_sequences)
alignment = parser.parse('analysis/iran-ref.txt', reference='NC_045512.2')

print('Read {} bases'.format(len(alignment)))

diff_finder = AligmentDifferenceFinder()
groups = diff_finder.analyze(alignment)

print('Diff ranges: {}'.format(utils.group_ranges(groups)))

for group in utils.group_ranges(groups):
    start, end = group

    print('> Range: {}-{}'.format(start, end))
    print('Reference: {}'.format(alignment.peek_reference(start, end)))
    print('Others: {}'.format(alignment.peek_others(start, end)))
    # print(alignment[start:end+1])
    # print(alignemnt.reference[])
    print()


# blocks = []
# for line in chunk(lines, aligned_sequences + 1):
#     blocks.append(ClustalAlignmentBlock.fromRaw(lines))

# print(blocks[0])
# print(blocks[0].count_differences())

# lines = pd.DataFrame(map(clean, lines))

####IDEA: Dzionario chiave = ID line, valore = contatore della differenza con refSeqId
#(tieni al momento anche refSeqId per sicurezza)
# sequences = dict()

# #print(line[0].split(' ')[0])   DEVO RIMUOVER
# #print(line[0].split(' ')[1])   elementi =
# #print(line[2].split(' ')[0])   ['']
# #print(line[2].split(' ')[1])
# cont = 0
# for l in line:
#     if (l.split(' ')[0] != ''):
#         #key = l.split(' ')[0]
#         key, seq = itemgetter(0, 1)(l.split(' '))
#         if key not in sequences.keys():
#             sequences[key] = 0
# ###take seq and compare with sequences[refSeqId], sequences[key] = sequences[key] ++ fpr each mismatch
#         refSeq = findSequence(l, refSeqId, cont)
#         print(seq)
#         for i in range(seq):
#             if seq[i] != refSeq[i]:     #give unknown problems
#                 sequences[key] = sequences[key] +1

#     cont = cont + 1 # why no ++ operand in python :(
#     #print(l.split(' '))
#     #print(l.split(' ')[1])#error, remove spaces doesn't work
# ###TODO contare mismetch rispetto refSeqId

# print(sequences)
