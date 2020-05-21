# Differenze tra reference e altri allineamenti
# conta le basi distinte
# python ad-hoc (secondo anno)
# formato: inizio, lunghezza (compatta mismatch consecutivi), tipo di match tra ref e sequenza allineata, dove per sequenza singola
import re

import utils
from parsers import ClustalParser, MuscleParser
from matcher import AligmentDifferenceFinder


reference_id = open("input/reference.fasta", "r").readline().split(' ')[0][1:] # NC_045512.2

# aligned_sequences = 8 #numero frequenze allineate (costante globale)

muscle_parser = MuscleParser(nseq=8)  #parser.py #(quante sequenze+file) --> sequenze lette
muscle_alignment = muscle_parser.parse('analysis/muscle-I20200512-170208-0225-69454386-p1m.clw', reference=reference_id)

parser = ClustalParser(nseq=3)  #parser.py #(quante sequenze+file) --> sequenze lette
alignment = parser.parse('analysis/iran-ref.txt', reference=reference_id)
# alignment = parser.parse('analysis/iran-ref.txt', reference='NC_045512.2')  #alignment == Ref [12]

print('Read {} bases'.format(len(alignment)))

# return Array con tutti i singoli indici con mismatch nelle stringhe.
diff_finder = AligmentDifferenceFinder() # matcher.py
groups = diff_finder.analyze(alignment) # vettori indici mismatch almeno 2 sequenze
# print(groups)

print('Diff ranges: {}'.format(utils.group_ranges(groups)))
# print(utils.group_ranges(groups)) # utils.group_ranges(groups) = gruppi mismatch contigui

# ritornare i range dei singoli gruppi di numeri contigui
for group in utils.group_ranges(groups):
    start, end = group # inizio e fine mismatch

    print('> Range: {}-{}'.format(start, end))
    print('Reference: {}'.format(alignment.peek_reference(start, end)))
    print('Others: {}'.format(alignment.peek_others(start, end)))
    # print(alignment[start:end+1])
    # print(alignemnt.reference[])
    print()

# In questo modo hai un Array nel quale ogni elemento è a sua volta un Array che corrisponde
# agli indici di inizio e fine della zona che non matcha con la reference

# possibile futuro lavoro: stampare mismatch solo ref e sequenza con mismatch
