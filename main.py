# Differenze tra reference e altri allineamenti
# conta le basi distinte
# python ad-hoc (secondo anno)
# formato: inizio, lunghezza (compatta mismatch consecutivi), tipo di match tra ref e sequenza allineata, dove per sequenza singola
#usa GISAID (come si scrive?), specifica tessuto da cui è prelevato il virus + età paziente;
#almeno 16 sequenze; specifica e confronta i tools; analizza albero in output (con date)
#NCBI gene start..end ha cds da tradurre; vedere se esistono mismatch tra frequenze nei punti in cui esistono geni
#annotare dove e in che gene(se in gene) esistono mutazioni
#sostituire simbolo != ACGT con parte nella ref [vedi IUPAC code]
#Tools: Muscle, Clustal; identificare origine:
#NCBI
import re

import utils
from parsers import ClustalParser
from matcher import AligmentDifferenceFinder
import json #for JSON


reference_id = open("input/reference.fasta", "r").readline().split(' ')[0][1:] # NC_045512.2
#analysis = open("analysis/iran-ref.txt", "r")

aligned_sequences = 14 #numero frequenze allineate (costante globale)
sequences=['MT276598.1', 'MT320891.2'] #

parser = ClustalParser(nseq=aligned_sequences)  #parser.py #(quante sequenze+file) --> sequenze lette
alignment = parser.parse('analysis/Global.txt', reference = reference_id, list = []) #the other is all.txt
#alignment = parser.parse('analysis/iran-ref.txt', reference='NC_045512.2')  #alignment == Ref [12]

print('Read {} bases'.format(len(alignment)))

#return Array con tutti i singoli indici con mismatch nelle stringhe.
diff_finder = AligmentDifferenceFinder()    #matcher.py
groups = diff_finder.analyze(alignment)     #vettori indici mismatch almeno 2 sequenze
#print(groups)

print('Diff ranges: {}'.format(utils.group_ranges(groups)))
#print(utils.group_ranges(groups)) ###utils.group_ranges(groups) = gruppi mismatch contigui
#ritornare i range dei singoli gruppi di numeri contigui
with open('data.json', 'w') as f:
    for group in utils.group_ranges(groups):
        start, end = group  #inizio e fine mismatch

        print('> Range: {}-{}'.format(start, end))
        print('Reference: {}'.format(alignment.peek_reference(start, end)))
        print('Others: {}'.format(alignment.peek_others(start, end)))
    # print(alignment[start:end+1])
    # print(alignemnt.reference[])
        print()
        #ids = [reference_id, alignment.peek_others(start, end)[1]['sequence']]
        #print(type(alignment.peek_others(start, end)[1]))
        #print(alignment.peek_others(start, end))
        #print(len(alignment.peek_others(start, end)))
        ids = [reference_id]
        for i in range(len(alignment.peek_others(start, end))): #i = classical counter
            if(not alignment.peek_others(start, end)[i]['value'] == alignment.peek_reference(start, end)):
                ids.append(alignment.peek_others(start, end)[i]['sequence'])

        data = {'from': start,
            'to': end,
            'sequences': ids, #None
            }

        json.dump(data, f)
        if group != utils.group_ranges(groups)[-1]: #do not add extra \n at the end
            f.write('\n')   #\n after every entry (add extra empty line)
### END FOR
f.close()
### In questo modo hai un Array nel quale ogni elemento è a sua volta un Array che corrisponde agli indici di inizio e fine della zona che non matcha con la reference
"""possibile futuro lavoro: stampare mismatch solo ref e sequenza con mismatch
#parsers.parseLines() only if id = ....."""
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
