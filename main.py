import utils
from parsers import ClustalParser, MuscleParser
from matcher import AligmentDifferenceFinder

def runClustal(inputFile, reference_id, nseq = 3):
    """function that contain algorithm to compute mismatches Clustal
    clw-->txt file and return the output name"""
    parser = ClustalParser(nseq) #parser.py #(quante sequenze+file) --> sequenze lette
    alignment = parser.parse(inputFile, reference=reference_id)
    print('Read {} bases'.format(len(alignment)))
    analyzer = AligmentDifferenceFinder()
    groups = analyzer.analyze(alignment)
    return utils.save(alignment, analyzer, path='output', reference_id = reference_id)

def runMuscle(inputFile, reference_id, nseq = 3):
    """function that contain algorithm to compute mismatches Muscle
    clw file and return the output name"""
    muscle_parser = MuscleParser(nseq)  #parser.py #(quante sequenze+file) --> sequenze lette
    muscle_alignment = muscle_parser.parse(inputFile, reference=reference_id)
    print('Read {} bases'.format(len(muscle_alignment)))
    analyzer = AligmentDifferenceFinder() # matcher.py
    groups = analyzer.analyze(muscle_alignment) # vettori indici mismatch rispetto al reference per sequenza
    return utils.save(muscle_alignment, analyzer, path='output', reference_id=reference_id)


def main():
    #get id of reference sequence
    reference_id = open("input/reference.fasta", "r").readline().split(' ')[0][1:] # NC_045512.2

    # Clustal Parser iran
    print(runClustal('analysis/iran-ref.txt', reference_id, 3))
""""BUG: same name output with same sequences and 2 pairs of sequences with different alignments"""
    # Muscle Parser iran
    print(runMuscle('analysis/muscle-I20200523-084930-0610-44096621-p1m.clw', reference_id, 3))

    # Clustal Parser israel
    runClustal('analysis/israel-ref.txt', reference_id, 3)

    # Muscle Parser israel
    runMuscle('analysis/muscle-I20200523-085708-0753-28920419-p1m.clw', reference_id, 4)

    #Clustal Parser ncbi
    runClustal('analysis/all.txt', reference_id, 3)

    # Muscle Parser ncbi
    runMuscle('analysis/muscle-I20200512-170208-0225-69454386-p1m.clw', reference_id, 8)

    # Clustal Global Parser
    runClustal('analysis/global.txt', reference_id, 14)

    # Muscle Global Parser
    runMuscle('analysis/muscle-I20200523-090837-0910-95910164-p1m.clw', reference_id, 15)


if __name__ == "__main__":
    main()
