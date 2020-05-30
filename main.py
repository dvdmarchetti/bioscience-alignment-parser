import utils
from parsers import ClustalParser, MuscleParser
from matcher import AligmentDifferenceFinder
import os #to remove file in folder

def runClustal(inputFile, reference_id, nseq = 3):
    """Call functions to parse a ClustalW alignment and produce a custom output structure"""
    parser = ClustalParser(nseq) #parser.py #(quante sequenze+file) --> sequenze lette
    alignment = parser.parse(inputFile, reference=reference_id)
    analyzer = AligmentDifferenceFinder()
    groups = analyzer.analyze(alignment)

    return utils.save(alignment, analyzer, reference_id=reference_id, tool='Muscle', path='output')


def runMuscle(inputFile, reference_id, nseq = 3):
    """Call functions to parse a MUSCLE alignment and produce a custom output structure"""
    muscle_parser = MuscleParser(nseq)  #parser.py #(quante sequenze+file) --> sequenze lette
    muscle_alignment = muscle_parser.parse(inputFile, reference=reference_id)
    analyzer = AligmentDifferenceFinder() # matcher.py
    groups = analyzer.analyze(muscle_alignment) # vettori indici mismatch rispetto al reference per sequenza

    return utils.save(muscle_alignment, analyzer, reference_id=reference_id, tool='Clustal', path='output')


def cleanUpOutputDir(folder=None):
    """clean folder output to make spaces for new output files"""
    dirname = os.path.dirname(__file__)
    path_to_dir = os.path.join(dirname, folder)  # path to directory you wish to remove
    files_in_dir = os.listdir(path_to_dir)     # get list of files in the directory

    for file in files_in_dir:# loop to delete each file in folder
        if os.name == 'nt': #if Windows
            os.remove(f'{path_to_dir}/{file}')     # delete file
        else: #if linux
            import sh
            sh.rm(f)

    #create new empty file for differences
    file = open(os.path.join(path_to_dir, "differences.txt"),"w+")
    file.write("Differences Clustal-Muscle alignment" + '\n\n')
    file.close()


def main():
    reference_id = open("input/reference.fasta", "r").readline().split(' ')[0][1:] # NC_045512.2

    cleanUpOutputDir(folder='output')

    # # ClustalJ = file json from clustal alignment; MuscleJ = file json from muscle alignment
    # # Iran
    # ClustalJ = runClustal('analysis/iran-ref.txt', reference_id, 3)
    # MuscleJ = runMuscle('analysis/muscle-I20200523-084930-0610-44096621-p1m.clw', reference_id, 3)
    # diff = utils.compare_outputs(clustal='output/'+ ClustalJ, muscle='output/'+ MuscleJ)
    # utils.saveCompareFile("differences.txt", "Iran", diff)

    # # Israel
    # ClustalJ = runClustal('analysis/israel-ref.txt', reference_id, 4)
    # MuscleJ = runMuscle('analysis/muscle-I20200523-085708-0753-28920419-p1m.clw', reference_id, 4)
    # diff = utils.compare_outputs(clustal='output/'+ ClustalJ, muscle='output/'+ MuscleJ)
    # utils.saveCompareFile("differences.txt", "Israel ", diff)

    # # GISAID only
    # ClustalJ = runClustal('analysis/GISAID-all.txt', reference_id, 7)
    # MuscleJ = runMuscle('analysis/muscle-I20200523-090216-0023-41230765-p1m.clw', reference_id, 7)
    # diff = utils.compare_outputs(clustal='output/'+ ClustalJ, muscle='output/'+ MuscleJ)
    # utils.saveCompareFile("differences.txt", "GISAID ", diff)

    # # NCBI only
    # ClustalJ = runClustal('analysis/all.txt', reference_id, 3)
    # MuscleJ = runMuscle('analysis/muscle-I20200512-170208-0225-69454386-p1m.clw', reference_id, 8)
    # diff = utils.compare_outputs(clustal='output/'+ ClustalJ, muscle='output/'+ MuscleJ)
    # utils.saveCompareFile("differences.txt", "ncbi ", diff)

    # # All sequences
    # ClustalJ = runClustal('analysis/global.txt', reference_id, 14)
    # MuscleJ = runMuscle('analysis/muscle-I20200523-090837-0910-95910164-p1m.clw', reference_id, 15)
    # diff = utils.compare_outputs(clustal='output/'+ ClustalJ, muscle='output/'+ MuscleJ)
    # utils.saveCompareFile("differences.txt","Global ", diff)

    # Compare global alignments
    print('[...]  Parsing ClustalW analysis...')
    ClustalJ = runClustal('analysis/clustal-global.clustal_num', reference_id, 14)
    print('[DONE] Parsing ClustalW analysis')

    print('[...]  Parsing Muscle analysis...')
    MuscleJ = runMuscle('analysis/muscle-global.clw', reference_id, 14)
    print('[DONE] Parsing Muscle analysis ')

    print('[...]  Comparing alignments...')
    diff = utils.compare_outputs(clustal='output/'+ ClustalJ, muscle='output/'+ MuscleJ, path='output')
    print('[DONE] Comparing alignments')
    # utils.saveCompareFile("differences.txt","Global ", diff)


if __name__ == "__main__":
    main()
