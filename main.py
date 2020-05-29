import utils
from parsers import ClustalParser, MuscleParser
from matcher import AligmentDifferenceFinder
import os#to remove file in folder

def runClustal(inputFile, reference_id, nseq = 3):
    """function that contain algorithm to compute mismatches Clustal
    clw-->txt file and return the output name"""
    parser = ClustalParser(nseq) #parser.py #(quante sequenze+file) --> sequenze lette
    alignment = parser.parse(inputFile, reference=reference_id)
    #print('Read {} bases'.format(len(alignment)))
    analyzer = AligmentDifferenceFinder()
    groups = analyzer.analyze(alignment)
    return utils.save(alignment, analyzer, reference_id=reference_id, tool='Muscle', path='output')

def runMuscle(inputFile, reference_id, nseq = 3):
    """function that contain algorithm to compute mismatches Muscle
    clw file and return the output name"""
    muscle_parser = MuscleParser(nseq)  #parser.py #(quante sequenze+file) --> sequenze lette
    muscle_alignment = muscle_parser.parse(inputFile, reference=reference_id)
    #print('Read {} bases'.format(len(muscle_alignment)))
    analyzer = AligmentDifferenceFinder() # matcher.py
    groups = analyzer.analyze(muscle_alignment) # vettori indici mismatch rispetto al reference per sequenza
    return utils.save(muscle_alignment, analyzer, reference_id=reference_id, tool='Clustal', path='output')


def main():
    #get id of reference sequence
    reference_id = open("input/reference.fasta", "r").readline().split(' ')[0][1:] # NC_045512.2

    """clean folder output to make spaces for new output files"""
    dirname = os.path.dirname(__file__)
    path_to_dir = os.path.join(dirname, 'output')  # path to directory you wish to remove
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
    # #ClustalJ = file jason from clustal alignment; MuscleJ = file jason from muscle alignment
    # # Clustal Parser iran
    ClustalJ = runClustal('analysis/iran-ref.txt', reference_id, 3)

    # # Muscle Parser iran
    MuscleJ = runMuscle('analysis/muscle-I20200523-084930-0610-44096621-p1m.clw', reference_id, 3)

    #compare them
    diff = utils.jsonComp('output/'+ ClustalJ, 'output/'+ MuscleJ)
    #print("Iran " + str(diff))
    utils.saveCompareFile("differences.txt", "Iran ", diff)

    # Clustal Parser israel
    ClustalJ = runClustal('analysis/israel-ref.txt', reference_id, 4)

    # Muscle Parser israel
    MuscleJ = runMuscle('analysis/muscle-I20200523-085708-0753-28920419-p1m.clw', reference_id, 4)

    diff = utils.jsonComp('output/'+ ClustalJ, 'output/'+ MuscleJ)
    #print("Israel " + str(diff))    # print("Israel " + str(data))
    utils.saveCompareFile("differences.txt", "Israel ", diff)

    # Clustal Parser GISAID
    ClustalJ = runClustal('analysis/GISAID-all.txt', reference_id, 7)

    # Muscle Parser GISAID
    MuscleJ = runMuscle('analysis/muscle-I20200523-090216-0023-41230765-p1m.clw', reference_id, 7)

    diff = utils.jsonComp('output/'+ ClustalJ, 'output/'+ MuscleJ)
    #print("GISAID " + str(diff))
    utils.saveCompareFile("differences.txt", "GISAID ", diff)

    #Clustal Parser ncbi
    ClustalJ = runClustal('analysis/all.txt', reference_id, 3)

    # Muscle Parser ncbi
    MuscleJ = runMuscle('analysis/muscle-I20200512-170208-0225-69454386-p1m.clw', reference_id, 8)

    #compare them
    data = utils.jsonComp('output/'+ ClustalJ, 'output/'+ MuscleJ)
    #print("ncbi " + str(data))
    utils.saveCompareFile("differences.txt", "ncbi ", diff)

    # Clustal Global Parser
    ClustalJ = runClustal('analysis/global.txt', reference_id, 14)

    # Muscle Global Parser
    MuscleJ = runMuscle('analysis/muscle-I20200523-090837-0910-95910164-p1m.clw', reference_id, 15)

    #compare them
    data = utils.jsonComp('output/'+ ClustalJ, 'output/'+ MuscleJ)
    #print("global " + str(data))
    utils.saveCompareFile("differences.txt","Global ", diff)


if __name__ == "__main__":
    main()
