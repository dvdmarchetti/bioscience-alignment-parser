import os   #IO
import json #Differences.json
import pandas as pd  #to read excels

def load_fasta(dir):
    with open(dir, 'r') as f:
        lines = f.readlines()
        id = lines[0].split(' ')[0][1:]
        seqlist = []
        for line in lines[1:]: #linee file
            if not (line.startswith(">") and line.startswith('\n')): #same gene
                seqlist.append(line.rstrip()) #remove \n for each line
        sequence = ''.join(seqlist[0:])

    return str(id + '\n' + sequence)

def main():
    with open(os.path.join('..', '..', 'project-1', 'input', 'reference.fasta'), 'r') as f:
        lines = f.readlines()
        reference_id = lines[0].split(' ')[0][1:] # NC_045512.2
        seqlist = []
        for line in lines[1:]: #linee file
            if not (line.startswith(">") and line.startswith('\n')): #same gene
                seqlist.append(line.rstrip()) #remove \n for each line
        reference = ''.join(seqlist[0:])

    with open(os.path.join(os.getcwd(), '..', '..', 'project-1', 'output', 'Clustal-NC_045512.2_2020-05-30_16-51.json')) as json_file:
        clustal_output = json.load(json_file)

    #sequence of sequences
    sequences = []
    #GISAID ones
    mydir = os.path.join('..', '..', 'project-1', 'input', 'GISAID')
    for file in os.listdir(mydir):
        if file.endswith(".fasta"):
            sequences.append(load_fasta(os.path.join(mydir, file)))
    #ncbi ones
    mydir = os.path.join('..', '..', 'project-1', 'input', 'ncbi')
    for file in os.listdir(mydir):
        if file.endswith(".fasta"):
            sequences.append(load_fasta(os.path.join(mydir, file)))
    
    print(len(sequences))
    #print(sequences[0][0:50])
    #print(sequences[7][0:50])
    #print(sequences[6][0:500])

    """se invece di controllare ogni singola sequenza ci segnassimo in una tabella da clustal_output
    quali id hanno variazioni con 0,1 come da descrizione e ...
    cercare algoritmo alberi filogenetici appena finisco gli altri esami"""

if __name__ == "__main__":
    main()
    print("not crashed")