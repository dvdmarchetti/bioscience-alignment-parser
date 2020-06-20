import os   #IO
import json #Differences.json
import pandas as pd  #to read excels

def in_CDS(_from, _to, data):
    mismatches = []
    for key, value in data['unmatches'].items():    #iterate and divide key, value of unmatche
        if value['from'] > _from and value['to'] < _to:
            mismatches.append(key)
    
    return mismatches

def main():
    reference_id = open("../../project-1/input/reference.fasta", "r").readline().split(' ')[0][1:] # NC_045512.2
    #print(reference_id)
    with open("../../project-1/output/Muscle_Clustal-NC_045512.2_2020-05-30_16-51.json")as json_file:
        data = json.load(json_file) #json.load reads the string from the file, parses the JSON data, populates a Python dict with the data and returns it back to you.
        #print(data)
    codons = [] #list of codons with unmatches

    dataframe_location = "../Genes-CDS.xlsx"; #location xlsx 
    df = pd.read_excel(dataframe_location)    #read file
    #print(type(df))
    #print(df)

    #append key unmatches inside codons
    for index, row in df.iterrows():
        codons.append(in_CDS(row['Start'], row['End'], data))

    print(codons)   #mismatches for each pair in df

    #dict RNA Translation
    Amino_acids = {
        "START" : 'AUG',
        "STOP" : ["UAA", "UAG", "UGA"],
        'F' : ['UUU', 'UUC'],
        'L' : ['UUA', 'UUG','CUU', 'CUA', 'CUC', 'CUG'],
        'I' : ['AUU', 'AUC', 'AUA'],
        'M' : ['AUG'],
        'V' : ['GUU', 'GUA', 'GUC', 'GUG'],
        'S' : ['UCU', 'UCA', 'UCC', 'UCG'],
        'P' : ['CCU', 'CCA', 'CCC', 'CCG'],
        'T' : ['ACU', 'ACA', 'ACC', 'ACG'],
        'A' : ['GCU', 'GCA', 'GCC', 'GCG'],
        'Y' : ['UAU', 'UAC'],
        'H' : ['CAU', 'CAC'],
        'Q' : ['CAA', 'CAG'],
        'N' : ['AAU', 'AAC'],
        'K' : ['AAA', 'AAG'],
        'D' : ['GAU', 'GAC'],
        'E' : ['GAA', 'GAG'],
        'C' : ['UGU', 'UGC'],
        'W' : ['UGG'],
        'R' : ['CGU', 'CGA', 'CGC', 'CGG', 'AGA', 'AGG'],
        'S' : ['AGU', 'AGC'],
        'G' : ['GGU', 'GGA', 'GGC', 'GGG']
    }
    #print(Amino_acids)

if __name__ == "__main__":
    main()