import os   #IO
import json #Differences.json
import pandas as pd  #to read excels


def load_output(filename):
    """
    Load a json file from output folder of project 1 given its filename
    """
    path = os.path.join(os.getcwd(), '..', '..', 'project-1', 'output', filename)
    with open(path) as json_file:
        return json.load(json_file) #json.load reads the string from the file, parses the JSON data, populates a Python dict with the data and returns it back to you.

def load_excel(path, sheets=None):
    """
    Load and excel file given its path. If sheets param is specified load only these specific sheets.
    """
    xls = pd.ExcelFile(os.path.join(os.getcwd(), path))

    if sheets is None:
        return pd.read_excel(path)

    for sheet in sheets:
        yield pd.read_excel(xls, sheet)

#dict RNA Translation as global variable
Amino_acids = {
    "START" : 'ATG',
    "STOP" : ["TAA", "TAG", "TGA"],
    'F' : ['TTT', 'TTC'],
    'L' : ['TTA', 'TTG','CTT', 'CTA', 'CTC', 'CTG'],
    'I' : ['ATT', 'ATC', 'ATA'],
    'M' : ['ATG'],
    'V' : ['GTT', 'GTA', 'GTC', 'GTG'],
    'S' : ['TCT', 'TCA', 'TCC', 'TCG'],
    'P' : ['CCT', 'CCA', 'CCC', 'CCG'],
    'T' : ['ACT', 'ACA', 'ACC', 'ACG'],
    'A' : ['GCT', 'GCA', 'GCC', 'GCG'],
    'Y' : ['TAT', 'TAC'],
    'H' : ['CAT', 'CAC'],
    'Q' : ['CAA', 'CAG'],
    'N' : ['AAT', 'AAC'],
    'K' : ['AAA', 'AAG'],
    'D' : ['GAT', 'GAC'],
    'E' : ['GAA', 'GAG'],
    'C' : ['TGT', 'TGC'],
    'W' : ['TGG'],
    'R' : ['CGT', 'CGA', 'CGC', 'CGG', 'AGA', 'AGG'],
    'S' : ['AGT', 'AGC'],
    'G' : ['GGT', 'GGA', 'GGC', 'GGG']
}


def main():
    reference_id = open(os.path.join('..', '..', 'project-1', 'input', 'reference.fasta'), 'r').readline().split(' ')[0][1:] # NC_045512.2
    #print(reference_id)
    reference = open(os.path.join('..', '..', 'project-1', 'input', 'reference_ncbi_gisaid_sequences.fasta'), 'r').readlines()[1:]
    seqlist=[]
    for line in reference: #linee file
        if not (line.startswith(">") and line.startswith('\n')): #same gene
            seqlist.append(line.rstrip()) #remove \n for each line
    #unify lines in single string
    reference = seqlist[0].join(seqlist[0:])
    #print(seqlist[0:20])
    del(seqlist) #usless
    #print(reference[266:13483])
    #print(reference[1135])
    #print(reference[0:30])

    # Read input files (Json + Excel)
    muscle_output = load_output('Muscle-NC_045512.2_2020-05-30_16-51.json')
    alignment_output = load_output('Muscle_Clustal-NC_045512.2_2020-05-30_16-51.json')
    clustal_output = load_output('Clustal-NC_045512.2_2020-05-30_16-51.json')

    variations = clustal_output['unmatches'].items()
    #variations = alignment_output['unmatches'].items()
    #variations = clustal_output['unmatches'].items()

    [df_ref_genes, df_ref_cds] = [x for x in load_excel(os.path.join('..', 'Genes-CDS.xlsx'), ['Geni', 'CDS'])]
    df_ref_genes.set_index('Gene ID', inplace=True)
    #print(df_ref_genes)
    #print()
    #print(df_ref_cds)
    #print()

    # For each variation:
    # 1) Find the CDS where it fall into
    # 2) Find the gene relative to the CDS
    # 3) Append informations to list used to build dataframe (more efficient than df.append for each row)
    variations_to_genes = []
    for key, value in variations:
        affected_cds = df_ref_cds.loc[(value['from'] > df_ref_cds['from']) & (value['to'] < df_ref_cds['to'])]

        #if affected_cds.empty:
            #continue
        ### Break & continue non mi piacciono
        if not affected_cds.empty:
            gene = df_ref_genes.loc[affected_cds['GeneID']]
            gene_id = affected_cds.iloc[0]['GeneID']
            gene_start = gene.iloc[0]['Start']
            gene_end = gene.iloc[0]['End']
            cds_start = affected_cds.iloc[0]['from']
            cds_end = affected_cds.iloc[0]['to']
          ###MOVED altered_codone at the end
            sequence = value['alt']
            relative_start = value['from'] - gene_start + 1
            relative_end = relative_start + len(sequence) - 1
            ###MOVED encoded_aminoacid after altered_codone       
            altered = ((relative_start - cds_start)%3) #closest multiple by 3 to relative_start - cds_start
            original = reference[relative_start:relative_end]
            reference = reference[:relative_start] + str(sequence[0:relative_end-relative_start]) + reference[relative_start+(relative_end-relative_start)+1:] 
            altered_codone = reference[relative_start - altered: relative_start - altered + 3*(relative_start-relative_end + 1)]#'' #TODO
            reference = reference[:relative_start] + str(original[0:relative_end-relative_start]) + reference[relative_start+(relative_end-relative_start)+1:]  #to not change original reference 
            encode = [] #if multiple amminoacids
            for key, value in Amino_acids.items():
                if altered_codone in value:
                    encode.append(key)
            encoded_aminoacid = ''.join(encode) #TODO
            
            #DEBUG# Amino_acids non legge 'S' WTF?
            if encoded_aminoacid == '':
                print(relative_end-relative_start, altered_codone)  
            if altered_codone == '':
                print(relative_end-relative_start, sequence)   
            #END DEBUG 
            
            variations_to_genes.append({
                'gene_id': gene_id,
                'gene_start': gene_start,
                'gene_end': gene_end,
                'cds_start': cds_start,
                'cds_end': cds_end,
                'altered_codone': altered_codone,
                'relative_start': relative_start,
                'relative_end': relative_end,
                'sequence': sequence,
                'encoded_aminoacid': encoded_aminoacid,
            })

    # Convert the list to dataframe
    columns = ['gene_id', 'gene_start', 'gene_end', 'cds_start', 'cds_end', 'altered_codone', 'relative_start', 'relative_end', 'sequence', 'encoded_aminoacid']
    df_variations_to_genes = pd.DataFrame(variations_to_genes, columns=columns)
    print(df_variations_to_genes)

    #print(Amino_acids)

if __name__ == "__main__":
    main()
