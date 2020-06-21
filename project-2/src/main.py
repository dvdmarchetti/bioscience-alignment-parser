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


def main():
    reference_id = open(os.path.join('..', '..', 'project-1', 'input', 'reference.fasta'), 'r').readline().split(' ')[0][1:] # NC_045512.2
    #print(reference_id)

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
            altered_codone = '' #TODO
            sequence = value['alt']
            relative_start = value['from'] - gene_start + 1
            relative_end = relative_start + len(sequence) - 1
            encoded_aminoacid = ''  #TODO

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
