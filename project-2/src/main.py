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
aminoacids_lookup_table = {
    "START" : 'ATG',
    "STOP" : ["TAA", "TAG", "TGA"],
    'F' : ['TTT', 'TTC'],
    'L' : ['TTA', 'TTG', 'CTT', 'CTA', 'CTC', 'CTG'],
    'I' : ['ATT', 'ATC', 'ATA'],
    'M' : ['ATG'],
    'V' : ['GTT', 'GTA', 'GTC', 'GTG'],
    'S' : ['TCT', 'TCA', 'TCC', 'TCG', 'AGT', 'AGC'],
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
    'G' : ['GGT', 'GGA', 'GGC', 'GGG']
}


def translate_into_aminoacid(rna):
    for key, value in aminoacids_lookup_table.items():
        print (rna, key, value)
        if rna in value:
            return key

    raise Exception('Invalid protein: {}'.format(rna))

def main():
    reference_id = None
    with open(os.path.join('..', '..', 'project-1', 'input', 'reference.fasta'), 'r') as f:
        reference_id = f.readline().split(' ')[0][1:] # NC_045512.2

    reference = None
    with open(os.path.join('..', '..', 'project-1', 'input', 'reference_ncbi_gisaid_sequences.fasta'), 'r') as f:
        seqlist = []
        for line in f.readlines()[1:]: #linee file
            if not (line.startswith(">") and line.startswith('\n')): #same gene
                seqlist.append(line.rstrip()) #remove \n for each line
        reference = seqlist[0].join(seqlist[0:])

    # Read input files (Json + Excel)
    muscle_output = load_output('Muscle-NC_045512.2_2020-05-30_16-51.json')
    alignment_output = load_output('Muscle_Clustal-NC_045512.2_2020-05-30_16-51.json')
    clustal_output = load_output('Clustal-NC_045512.2_2020-05-30_16-51.json')

    #variations = muscle_output['unmatches'].items()
    #variations = alignment_output['unmatches'].items()
    variations = clustal_output['unmatches'].items()

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

        if not affected_cds.empty:
            gene = df_ref_genes.loc[affected_cds['GeneID']]
            gene_id = affected_cds.iloc[0]['GeneID']
            gene_start = gene.iloc[0]['Start']
            gene_end = gene.iloc[0]['End']
            cds_start = affected_cds.iloc[0]['from']
            cds_end = affected_cds.iloc[0]['to']
            sequence = value['alt']
            relative_start = value['from'] - gene_start
            relative_end = relative_start + len(sequence)

            # CAG AAG CTA
            #      ^^ ^ alteration
            #     ^ altered group start
            alteration_start = relative_start - ((relative_start + 1) % 3)
            global_alteration_start = alteration_start + gene_start

            alteration_end = relative_end + (3 - ((relative_end + 1) % 3))
            if (relative_end + 1) % 3 == 0:
                alteration_end -= 3
            global_alteration_end = alteration_end + gene_start

            # print(global_alteration_start, value['from'], value['to'], global_alteration_end)
            # print(reference[global_alteration_start:global_alteration_end])
            # print(reference[global_alteration_start:value['from']] + str(sequence) + reference[value['to']:global_alteration_end])
            # print()

            encoded_aminoacid = ''
            altered_codone = reference[global_alteration_start:value['from']] + str(sequence) + reference[value['to']:global_alteration_end]
            for i in range(0, len(altered_codone), 3):
                group = altered_codone[i:i+3]
                if '-' not in group:
                    encoded_aminoacid += translate_into_aminoacid(group)
            # print (len(altered_codone), len(altered_codone) / 3, encoded_aminoacid)

            # #DEBUG# Amino_acids non legge 'S' WTF?
            # if encoded_aminoacid == '':
            #     print(relative_end-relative_start, altered_codone)
            # if altered_codone == '':
            #     print(relative_end-relative_start, sequence)
            # #END DEBUG

            variations_to_genes.append({
                'gene_id': gene_id,
                'gene_start': gene_start,
                'gene_end': gene_end,
                'cds_start': cds_start,
                'cds_end': cds_end,
                'altered_codone': altered_codone,
                'relative_start': relative_start + 1,
                'relative_end': relative_end,
                'sequence': sequence,
                'encoded_aminoacid': encoded_aminoacid,
            })

    # Convert the list to dataframe
    columns = ['gene_id', 'gene_start', 'gene_end', 'cds_start', 'cds_end', 'altered_codone', 'relative_start', 'relative_end', 'sequence', 'encoded_aminoacid']
    df_variations_to_genes = pd.DataFrame(variations_to_genes, columns=columns)
    print(df_variations_to_genes)

    #print(Amino_acids)
    print("not crashed")

if __name__ == "__main__":
    main()
