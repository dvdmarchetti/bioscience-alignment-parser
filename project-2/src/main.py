import os   #IO
import json #Differences.json
import pandas as pd  #to read excels
import phylogeny

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
        #print (rna, key, value)
        if rna in value:
            return key

    raise Exception('Invalid protein: {}'.format(rna))

def load_fasta_id(dir):
    with open(dir, 'r') as f:
        line = f.readline()
        id = line.split(' ')[0][1:]

    return str(id)


def main():
    reference_id = None
    reference = None
    with open(os.path.join('..', '..', 'project-1', 'input', 'reference.fasta'), 'r') as f:
        lines = f.readlines()
        reference_id = lines[0].split(' ')[0][1:] # NC_045512.2
        seqlist = []
        for line in lines[1:]: #linee file
            if not (line.startswith(">") and line.startswith('\n')): #same gene
                seqlist.append(line.rstrip()) #remove \n for each line
        reference = ''.join(seqlist[0:])

    # Read input files (Json + Excel)
    # muscle_output = load_output('Muscle-NC_045512.2_2020-05-30_16-51.json')
    # alignment_output = load_output('Muscle_Clustal-NC_045512.2_2020-05-30_16-51.json')
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
        value['from'] -= 4
        value['to'] -= 4
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

            original_aminoacid = ''
            encoded_aminoacid = ''
            original_codone = reference[global_alteration_start: global_alteration_end]
            altered_codone = reference[global_alteration_start:value['from']] + str(sequence) + reference[value['to']:global_alteration_end]
            for i in range(0, len(altered_codone), 3):
                group = altered_codone[i:i+3]
                ref_group = reference[global_alteration_start+i:global_alteration_start+i+3]
                if '-' not in group:
                    encoded_aminoacid += translate_into_aminoacid(group)
                    original_aminoacid += translate_into_aminoacid(ref_group)

            variations_to_genes.append({
                'gene_id': gene_id,
                'gene_start': gene_start,
                'gene_end': gene_end,
                'cds_start': cds_start,
                'cds_end': cds_end,
                'original_codone': original_codone,
                'altered_codone': altered_codone,
                'relative_start': relative_start + 1,
                'relative_end': relative_end,
                'sequence': sequence,
                'original_aminoacid': original_aminoacid,
                'encoded_aminoacid': encoded_aminoacid
            })

    # Convert the list to dataframe
    columns = ['gene_id', 'gene_start', 'gene_end', 'cds_start', 'cds_end', 'original_codone', 'altered_codone', 'relative_start', 'relative_end', 'sequence', 'original_aminoacid', 'encoded_aminoacid']
    df_variations_to_genes = pd.DataFrame(variations_to_genes, columns=columns)
    df_variations_to_genes.to_csv(os.path.join('..', 'output', 'out.csv'))
    print(df_variations_to_genes)

### END PART 2 ###############################################################
###PART 3 ####################################################################
    #print(variations)
    sequences = []
    #GISAID ones
    mydir = os.path.join('..', '..', 'project-1', 'input', 'GISAID')
    for file in os.listdir(mydir):
        if file.endswith(".fasta"): #.split('|')[1] perchè per qualche motivo nel file di output è tagliato in 
            sequences.append(load_fasta_id(os.path.join(mydir, file)).split('|')[1]) #questo modo e deve combaciare
    #ncbi ones
    mydir = os.path.join('..', '..', 'project-1', 'input', 'ncbi')
    for file in os.listdir(mydir):
        if file.endswith(".fasta"):
            sequences.append(load_fasta_id(os.path.join(mydir, file)))

    #print(len(sequences)) #must be 13
    #print(sequences)

    ###to put in already existing cicle in part 2
    indexes = []
    data3 = []
    for key, value in variations:
        #print(key, value)
        row = [0] * (len(sequences)) 
        indexes.append(key)
        for v in value['sequences']:
            row[sequences.index(v)]  = 1
        data3.append(row)

    #print(data3)
    """NEED NEW NAME"""
    df3 = pd.DataFrame(data3, index = indexes, columns=sequences)
    df3 = phylogeny.SortDF(df3) #sorting function
    #df3.to_csv(os.path.join('..', 'output', 'table.csv'))
    print(df3)
    phylogeny.test(df3, reference)

### END PART 3 ###############################################################

if __name__ == "__main__":
    main()  

