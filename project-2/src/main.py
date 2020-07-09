import os
import json
import pandas as pd


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
    clustal_output = load_output('Clustal-NC_045512.2_2020-05-30_16-51.json')
    variations = clustal_output['unmatches'].items()

    [df_ref_genes, df_ref_cds] = [x for x in load_excel(os.path.join('..', 'Genes-CDS.xlsx'), ['Geni', 'CDS'])]
    df_ref_genes.set_index('Gene ID', inplace=True)
    df_ref_genes['Start'] -= 1
    # df_ref_genes['End'] -= 1
    df_ref_cds['from'] -= 1
    # df_ref_cds['to'] -= 1

    # For each variation:
    # 1) Find the CDS where it fall into
    # 2) Find the gene relative to the CDS
    # 3) Append informations to list used to build dataframe (more efficient than df.append for each row)
    variations_to_genes = []
    for key, value in variations:
        value['from'] -= 4
        value['to'] -= 4
        affected_cdses = df_ref_cds.loc[(value['from'] > df_ref_cds['from']) & (value['to'] < df_ref_cds['to'])]

        for index, cds in affected_cdses.iterrows():
            cds_sequence = reference[cds['from']:cds['to'] + 1]
            if pd.notna(cds['join_at']):
                join_at = int(cds['join_at'])
                cds_sequence = reference[cds['from']:join_at + 1] + reference[join_at:cds['to'] + 1]

            gene = df_ref_genes.loc[cds['GeneID']]
            gene_id = cds['GeneID']
            gene_start = gene['Start']
            gene_end = gene['End']
            cds_start = cds['from']
            cds_end = cds['to']
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

            original_aminoacid = ''
            encoded_aminoacid = ''
            original_codone = cds_sequence[alteration_start:alteration_end]
            altered_codone = cds_sequence[alteration_start:relative_start] + str(sequence) + cds_sequence[relative_end:alteration_end]
            # original_codone = reference[global_alteration_start: global_alteration_end]
            # altered_codone = reference[global_alteration_start:value['from']] + str(sequence) + reference[value['to']:global_alteration_end]
            for i in range(0, len(altered_codone), 3):
                group = altered_codone[i:i+3]
                # ref_group = reference[global_alteration_start+i:global_alteration_start+i+3]
                ref_group = cds_sequence[alteration_start+i:alteration_start+i+3]
                if '-' not in group:
                    encoded_aminoacid += translate_into_aminoacid(group)
                    original_aminoacid += translate_into_aminoacid(ref_group)

            variations_to_genes.append({
                'gene_id': gene_id,
                'gene_start': gene_start + 1,
                'gene_end': gene_end,
                'cds_start': cds_start + 1,
                'cds_end': cds_end,
                'original_codone': original_codone,
                'altered_codone': altered_codone,
                'relative_start': relative_start + 1,
                'relative_end': relative_end,
                'alteration': sequence,
                'original_aminoacid': original_aminoacid,
                'encoded_aminoacid': encoded_aminoacid
            })

    # Convert the list to dataframe
    columns = ['gene_id', 'gene_start', 'gene_end', 'cds_start', 'cds_end', 'original_codone', 'altered_codone', 'relative_start', 'relative_end', 'alteration', 'original_aminoacid', 'encoded_aminoacid']
    df_variations_to_genes = pd.DataFrame(variations_to_genes, columns=columns)
    df_variations_to_genes.to_csv(os.path.join('..', 'output', 'out.csv'))
    # print(df_variations_to_genes)


def load_output(filename):
    """
    Load a json file from output folder of project 1 given its filename
    """
    path = os.path.join(os.getcwd(), '..', '..', 'project-1', 'output', filename)
    with open(path) as json_file:
        return json.load(json_file)


def load_excel(path, sheets=None):
    """
    Load and excel file given its path. If sheets param is specified load only these specific sheets.
    """
    xls = pd.ExcelFile(os.path.join(os.getcwd(), path))

    if sheets is None:
        return pd.read_excel(path)

    for sheet in sheets:
        yield pd.read_excel(xls, sheet)


# Dict RNA Translation Table
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
        if rna in value:
            return key

    raise Exception('Invalid protein: {}'.format(rna))


if __name__ == "__main__":
    main()
