import os
import json
import pandas as pd
import matplotlib.pyplot as plt


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
    clustal_output = load_output('Clustal-NC_045512.2_2020-07-21_15-14.json')
    variations = clustal_output['unmatches'].items()

    [df_ref_genes, df_ref_cds] = [x for x in load_excel(os.path.join('..', 'Genes-CDS.xlsx'), ['Geni', 'CDS'])]
    df_ref_genes.set_index('Gene ID', inplace=True)
    # df_ref_genes['Start'] -= 1
    # df_ref_genes['End'] -= 1
    # df_ref_cds['from'] -= 1
    # df_ref_cds['to'] -= 1

    # For each variation:
    # 1) Find the CDS where it fall into
    # 2) Find the gene relative to the CDS
    # 3) Append informations to list used to build dataframe (more efficient than df.append for each row)
    variations_to_genes = []
    gene_sequences = []
    var_not_var = []
    for key, value in variations:
        value['from'] -= 4
        value['to'] -= 4
        value['alt'] = value['alt'].replace('T', 'U')
        affected_cdses = df_ref_cds.loc[(value['from'] > df_ref_cds['from']) & (value['to'] < df_ref_cds['to'])]
        in_cds = False

        for index, cds in affected_cdses.iterrows():
            in_cds = True
            cds_sequence = reference[cds['from']:cds['to'] + 1]
            if pd.notna(cds['join_at']):
                join_at = int(cds['join_at'])
                cds_sequence = reference[cds['from']:join_at + 1] + reference[join_at:cds['to'] + 1]

            cds_sequence = cds_sequence.replace('T', 'U')
            gene = df_ref_genes.loc[cds['GeneID']]
            gene_id = cds['gene']
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
                original_aminoacid += translate_into_aminoacid(ref_group)
                if '-' not in group:
                    encoded_aminoacid += translate_into_aminoacid(group)

            variations_to_genes.append({
                'gene_id': gene_id,
                'gene_start': gene_start,
                'gene_end': gene_end,
                'cds_start': cds_start,
                'cds_end': cds_end,
                'original_codone': original_codone,
                'altered_codone': altered_codone,
                'relative_start': relative_start,
                'relative_end': relative_end,
                'alteration': sequence,
                'original_aminoacid': original_aminoacid,
                'encoded_aminoacid': encoded_aminoacid,
            })

            for sequence in value['sequences']:
                gene_sequences.append({
                    'gene_id': gene_id,
                    'sequence': sequence
                })

        var_not_var.append({
            'variation': key,
            'in_cds': in_cds
        })

    # Convert the list to dataframe
    columns = ['gene_id', 'gene_start', 'gene_end', 'cds_start', 'cds_end', 'relative_start', 'relative_end', 'alteration', 'original_codone', 'original_aminoacid', 'altered_codone', 'encoded_aminoacid']
    df_variations_to_genes = pd.DataFrame(variations_to_genes, columns=columns).set_index('gene_id')
    df_variations_to_genes.to_csv(os.path.join('..', 'output', 'alteration-table.csv'))

    colors = ['mediumvioletred', 'darkorange', 'dodgerblue', 'green', 'coral', 'mediumpurple', 'gold', 'm', 'olivedrab', 'salmon','yellowgreen', 'c', 'lightgray']
    chars_per_var = df_variations_to_genes.groupby(['gene_id']).count().sort_values(by='alteration')
    ax = chars_per_var.plot.pie(y='alteration', autopct='%.2f%%', colors=colors, legend=False)
    # ax.set_title('Variations per gene')
    ax.set_ylabel('')
    plt.tight_layout()
    # plt.show()
    plt.savefig(os.path.join('..', 'relazione', 'images', 'plot-variations-per-gene'))

    # plot
    columns = ['variation', 'in_cds']
    df_in_cds = pd.DataFrame(var_not_var, columns=columns)

    chars_per_var = df_in_cds.groupby(['in_cds']).count()
    ax = chars_per_var.plot.pie(y='variation', autopct='%.2f%%', labels=['Outside all CDS', 'Inside at least one CDS'], colors=colors, legend=False)
    ax.set_ylabel('')
    plt.tight_layout()
    # plt.show()
    plt.savefig(os.path.join('..', 'relazione', 'images', 'plot-alteration-in-cds'))

    # plot
    columns = ['gene_id', 'sequence']
    df_gene_sequences = pd.DataFrame(gene_sequences, columns=columns)

    chars_per_var = df_gene_sequences.groupby(['sequence']).count()
    ax = chars_per_var.plot.pie(y='gene_id', autopct='%.2f%%', colors=colors, legend=False, explode=[0.1] * 13)
    ax.set_ylabel('')
    plt.tight_layout()
    # plt.show()
    plt.savefig(os.path.join('..', 'relazione', 'images', 'plot-alteration-in-gene-per-sequence'))

    plt.close()


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
    'F': ['UUU', 'UUC'],
    'L': ['UUA', 'UUG', 'CUU', 'CUA', 'CUC', 'CUG'],
    'I': ['AUU', 'AUC', 'AUA'],
    'M': ['AUG'],
    'V': ['GUU', 'GUA', 'GUC', 'GUG'],
    'S': ['UCU', 'UCA', 'UCC', 'UCG', 'AGU', 'AGC'],
    'P': ['CCU', 'CCA', 'CCC', 'CCG'],
    'T': ['ACU', 'ACA', 'ACC', 'ACG'],
    'A': ['GCU', 'GCA', 'GCC', 'GCG'],
    'Y': ['UAU', 'UAC'],
    'H': ['CAU', 'CAC'],
    'Q': ['CAA', 'CAG'],
    'N': ['AAU', 'AAC'],
    'K': ['AAA', 'AAG'],
    'D': ['GAU', 'GAC'],
    'E': ['GAA', 'GAG'],
    'C': ['UGU', 'UGC'],
    'W': ['UGG'],
    'R': ['CGU', 'CGA', 'CGC', 'CGG', 'AGA', 'AGG'],
    'G': ['GGU', 'GGA', 'GGC', 'GGG'],
    'START': ['AUG'],
    'STOP': ['UAA', 'UAG', 'TGA']
}


def translate_into_aminoacid(rna):
    for key, value in aminoacids_lookup_table.items():
        if rna in value:
            return key

    raise Exception('Invalid protein: {}'.format(rna))


if __name__ == "__main__":
    main()
