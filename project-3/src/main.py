import os
import json
import numpy as np
import pandas as pd
import phylogeny


def main():
    # 1. Read fasta sequences
    reference_id = load_fasta_id(os.path.join('..', '..', 'project-1', 'input', 'reference.fasta'))
    sequence_ids = read_sequence_ids(paths=[
        os.path.join('..', '..', 'project-1', 'input', 'GISAID'),
        os.path.join('..', '..', 'project-1', 'input', 'ncbi'),
    ])
    sequence_ids.insert(0, reference_id)    #insert reference no variations

    assert len(sequence_ids) == 14

    # 2. Read variations (project 1 output)
    clustal_output = load_output('Clustal-NC_045512.2_2020-07-21_15-14.json')
    variations = clustal_output['unmatches'].items()

    # 3. Build trait matrix
    indexes = []
    rows = []
    counter = 1
    for key, value in variations:
        row = np.zeros(len(sequence_ids))
        indexes.append('C{}'.format(counter))
        for sequence in value['sequences']:
            row[sequence_ids.index(sequence)] = 1
        rows.append(row)
        counter += 1

    # 4. Build tree
    trait_matrix = pd.DataFrame(rows, index=indexes, columns=sequence_ids, dtype=np.uint8).transpose()
    trait_matrix = phylogeny.reorder_columns(trait_matrix, axis=0, ascending=False)
    trait_matrix.to_csv(os.path.join('..', 'output', 'character_table.csv'))

    candidate_matrix = get_perfect_phylogeny_character_matrix(trait_matrix)
    if phylogeny.is_forbidden_matrix(candidate_matrix):
        raise Exception('Invalid perfect phylogeny matrix')
    #print(phylogeny.is_forbidden_matrix(candidate_matrix))
    candidate_matrix.to_csv(os.path.join('..', 'output', 'perfect_phylogeny_table.csv'))

    phylogeny.build_tree(candidate_matrix)


def get_perfect_phylogeny_character_matrix(df):
    """parameters = boolean dataframe
    for each column, add dataframe if only doesn't create prohibited matrix
    return resulted matrix (with forbidden column removed)
    """
    columns = df.columns
    candidate_matrix = df[columns[0:1]]
    for i in range(1, len(columns)):
        candidate_matrix = candidate_matrix.join(df[columns[i:i+1]])

        if phylogeny.is_forbidden_matrix(candidate_matrix):
            candidate_matrix = candidate_matrix.drop(labels=candidate_matrix.columns[-1], axis=1)

    return candidate_matrix


def read_sequence_ids(paths=[]):
    ids = []

    for path in paths:
        for f in os.listdir(path):
            if f.endswith('.fasta'):
                ids.append(load_fasta_id(os.path.join(path, f)))

    return ids


def load_fasta_id(path):
    with open(path, 'r') as f:
        line = f.readline()
        return line.split('|')[1] if '|' in line else line.split(' ')[0][1:]


def load_output(filename):
    """
    Load a json file from output folder of project 1 given its filename
    """
    path = os.path.join(os.getcwd(), '..', '..', 'project-1', 'output', filename)
    with open(path) as json_file:
        return json.load(json_file)


if __name__ == "__main__":
    main()
