import io
import os

from anytree import Node, RenderTree, AsciiStyle
from Bio import Phylo
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd


def build_tree(df):
    #create dictionary datas for rendering tree
    sequences_data = {
        'NC_045512.2': {
            'date': '17/01/2020',
            'location': 'China'
        },
        'MT320891.2': {
            'date': '10/04/2020',
            'location': 'Iran'
        },
        'MT281530.2': {
            'date': '04/04/2020',
            'location': 'Iran'
        },
        'EPI_ISL_442523': {
            'date': '09/03/2020',
            'location': 'Iran'
        },
        'EPI_ISL_437512': {
            'date': '26/03/2020',
            'location': 'Iran'
        },
        'MT276598.1': {
            'date': '02/04/2020',
            'location': 'Israel'
        },
        'MT276597.1': {
            'date': '02/04/2020',
            'location': 'Israel'
        },
        'EPI_ISL_447469': {
            'date': '14/04/2020',
            'patient_age': 47,
            'location': 'Israel',
            'source': 'Naso-pharyngeal'
        },
        'MT262993.1': {
            'date': '25/03/2020',
            'location': 'Pakistan'
        },
        'MT240479.1': {
            'date': '25/03/2020',
            'location': 'Pakistan'
        },
        'EPI_ISL_417444': {
            'date': '25/03/2020',
            'location': 'Pakistan'
        },
        'MT327745.1': {
            'date': '13/04/2020',
            'location': 'Turkey'
        },
        'EPI_ISL_437334': {
            'date': '24/04/2020',
            'location': 'Turkey'
        },
        'EPI_ISL_437317': {
            'date': '27/04/2020',
            'location': 'Turkey'
        }
    }

    # Algorithm on slide
    root = Node('root', edges={})
    # For each sequence
    for i, row in df.iterrows():
        current_node = root

        # For each alteration
        for j in range(len(row)):
            # If alteration is present in the current sequence
            if row.iloc[j]: #true in cell
                key = j
                # If current_node is already linked to the j-th variation
                # edges = dictionary with links from (current_node, u)
                if key in current_node.edges:
                    # Update the current node
                    current_node = current_node.edges[key]
                else:
                    # Create a new node u and link it with the last node of this sequence
                    u = Node('U-{}'.format(row.index[j]), edges={}, parent=current_node)
                    current_node.edges[key] = u
                    current_node = u
        # We looped all the variations, insert the sequence node at the end of the chain
        Node(i, parent=current_node)

    # Render tree with characters (for debug only)
    # print(RenderTree(root, style=AsciiStyle()).by_attr())

    # Bulid Newick tree for visualization
    newick_tree = to_newick_tree(root, intermediates=False)
    newick_tree = Phylo.read(io.StringIO(newick_tree), 'newick')

    # Display tree in ascii mode
    # newick_tree.ladderize() # to put leaf on the right
    # Phylo.draw_ascii(newick_tree)

    # Figure aesthetics
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(1, 1, 1)

    # Assign sequences data to the leaf
    root = newick_tree.clade
    merge_sequences_data(root, sequences_data) # Add sequence data to sequence leaves
    # Render with plt
    Phylo.draw(newick_tree, do_show=False, axes=ax)

    ax.set_xlabel('Number of alterations')
    ax.set_ylabel('Sequences')
    plt.tight_layout()
    plt.savefig(os.path.join('..', 'output', 'phylogenetic-tree.png'))
    plt.show()


def to_newick_tree(node, intermediates=False):
    """recursively create tree in order to format it for Phylo (newick format)
    choose if intermediate node have to be seen
    """
    if node.is_leaf:
        return node.name

    if intermediates:
        return '({}, {})'.format(node.name, ','.join([to_newick_tree(child, intermediates) for child in node.children]))

    return '({})'.format(','.join([to_newick_tree(child, intermediates) for child in node.children]))


def merge_sequences_data(node, sequences_data):
    if node.is_terminal():
        data = sequences_data[node.name]
        node.name = '{}\n {} in {}'.format(node.name, data['date'], data['location'])
        return

    [ merge_sequences_data(child, sequences_data) for child in node.clades ]


def is_forbidden_matrix(df):
    # Build auxiliary matrix
    M = reorder_columns(df).to_numpy() # sort the columns by decreasing number of 1’s
    rows, columns = M.shape

    L = np.zeros(M.shape)   # build the auxiliary matrix L defined as follows:
    for i in range(rows):
        k = -1
        for j in range(columns):    # if Mij = 0, then Lij = 0
            if M[i,j] == 1:         # if Mij = 1 then Lij = k,
                L[i,j] = k          # being k the rightmost column to the left of j such that Mik = 1,
                k = j + 1           # otherwise Lij = −1

    # Check if we have a forbidden matrix
    for j in range(columns):
        for i in range(rows):
            if L[i,j] != 0:
                for l in range(rows):
                    if L[l,j] != 0 and L[i,j] != L[l,j]:
                        return True # laminar --> support perfect phylogeny

    return False # not laminar --> no perfect phylogeny


def reorder_columns(df, axis=0, ascending=False):
    """Sort a dataframe in descending number of ones by the specified axis (default: by columns).

    The sorting is not done in place, you must reassign the dataframe
    with the return value of this function call.
    """
    sorted_axis = df.sum(axis=axis).sort_values(ascending=ascending)
    if axis == 0:
        return df[sorted_axis.index]

    return df.reindex(sorted_axis.index)


def test_forbidden_matrix():
    test_case_1 = pd.DataFrame(
        data=[
            [1,1,0,0,0],
            [0,0,1,0,1],
            [1,1,0,0,1],
            [0,0,1,1,0],
            [0,1,0,0,1],
        ],
        index=['S1', 'S2', 'S3', 'S4', 'S5'],
        columns=['C1', 'C2', 'C3', 'C4', 'C5']
    )
    assert is_forbidden_matrix(test_case_1) == True

    test_case_2 = pd.DataFrame(
        data=[
            [0,1],
            [1,0],
            [1,1],
        ],
        index=['S1', 'S2', 'S3'],
        columns=['C1', 'C2']
    )
    assert is_forbidden_matrix(test_case_2) == True

    test_case_3 = pd.DataFrame(
        data=[
            [1,1,0,0,0],
            [0,0,1,0,0],
            [1,1,0,0,1],
            [0,0,1,1,0],
            [0,1,0,0,0],
        ],
        index=['S1', 'S2', 'S3', 'S4', 'S5'],
        columns=['C1', 'C2', 'C3', 'C4', 'C5']
    )
    assert is_forbidden_matrix(test_case_3) == False
