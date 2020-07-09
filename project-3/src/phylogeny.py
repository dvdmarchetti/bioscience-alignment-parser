import io

from anytree import Node, RenderTree, AsciiStyle
from Bio import Phylo
import numpy as np
import pandas as pd


# def build_tree_recursive(root, df):
#     #https://pandas.pydata.org/pandas-docs/stable/getting_started/10min.html
#     #print(type(root.data), root.data)
#     if df.empty or len(root.data) <= 1:
#         return
#     else:
#         left_data = []
#         right_data = []
#         #split left and right if have mutation
#         id = df.index[0] #name mutation
#         i = 0
#         #print("id = " + id)
#         for element in df.iloc[0]: #guarda df.iloc and df.loc:  By integer slices, acting similar to numpy/python;  Selecting on a multi-axis by label.
#             if element == 1:
#                 right_data.append(root.data[i]) #column name with mutation = 1
#             else:
#                 left_data.append(root.data[i])  #column name without mutation = 0
#             i = i+1
#         right = Node(id, parent = root, data = right_data)
#         left = Node("!" + id, parent = root, data = left_data)
#         #root.setLeft(left)
#         #root.setRight(right)
#         build_tree_recursive(right, df[1:])
#         build_tree_recursive(left, df[1:])

#         #for node in root.iter_path_reverse():
#             #print(node)


def build_tree(df, reference_id):
    """
    #custom node test#
    nr = Node("test1", data = [0,1,0,1])
    nl = Node("test2", data = [0,0,0,0])
    root = Node(data = [0,1,0,1,0,0,0,0], left= nl, right=nr )
    nr.setLeft(Node("test3", data = [1,1,1,1], right = (Node("test4", data = [1,0]))))
    #print(root.data)
    """

    # Algorithm on slide
    root = Node('root', edges={})
    # For each sequence
    for i, row in df.iterrows():
        current_node = root

        # For each alteration
        for j in range(len(row)):
            # If alteration is present in the current sequence
            if row.iloc[j]:
                # If current_node is already linked to the j-th variation
                if j in current_node.edges:
                    # Update the current node
                    current_node = current_node.edges[j]
                else:
                    # Create a new node u and link it with the last node of this sequence
                    u = Node('U-{}'.format(row.index[j]), edges={})
                    u.parent = current_node
                    current_node.edges[j] = u
                    current_node = u

        # We looped all the variations, insert the sequence node at the end of the chain
        current_node.name = i
        current_node.sequence = i

    print(RenderTree(root, style=AsciiStyle()).by_attr())

    with_intermediate_alteration_nodes = False
    if with_intermediate_alteration_nodes:
        newick_tree = to_newick_tree(root, intermediates=True)
    else:
        newick_tree = '({},{})'.format(root.sequence, to_newick_tree(root, intermediates=False))

    newick_tree = Phylo.read(io.StringIO(newick_tree), 'newick')
    newick_tree.ladderize()
    Phylo.draw_ascii(newick_tree)
    # Phylo.draw(newick_tree)



def to_newick_tree(node, intermediates=False):
    if node.is_leaf:
        return node.sequence

    if intermediates:
        return '({}, {})'.format(node.name, ','.join([to_newick_tree(child, intermediates) for child in node.children]))

    return '({})'.format(','.join([to_newick_tree(child, intermediates) for child in node.children]))


def is_forbidden_matrix(df):
    # Build auxiliary matrix
    M = reorder_columns(df).to_numpy()
    rows, columns = M.shape

    L = np.zeros(M.shape)
    for i in range(rows):
        k = -1
        for j in range(columns):
            if M[i,j] == 1:
                L[i,j] = k
                k = j + 1

    # Check if we have a forbidden matrix
    for j in range(columns):
        for i in range(rows):
            if L[i,j] != 0:
                for l in range(rows):
                    if L[l,j] != 0 and L[i,j] != L[l,j]:
                        return True

    return False


def reorder_columns(df, axis=0):
    """Sort a dataframe in descending number of ones by the specified axis (default: by columns).

    The sorting is not done in place, you must reassign the dataframe
    with the return value of this function call.
    """
    sorted_axis = df.sum(axis=axis).sort_values(ascending=False)
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
