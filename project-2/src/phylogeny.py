"""from https://github.com/mcmero/perfect_phylogeny"""
#download code and test prints
    #requrements error

"""from https://biopython.org/wiki/Phylo
from Bio import Phylo

def test(file)
    trees = Phylo.parse('phyloxml_examples.xml', 'phyloxml')
        for tree in trees:
            print(tree.name)"""

"""custom with https://pypi.org/project/anytree/ , https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.dendrogram.html"""
from anytree import Node, RenderTree
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np

### form https://www.geeksforgeeks.org/check-if-two-arrays-are-equal-or-not/ ###
#check if 2 arrays are identical [used to check index sorted correctly]
def areEqual(arr1, arr2, n, m): 
  
    # If lengths of array are not  
    # equal means array are not equal 
    if (n != m): 
        return False; 
  
    # Linearly compare elements 
    for i in range(0, n - 1): 
        if (arr1[i] != arr2[i]): 
            return False; 
  
    # If all elements were same. 
    return True; 

def SortDF(dataframe, ascending=False):
    """sort pandas.dataframe by row (descending by value) and return sorted one
    (value = number of 1s)"""
    ### array counting '1' for row ###
    #count = dataframe.astype(bool).sum(axis=0) ##questo lo fa per colonna
    count = dataframe.astype(bool).sum(axis=1) ##questo lo fa per riga
    #print(count)
    ### sort dataframe by number of 1s """https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.Series.sort_values.html"""###
    count.sort_values(ascending=True, inplace=True, ignore_index=False) #inplace or doesn't work
    #print('Sorted count:\n', count)
    #print(count.index)
    #print(dataframe.reindex(count.index).index)
    #print(areEqual(count.index, dataframe.reindex(count.index).index, len(count.index), len(dataframe.reindex(count.index).index)))
    #dataframe.sort_index(level = count, ascending=True) ###SOMETHING WRONG HERE
    if areEqual(count.index, dataframe.reindex(count.index).index, len(count.index), len(dataframe.reindex(count.index).index)):
        return dataframe.reindex(count.index)
    else:
        print("error")
        return dataframe

class Node:
    def __init__(self, name = "root", data = [], left = None, right = None):
        self.left = left
        self.right = right
        self.name = name
        self.data = data
    def __str__(self):  #print/tostring
        string = str(self.name) + ":" + str(' '.join([str(elem) for elem in self.data]))
        if (not self.isLeaf()):
            if (not self.left == None):
                string = '(' + self.left.__str__() + ')' + string 
            if (not self.right == None):
                string = string + '(' + self.right.__str__() + ')'
        return string
    def isLeaf(self):
        if self.left == None and self.right == None:
            return True
        else:
            return False
    def setLeft(self, node):
        self.left = node
    def setRight(self, node):
        self.right = node
#visit recursive
def visit(node):
    if node.isLeaf():
        return node
    if node.left is not None:
        visit(node.left)
    if node.right is not None:
        visit(node.right)

def test(dataframe, reference):
    """attempt anytree tre build"""
    #ref = Node(reference) #parent node tree

    """https://www.it-swarm.dev/it/python/come-posso-implementare-un-albero-python-ci-sono-delle-strutture-dati-integrate-python-come-java/967995793/
    For i in row
        if i == 1
            divide tree with iteritems() 
                if column has mutation [left]
                else rightt
        """
    #print(RenderTree(ref))
    """try normal tree
    devo trovare il modo di mettere sta parte in funzione ricorsiva
    e capire come stampare l'albero"""

    """root = Tree(name = reference, data = dataframe.head())
    for index, row in dataframe.iterrows():
        left = []
        right = []
        #print("index =", index)
        #print(row)
        for element in row:
            if element == 1:
                right.append(element)
            else:
                left.append(element)
        
        root.left = Tree(name = str('NOT' + index), data = left)
        root.right = Tree(name = index, data = right)"""
    """custom node test"""
    nr = Node("test1", data = [0,1,0,1])
    nl = Node("test2", data = [0,0,0,0])
    root = Node(data = [0,1,0,1,0,0,0,0], left= nl, right=nr )
    nr.setLeft(Node("test3", data = [1,1,1,1], right = (Node("test4", data = [1,0]))))

    print(root)


