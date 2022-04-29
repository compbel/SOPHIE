#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 14:45:53 2022

@author: fatemehmohebbi
"""

from Bio import Phylo
import networkx as nx
import numpy as np
import pandas as pd
from scipy import sparse
import scipy.special as sc
import matplotlib.pyplot as plt
import pydot
from networkx.drawing.nx_pydot import graphviz_layout
import math, statistics
from networkx.algorithms.tree import Edmonds


def find_parents(tree):
    parents = []
    for leaf in tree.find_clades(terminal=True, order="level"):
        parents.append(tree.get_path(leaf)[-2])
        


def to_adjacency_matrix(tree):
    """Create an adjacency matrix (NumPy array) from clades/branches in tree"""
    allclades = list(tree.find_clades(order="level"))
    lookup = {}
    for i, elem in enumerate(allclades):
        lookup[elem] = i
    adjmat = np.zeros((len(allclades), len(allclades)))
    for parent in tree.find_clades(terminal=False, order="level"):
        for child in parent.clades:
            adjmat[lookup[parent], lookup[child]] = 1
    if not tree.rooted:
        # Branches can go from "child" to "parent" in unrooted trees
        adjmat = adjmat + adjmat.transpose()
    return (allclades, np.matrix(adjmat))


def log_choose_k(n, k, flag):
    if flag == 'int':
        l = sum(np.log(np.array(range(1,n)))) - sum(np.log(np.array(range(1,n - k)))) - \
            sum(np.log(np.array(range(1,k))))
    if flag == 'real':
        l = sc.gammaln(n + 1) - sc.gammaln(k + 1) - sc.gammaln(n - k + 1)
    return l


def get_ids(patients):
    """the ids may be different for different trees, fix it"""
    ids = []
    for name in patients:
        try:
            ids.append(name.split('|')[1])
        except:
            ids.append(name)
    return ids


def read_tree(tree_, meta_file, tree_type):
    """read the tree file and return adjacency and weight matrices, 
    and the list of patients"""
    print('Reading data ...')
    metadata = pd.read_csv(meta_file)
    metadata = metadata.set_index('name')
    tree = Phylo.read(tree_, tree_type)
    graph = Phylo.to_networkx(tree)
    #adjacency_mat.toarray() to read the mat
    # weight_mat = nx.adjacency_matrix(graph, weight='weight')
    # adjacency_mat = sparse.csr_matrix(np.sign(weight_mat.toarray()))
    weight_mat = nx.to_numpy_matrix(graph, weight='weight')
    adjacency_mat = sparse.csr_matrix(np.sign(weight_mat))
    patients_ids= []
    # dates = []
    for node in graph.nodes:
        if node.name is not None:
            # dates.append(metadata['date'].loc[node.name])
            patients_ids.append(str(metadata['host'].loc[node.name]))
        else:
            # dates.append(np.nan)
            patients_ids.append('None')
    return adjacency_mat, weight_mat, patients_ids


def tree_reduction(adjacency_mat, weight_mat, patients):
    """the children of each node (leaves as children) are removed if they 
    have the same label"""
    # ids = get_ids(patients)
    print('Tree reduction ...')
    ids = patients[:]
    G = nx.from_numpy_matrix(np.triu(weight_mat), create_using=nx.DiGraph)
    for i in list(G.in_degree):
        if i[1] == 0:
            root = i[0]
    dfs_order = list(nx.dfs_preorder_nodes(G, source=root))
    dfs_order.reverse()
    #are labels matching nodes properly?
    nodes_removed = []
    for i in dfs_order:
        if G.out_degree(i) == 0:
            continue
        else:
            child = list(G.successors(i))
            if (ids[child[0]] == ids[child[1]]) & (ids[child[0]] != 'None'):
                ids[i] = ids[child[0]]
                nodes_removed.extend((child[0], child[1]))
                ids[child[0]] = '-1'
                ids[child[1]] = '-1'
    G.remove_nodes_from(nodes_removed)
    weight_mat = nx.to_numpy_matrix(G, weight='weight')
    adjacency_mat = sparse.csr_matrix(np.sign(weight_mat))
    return adjacency_mat, weight_mat, [x for x in ids if x != '-1']


def get_consensus_net(networks, networks_likelihood):
    """calculates consensus network for all networks given as input"""
    inf_nets = np.where(networks_likelihood == -math.inf)
    if len(inf_nets[0]) != 0:
        networks_likelihood = np.delete(networks_likelihood, inf_nets[0])
        networks = np.delete(networks, inf_nets[0], 0)
        
    mean_value = statistics.mean(networks_likelihood)
    avg_added_likel = networks_likelihood[:] - mean_value
    weight_mats = [i * math.exp(j) for i, j in zip(networks, avg_added_likel)]
    weight_sum = weight_mats[0]
    for net_i in weight_mats[1:]:
        weight_sum = np.add(weight_sum, net_i)
    G = nx.from_numpy_matrix(weight_sum, create_using=nx.DiGraph)
    edmonds = Edmonds(G)
    T = edmonds.find_optimum(attr='weight', kind='max')
    tree_adj_mat = nx.to_numpy_matrix(T) #nx.adjacency_matrix(T).toarray()
    return np.where(tree_adj_mat != 0, 1, 0)


def plot_network(the_consensus, patients, output_dir):
    """plot and save the transmission network"""
    G = nx.from_numpy_matrix(the_consensus, create_using=nx.DiGraph)
    pos = graphviz_layout(G, prog = 'dot') #"twopi")
    ids = sorted(list(set(patients)))[:-1]
    labels = {}
    count = 0
    for i in ids:
        labels[count] = str(i)
        count = count + 1
    nx.draw(G, pos)
    nx.draw_networkx_labels(G, pos, labels, font_size=7.5)
    plt.savefig(output_dir + 'transmission_network.svg')
    print('Transmission network is saved as transmission_network.svg')
    plt.show()
    