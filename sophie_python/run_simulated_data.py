#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 10:22:29 2022

@author: fatemehmohebbi
"""


import numpy as np
from Bio import Phylo
import networkx as nx
from scipy import sparse
import glob, os, sys
from utils import tree_reduction, get_consensus_net, plot_network
from likelihood import network_likelihood, network_likelihood_parallel
from felsenstein_pruning import felsenstein, felsenstein_parallel


def read_tree(tree_, tree_type):
    """read the tree file and return adjacency and weight matrices, 
    and the list of patients"""
    tree = Phylo.read(tree_, tree_type)
    graph = Phylo.to_networkx(tree)
    #adjacency_mat.toarray() to read the mat
    weight_mat = nx.to_numpy_matrix(graph, weight='weight')
    adjacency_mat = sparse.csr_matrix(np.sign(weight_mat))
    patients_ids= []
    for node in graph.nodes:
        if node.name is not None:
            id_ = str(node.name).split('|')[1]
            patients_ids.append(id_)
        else:
            patients_ids.append('None')
    return adjacency_mat, weight_mat, patients_ids


def infer_trans_network(adjacency_mat, weight_mat, patients):
    adjacency_reduced, weight_mat_reduced, patients = \
        tree_reduction(adjacency_mat, weight_mat, patients)
    num_of_patients = len(set(patients)) - 1
    patients_list = range(num_of_patients)
    if processes > 1:
        sampled_labels, sampled_likelihood = felsenstein_parallel(iterations, patients, weight_mat_reduced, mu, processes)
        matching_likelihood, networks = network_likelihood_parallel(distr_type, deg_distr, sampled_labels, \
                                                            weight_mat_reduced,patients_list, enforce_tree, processes)
    else:
        sampled_labels, sampled_likelihood = felsenstein(iterations, patients, weight_mat_reduced, mu)
        matching_likelihood, networks = network_likelihood(distr_type, deg_distr, sampled_labels, \
                                                        weight_mat_reduced,patients_list, enforce_tree)
        
    total_likelihood = sampled_likelihood + np.array(matching_likelihood)
    
    if cons_type == 'joint':
        consensus_likelihood = total_likelihood[:]
    if cons_type == 'network':
        consensus_likelihood = np.array(matching_likelihood)
    if cons_type == 'phylogenetic':
        consensus_likelihood = sampled_likelihood[:]
    the_consensus = get_consensus_net(networks, consensus_likelihood)
    plot_network(the_consensus, patients, output_directory)

    return the_consensus


if __name__ == '__main__':
    file_path = sys.argv[1]
    output_directory = sys.argv[2]
    
    tree_type = "newick" #"nexus" 
    mu = 0.05 #trans_rate
    iterations = int(sys.argv[3])
    cons_type = 'joint'
    distr_type = 'power law'
    deg_distr = 2.00
    enforce_tree = 0
    processes = int(sys.argv[4])
    
    if output_directory != "none":
        adjacency_mat, weight_mat, patients_ids = read_tree(file_path, tree_type)
        infer_trans_network(adjacency_mat, weight_mat, patients_ids)
    else:
        for folder in glob.glob(os.path.join(file_path, '*')):
            file = folder + "/error_free_files/phylogenetic_trees/tree_0.time.tre"
            print(file)
            output_directory = folder + "/error_free_files/phylogenetic_trees/"
            adjacency_mat, weight_mat, patients_ids = read_tree(file, tree_type)
            infer_trans_network(adjacency_mat, weight_mat, patients_ids)
        
        
