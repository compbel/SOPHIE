
"""
Created on Wed Mar  2 16:48:49 2022

@author: fatemehmohebbi
"""


import numpy as np
import sys
from utils import read_tree, tree_reduction, get_consensus_net, plot_network
from likelihood import network_likelihood, network_likelihood_parallel
from felsenstein_pruning import felsenstein, felsenstein_parallel


def infer_trans_network(tree_file, meta_file, tree_type):
    adjacency_mat, weight_mat, patients = read_tree(tree_file, meta_file, tree_type)
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
    """inputs are:
        tree_file: the input phylogenetic tree
        meta_file: metadata file including patients name, id, and collection dates as a csv file 
        with columns' names 'name', 'id', and 'date'
        output_directory: the direcotry the output should be saved in.
        tree_type: type of input tree file such as 'newick', 'nexus', 'phyloxml'.
        mu: tranmsmission rate.
        cons_type: the likelihood used to build a consensus transmission tree: 'joint' (default),
        'phylogenetic', 'network'
        distr_type: type of random contact network degree distribution: 'power law', 'custom'
        deg_distr: exponent of the power law contact network degree distribution 
        (if degDistrType = 'power law') or actual distribution (if degDistrType = 'custom')
        enforce_tree: whether to force algorithm to consider only sampled acyclic transmission networks: 0, 1
        processes: number of processes to run in parallel, if 0 or 1 it runs in non-parallel.
    """
        
    tree_file = sys.argv[1]
    meta_file = sys.argv[2]
    output_directory = sys.argv[3]
    tree_type = 'newick'
    mu = 0.05
    iterations = int(sys.argv[4])
    cons_type = 'joint'
    distr_type = 'power law'
    deg_distr = 2.5 
    enforce_tree = 0
    processes = int(sys.argv[5])
    trans_net = infer_trans_network(tree_file, meta_file, tree_type)
    
