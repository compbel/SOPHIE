#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 15:34:18 2022

@author: fatemehmohebbi
"""


import networkx as nx
import numpy as np
import pandas as pd
import random, math, pickle
from scipy.linalg import expm
from collections import Counter
import multiprocessing as mp
from likelihood import condition_likelihood
from functools import partial


def get_clades(G, labels, none_label):
    """none_label is None in list"""
    labels_count = []
    nodes = list(G.nodes)
    for i in nodes:
        if G.out_degree(i) == 0:
            labels_count.append(pd.DataFrame(np.zeros(len(set(labels))), index = list(set(labels))))
        else:
            descendants = nx.descendants(G, i)
            labels_ = np.array(labels)[list(descendants)]
            labels_count.append(pd.DataFrame.from_dict(Counter(labels_), orient='index'))
    clades = pd.concat(labels_count, axis=1)
    clades = clades.fillna(0)
    clades.columns = nodes
    clades = clades.T.drop(columns=[none_label])
    clades = clades.reindex(sorted(clades.columns), axis=1)
    return clades


def felsenstein(iterations, patients, weight_mat_reduced, mu):
    """implementation of felsenstein pruning algorithm, which samples labels for 
    tree nodes"""
    tree = nx.from_numpy_matrix(np.triu(weight_mat_reduced), create_using=nx.DiGraph)
    patients_list = sorted(list(set(patients) - {'None'}))
    patients_list.append('None')
    patients_index = [patients_list.index(i) for i in patients]
    all_patients = len(set(patients)) - 1
    prob_distr = [1/all_patients] * all_patients
    Q = (mu/all_patients) * np.ones((all_patients, all_patients))
    Q[np.diag_indices_from(Q)] = -(all_patients - 1) * (mu/all_patients)
    likelihood = condition_likelihood(patients, weight_mat_reduced, Q)
    
    root = list(nx.topological_sort(tree))[0]      
    bfs_order = [v for u, v in nx.bfs_edges(tree, root)]

    sampled_likelihood = [0] * iterations
    all_sampled_patient = []
    clades = get_clades(tree, patients_index, patients_list.index('None'))
    P_ = calculate_P(tree, weight_mat_reduced, Q)
    
    for iteration in range(iterations):
        print(iteration)
        parents_tnet = np.full(len(patients_list) - 1, -1)
        sampled_patient = patients_index[:]
        PP = prob_distr * likelihood[root,:] * clades.loc[root]
        CDF = pd.DataFrame(PP).div(sum(PP)).cumsum() >= random.uniform(0, 1)
        sampled_patient[root] = np.where(CDF == True)[0][0]
        sampled_likelihood[iteration] = sampled_likelihood[iteration] + math.log(prob_distr[sampled_patient[root]])
            
        for i in bfs_order:
            parent = list(tree.predecessors(i))[0]
            P = P_[i]
            if tree.out_degree(i) != 0:
                clade_i = clades.loc[i]
                PP = P[sampled_patient[parent],:] * likelihood[i,:] 
                no_samp = (parents_tnet > -1) & (parents_tnet != sampled_patient[parent])
                no_samp[sampled_patient[parent]] = 0
                PP = PP * np.logical_not(no_samp)
                if sum(PP) ==0:
                    sampled_likelihood[iteration] = -math.inf
                    break
                clade_i.loc[sampled_patient[parent]] = max(1, clade_i.loc[sampled_patient[parent]])
                
                PP = PP * clade_i
                CDF = pd.DataFrame(PP).div(sum(PP)).cumsum() >= random.uniform(0, 1)
                sampled_patient[i] = np.where(CDF == True)[0][0]
                if (sampled_patient[i] != sampled_patient[parent]):
                    parents_tnet[sampled_patient[i]] = sampled_patient[parent]
            sampled_likelihood[iteration] = sampled_likelihood[iteration] + \
                    math.log(P[sampled_patient[parent],sampled_patient[i]])
        all_sampled_patient.append(sampled_patient)
        
    isfinite = np.where(np.isfinite(sampled_likelihood))[0]
    sampled_likelihood = np.array(sampled_likelihood)[isfinite]
    all_sampled_patient = np.array(all_sampled_patient)[isfinite]
    
    all_sampled_patient, unique_indx = np.unique(all_sampled_patient, axis=0, return_index=True)
    sampled_likelihood = np.array(sampled_likelihood)[unique_indx]
    
    all_sampled_patient = pd.DataFrame(all_sampled_patient).T
    return all_sampled_patient.to_numpy(), sampled_likelihood


def calculate_P(tree, weight_mat_reduced, Q):
    root = list(nx.topological_sort(tree))[0]
    bfs_order = [v for u, v in nx.bfs_edges(tree, root)]
    P = [0] * len(tree.nodes)
    for i in bfs_order:
        parent = list(tree.predecessors(i))[0]
        P[i] = expm(weight_mat_reduced[parent, i] * Q)
    return P


def felsenstein_func(iteration, likelihood, P_, tree, clades, prob_distr, patients_list, patient_index):
    # print(iteration)
    
    parents_tnet = np.full(len(patients_list) - 1, -1)
    sampled_patient = patient_index[:]
    root = list(nx.topological_sort(tree))[0]      
    bfs_order = [v for u, v in nx.bfs_edges(tree, root)]
    rand = np.random.rand(len(bfs_order))
    PP = prob_distr * likelihood[root,:] * clades.loc[root]
    CDF = pd.DataFrame(PP).div(sum(PP)).cumsum() >= random.uniform(0, 1)
    sampled_patient[root] = np.where(CDF == True)[0][0]
    sampled_likelihood =  math.log(prob_distr[sampled_patient[root]])
    
    for i in bfs_order:
        parent = list(tree.predecessors(i))[0]
        P = P_[i]
        if tree.out_degree(i) != 0:
            clade_i = clades.loc[i]
            PP = P[sampled_patient[parent],:] * likelihood[i,:] 
            no_samp = (parents_tnet > -1) & (parents_tnet != sampled_patient[parent])
            no_samp[sampled_patient[parent]] = 0
            PP = PP * np.logical_not(no_samp)
            if sum(PP) ==0:
                sampled_likelihood = -math.inf
                break
            clade_i.loc[sampled_patient[parent]] = max(1, clade_i.loc[sampled_patient[parent]])
            
            PP = PP * clade_i
            CDF = pd.DataFrame(PP).div(sum(PP)).cumsum() >= rand[i] #random.uniform(0, 1)
            sampled_patient[i] = np.where(CDF == True)[0][0]
            if (sampled_patient[i] != sampled_patient[parent]):
                parents_tnet[sampled_patient[i]] = sampled_patient[parent]
        sampled_likelihood = sampled_likelihood + math.log(P[sampled_patient[parent],sampled_patient[i]])
        
    return sampled_patient, sampled_likelihood


def felsenstein_parallel(iterations, patients, weight_mat_reduced, mu, processes):
    """implementation of a parallel felsenstein pruning algorithm, which samples labels for 
    tree nodes"""
    print('Felsenstein with ', iterations, ' number of iterations ...')
    a_pool = mp.Pool(processes = processes)
    all_sampled_patient = []
    all_sampled_likelihood = []
    all_patients = len(set(patients)) - 1
    Q = (mu/all_patients) * np.ones((all_patients, all_patients))
    Q[np.diag_indices_from(Q)] = -(all_patients - 1) * (mu/all_patients)
    
    likelihood = condition_likelihood(patients, weight_mat_reduced, Q)
    patients_list = sorted(list(set(patients) - {'None'}))
    patients_list.append('None')
    patient_index = [patients_list.index(i) for i in patients]
    all_patients = len(set(patients)) - 1
    
    prob_distr = [1/all_patients] * all_patients
    parents_tnet = np.full(all_patients, -1)
    tree = nx.from_numpy_matrix(np.triu(weight_mat_reduced), create_using=nx.DiGraph)
    clades = get_clades(tree, patient_index, patients_list.index('None'))
    P = calculate_P(tree, weight_mat_reduced, Q)
    
    partial_function = partial(felsenstein_func, likelihood=likelihood, P_=P,\
                               tree=tree, clades=clades, prob_distr=prob_distr, patients_list=patients_list, patient_index=patient_index)
    for samples, sampled_likelihood in a_pool.map(partial_function, range(iterations)):
        all_sampled_patient.append(samples)
        all_sampled_likelihood.append(sampled_likelihood)
    print(len(all_sampled_patient))
        
    isfinite = np.where(np.isfinite(all_sampled_likelihood))[0]
    all_sampled_likelihood = np.array(all_sampled_likelihood)[isfinite]
    all_sampled_patient = np.array(all_sampled_patient)[isfinite]
    print(all_sampled_patient.shape)
    
    all_sampled_patient, unique_indx = np.unique(all_sampled_patient, axis=0, return_index=True)
    all_sampled_likelihood = np.array(all_sampled_likelihood)[unique_indx]
    print(all_sampled_patient.shape)
    
    all_sampled_patient = pd.DataFrame(all_sampled_patient).T
    
    if len(all_sampled_likelihood) == 0:
        raise ValueError("No network has been sampled.")
    
    return all_sampled_patient.to_numpy(), all_sampled_likelihood


