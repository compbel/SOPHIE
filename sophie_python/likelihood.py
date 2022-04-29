#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 15:06:26 2022

@author: fatemehmohebbi
"""


import networkx as nx
import numpy as np
import pandas as pd
from scipy.linalg import expm
from collections import Counter
from scipy.special import zeta
import math
import multiprocessing as mp
from utils import log_choose_k
from scipy.optimize import linear_sum_assignment
from functools import partial


def get_netwok(sampled_label_arr, tree, patients_list):
    """convert a sampled labeling to a network"""
    network = np.zeros((len(patients_list), len(patients_list)))
    for edge in list(tree.edges):
        i = sampled_label_arr[edge[0]]
        j = sampled_label_arr[edge[1]]
        if i != j:
            network[i, j] = 1
    return network


def matching(adj_mat, degree_distribution, N, enforce_tree):
    """find the maximum matching cost for the network"""
    flag = 'real'
    num_edges = np.sum(adj_mat)
    num_nodes = len(adj_mat)
    if (enforce_tree == 1) & (num_edges != num_nodes - 1):
        matching_cost =  -math.inf
    else:
        c_total = N * degree_distribution
        c = np.round(c_total)
        if sum(c) < num_nodes:
            ind = np.where(c == 0)
            for i in ind[:num_nodes - sum(c)]:
                c[i] = 1
        delta = np.where(c == 0)[0][0] 
        M = sum(np.array(range(1, int(N))) * c_total) / 2
        w = np.zeros((num_nodes,delta))
        d = np.sum(adj_mat, axis = 0) + np.sum(adj_mat, axis = 1)
        for i in range(num_nodes):
            for j in range(delta):
                if d[i] <= (j + 1):
                    w[i,j] = d[i] * math.log(j + 1)
                else:
                    w[i,j] = -math.inf
        u = []
        aux_facilities = 0
        for j in range(delta):
            u_j = np.zeros(int(c[j]) + 1)
            for k in range(int(c[j])):
                u_j[k + 1] = log_choose_k(c[j], k + 1,flag)
            u.append(u_j)
            aux_facilities = aux_facilities + len(u_j) - 1
        
        cost_matrix = np.zeros((num_nodes, aux_facilities))
        for i in range(num_nodes):
            count = 0
            for j in range(delta):
                for k in range(len(u[j]) - 1):
                    cost_matrix[i,count] = w[i,j] + (u[j][k + 1] - u[j][k])
                    count = count + 1
        try:
            row_ind, col_ind = linear_sum_assignment(cost_matrix, maximize=True)
            matching_cost = cost_matrix[row_ind, col_ind].sum() -\
                    num_edges * math.log(2 * M) - log_choose_k(N, num_nodes, flag)
        except:
            matching_cost = -math.inf
    return matching_cost


def network_likelihood(eg_distr_type, deg_distr, sampled_label, weight_mat_reduced, patients_list, enforce_tree):
    print('Networks likelihood calculation ...')
    tree = nx.from_numpy_matrix(np.triu(weight_mat_reduced), create_using=nx.DiGraph)
    if eg_distr_type == 'power law':
        d = deg_distr
        r_zeta = zeta(d)
    num_of_patients = len(patients_list)
    N = []
    degree_distr_C = [0] * sampled_label.shape[1]
    network_likel = []
    networks = []
    for i in range(sampled_label.shape[1]):
        network = get_netwok(sampled_label[:,i], tree, patients_list)
        degree_counts = pd.DataFrame(np.zeros(num_of_patients), index = patients_list)
        if eg_distr_type == 'power law':
            degrees = np.sum(network, axis=1)
            counts = Counter(degrees)
            for k in counts.keys():
                degree_counts.iloc[int(k)] = counts[k]
            N.append(max(np.ceil(r_zeta * degree_counts[0].values * np.array(patients_list) ** d)))
            degree_distr_C[i] = np.array(range(1, int(N[i]))) ** -d / r_zeta
        if eg_distr_type == 'costum':
             degree_distr_C[i] = deg_distr
        network_likel.append(matching(network, degree_distr_C[i], N[i], enforce_tree))
        networks.append(network)
    return network_likel, networks


def network_function(i, eg_distr_type, deg_distr, all_sampled_label, tree, patients_list, enforce_tree):
    sampled_label = all_sampled_label[:,i]
    if eg_distr_type == 'power law':
        d = deg_distr
        r_zeta = zeta(d)
    N = []
    degree_distr_C = [0] * len(sampled_label)
    network = get_netwok(sampled_label, tree, patients_list)
    degree_counts = pd.DataFrame(np.zeros(len(patients_list)), index = patients_list)
    if eg_distr_type == 'power law':
        degrees = np.sum(network, axis=1)
        counts = Counter(degrees)
        for k in counts.keys():
            degree_counts.iloc[int(k)] = counts[k]
        N = max(np.ceil(r_zeta * degree_counts[0].values * np.array(patients_list) ** d))
        degree_distr_C = np.array(range(1, int(N))) ** -d / r_zeta
    if eg_distr_type == 'costum':
        degree_distr_C = deg_distr
    network_likel = matching(network, degree_distr_C, N, enforce_tree)
    return network_likel, network


def network_likelihood_parallel(eg_distr_type, deg_distr, all_sampled_label, weight_mat_reduced, patients_list, enforce_tree, processes):
    print('Networks likelihood calculation ...')
    a_pool = mp.Pool(processes = processes)
    tree = nx.from_numpy_matrix(np.triu(weight_mat_reduced), create_using=nx.DiGraph)
    network_likel = []
    networks = []

    pool_partial = partial(network_function, eg_distr_type = eg_distr_type, deg_distr = deg_distr,all_sampled_label = all_sampled_label,\
                           tree = tree, patients_list = patients_list, enforce_tree = enforce_tree)
    for net_likle, network in a_pool.map(pool_partial, range(all_sampled_label.shape[1])):
        network_likel.append(net_likle)
        networks.append(network)
    return network_likel, networks


def condition_likelihood(patients, weight_mat_reduced, Q):
    """calculates conditiona likelihood where p=exp(t*Q) and t is (distance) weight"""
    print('likelihood calculation ...')
    min_prob = 1.0e-10
    tree = nx.from_numpy_matrix(np.triu(weight_mat_reduced), create_using=nx.DiGraph)
    patients_list = sorted(list(set(patients)))
    patients_index = [patients_list.index(i) for i in patients]
    L = np.zeros((len(list(tree.nodes)), len(patients_list) - 1))
    for i in list(tree.in_degree):
        if i[1] == 0:
            root = i[0]
    dfs_order = list(nx.dfs_preorder_nodes(tree, source=root))
    dfs_order.reverse()
    for i in dfs_order:
        # print(i)
        if tree.out_degree(i) == 0:
            L[i][patients_index[i]] = 1
        else:
            child = list(tree.successors(i))
            if (weight_mat_reduced[i,child[0]] == 0) & (weight_mat_reduced[i,child[1]] == 0):
                L[i,:] = (L[child[1],:] + L[child[2],:]) / 2
                continue
            P_1 = expm(weight_mat_reduced[i,child[0]] * Q)
            P_2 = expm(weight_mat_reduced[i,child[1]] * Q)
            for j in range(len(patients_list) - 1): 
                L_1 = sum(P_1[j,:] * L[child[0],:])
                L_2 = sum(P_2[j,:] * L[child[1],:])
                L[i,j] = L_1 * L_2
        min_L = min(L[i,:][np.nonzero(L[i,:])])
        if min_L < min_prob:
            sc_factor = min_prob / min_L
            L[i,:] = L[i,:] * sc_factor
            sc_factor_1 = sum(L[i,:])
            L[i,:] = L[i,:] / sc_factor_1
    return L
