# SOPHIE
SOcial and PHilogenetic Investigation of Epidemics

## Tool
Main function: AMNetCons = inferTransNetSOPHIE(filePhylo,nSamp,degDistrType,degDistr,transRate,consType,enforceTree,visualize)

Input paramers:
- filePhylo - name of file with the phylogenetic tree (in Newick or Nexus formats)
- nSamp - number of samples from the joint distribution of ancestral label assignments
- degDistrType - type of random contact network degree distribution.
                possible variants: 'power law', 'custom'
- degDistr - exponent of the power law contact network degree distribution
        (if degDistrType = 'power law') or actual distribution (if degDistrType = 'custom')
- transRate - transmission rate (used in Jukes-Cantor trait substitution model)
- consType - the likelihood used to build a consensus transmission
                tree. Variants: 'joint' (default), 'phylogenetic', 'network'
- enforceTree - whether to force algorithm to consider only sampled
               acyclic transmission networks. Variants: 0 (default), 1
- visualize -  whether to output a figure with the inferred transmission
             network. Default value: 0

Output:

AMNetCons - adjacency matrix of the inferred transmission network

The example of tool execution can be found in exampleRun.m


Implementation of Edmonds algorithm for maximal arborescence is taken
from https://www.mathworks.com/matlabcentral/fileexchange/24899-edmonds-algorithm 
Copyright (c) 2009, Ashish Choudhary. All rights reserved.

# DOI
[![DOI](https://zenodo.org/badge/479482812.svg)](https://zenodo.org/badge/latestdoi/479482812)

# Python version of SOPHIE 
You can find the python version under "sophie_python" folder. 
