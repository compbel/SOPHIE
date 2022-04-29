
A python implementation of SOPHIE.
The main script is 'sophie.py'.

inputs are:
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
       Iterations: number of iterations.

How to run?
command line: "python3.8 sophie.py tree_0.time.tre tree_test_metadata.csv output_directory Iterations processes"
