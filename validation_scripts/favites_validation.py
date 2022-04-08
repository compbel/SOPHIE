import os
import sys
from os import listdir
from os.path import join, isfile
from ete3 import Tree
import re

metrics = ['sensitivity', 'specificity', 'f1']
tools = ['tnet', 'phylo']


def process_folder_phyloscanner(folder, suffix, tree_file_name):
    t = Tree(join(folder, 'error_free_files', tree_file_name))
    for node in t.traverse("postorder"):
        if node.name.startswith("N"):
            sp = node.name.split("|")
            host_id = sp[1]
            sequence_id = sp[0][1:]
            node.name = "H{}H_{}".format(host_id, sequence_id)

    output = open(join(folder, 'error_free_files/formatted_best_tree_scanner_' + suffix + '.raxmltree'), 'w')
    output.write(t.write())
    output.write("\n")
    output.close()


def process_folder_tnet(folder, suffix, tree_file_name):
    t = Tree(join(folder, 'error_free_files', tree_file_name))
    for node in t.traverse("postorder"):
        if node.name.startswith("N"):
            sp = node.name.split("|")
            host_id = sp[1]
            sequence_id = sp[0][1:]
            node.name = "{}_{}".format(host_id, sequence_id)

    output = open(join(folder, 'error_free_files/formatted_best_tree_tnet_' + suffix + '.raxmltree'), 'w')
    output.write(t.write())
    output.write("\n")
    output.close()


def run_phyloscanner(folder, suffix):
    input_path = join(folder, 'error_free_files/formatted_best_tree_scanner_' + suffix + '.raxmltree')
    outfolder = join('out_' + suffix + '/', folder.split('/')[-1])
    mkdir_command = 'mkdir -p {}'
    run_command = '~/phyloscanner/phyloscanner_analyse_trees.R {} network s,0 -od {}  --overwrite --tipRegex "(H[0-9]+H)_([0-9]+)" --verbose 2 --allClassifications --collapsedTrees'
    os.system(mkdir_command.format(outfolder))
    os.system(run_command.format(input_path, outfolder))


def run_tnet(folder, suffix):
    input_path = join(folder, 'error_free_files/formatted_best_tree_tnet.raxmltree')
    mkdir_command = 'mkdir -p {}'
    outfolder = 'tnet_out_' + suffix
    output_name = outfolder + '/' + folder.split('/')[-1] + '.txt'
    run_command = '~/projects/tnet_python/tnet.py {} {}'
    os.system(mkdir_command.format(outfolder))
    os.system(run_command.format(input_path, output_name))


def read_ground_truth_network(folder):
    path = join(folder, 'error_free_files/transmission_network.txt')
    lines = open(path).readlines()
    edges = []
    sp = lines[0].split('\t')
    # skip the root since phyloscanner doesn't give the root
    # edges.append((-1, int(sp[1])))
    for line in lines[1:]:
        sp = line.split('\t')
        u = int(sp[0])
        v = int(sp[1])
        if u == v:
            continue
        edges.append((u, v))
    return edges


def read_tnet_network(folder, suffix):
    outfolder = 'tnet_out_' + suffix
    path = outfolder + '/' + folder.split('/')[-1] + '.txt'
    lines = open(path).readlines()
    edges = []
    sp = lines[0].split('\t')
    # skip the root since phyloscanner doesn't give the root
    # edges.append((-1, int(sp[1])))
    for line in lines[1:]:
        sp = line.split('\t')
        edge = (int(sp[0]), int(sp[1]))
        if edge not in edges:
            edges.append(edge)
    return edges


def read_phyloscanner_network(folder, suffix):
    path = join('out_' + suffix + '/', folder.split('/')[-1], 'network_classification.csv')
    try:
        lines = open(path).readlines()
    except FileNotFoundError:
        # print("Didn't find " + path)
        return []
    edges = []
    for line in lines[1:]:
        sp = line.split(',')
        u = int(sp[0][1:-1])
        v = int(sp[1][1:-1])
        direction = sp[8]
        adjacent = sp[2]
        if adjacent == 'TRUE' and (direction == 'anc' or direction == 'desc'):
            edge = (u, v) if direction == 'anc' else (v, u)
            if edge not in edges:
                edges.append(edge)
    return edges

# calculates the metrics for reconstructed networks by comparing given edge lists
def calculate_stat(ground_truth, inferred):
    sensitivity = len([e for e in ground_truth if e in inferred]) / len(ground_truth)
    if len(inferred) == 0:
        specificity = 0
    else:
        specificity = len([e for e in inferred if e in ground_truth]) / len(inferred)
    if sensitivity < 0.00001 or specificity < 0.00001:
        f1 = 0
    else:
        f1 = 2 * (sensitivity * specificity) / (sensitivity + specificity)
    return sensitivity, specificity, f1

# prints statistics to console and save detailed statistics by categories in results.txt
def print_results(results, suffix):
    for category in results:
        n_samples = len(results[category][tools[0]]['sensitivity'])
        print("Category {}. Total samples {}".format(category, n_samples))

        for tool in tools:
            values = [sum(results[category][tool][metric]) / n_samples for metric in metrics]
            print(tool + " ", end='')
            print(",".join(["{} = {:.3f}".format(metric, value) for (metric, value) in zip(metrics, values)]))

    out = open('results_'+suffix+'.txt', 'w')
    for category in results:
        out.write("{}\n".format(category))
        column_names = ['sample_name']
        for tool in tools:
            for metric in metrics:
                column_names.append(tool + "_" + metric)

        out.write(",".join(column_names))
        out.write("\n")
        sample_n = 0
        for sample_name in results[category]['sample_names']:
            values = [sample_name]
            for tool in tools:
                for metric in metrics:
                    values.append("{:.3f}".format(results[category][tool][metric][sample_n]))
            sample_n += 1
            out.write(",".join(values))
            out.write("\n")


def process_data(folders, suffix, tree_file_name):
    for folder in folders:
        process_folder_phyloscanner(folder, suffix, tree_file_name)
        process_folder_tnet(folder, suffix, tree_file_name)


def run_tools(folders, suffix):
    for folder in folders:
        run_phyloscanner(folder, suffix)
        run_tnet(folder, suffix)


# reads the outputs for all tools and calculate metrics with regard to the ground truth
def run_analysis(folders, suffix):
    categories = {}
    for folder in folders:
        sample_name = folder.split('/')[-1]
        m = re.search('FAVITES_output_(\w+)_T.*', sample_name)
        category = m.group(1)
        ground_truth = read_ground_truth_network(folder)
        tnet = read_tnet_network(folder, suffix)
        phylo = read_phyloscanner_network(folder, suffix)
        reconstructed_networks = {'tnet': tnet, 'phylo': phylo}
        if len(ground_truth) == 0:
            continue
        if category not in categories:
            categories[category] = {tool: {metric: [] for metric in metrics} for tool in tools}
            categories[category]['sample_names'] = []
        for tool in tools:
            values = calculate_stat(ground_truth, reconstructed_networks[tool])
            for metric, value in zip(metrics, values):
                categories[category][tool][metric].append(value)
        categories[category]['sample_names'].append(folder.split('/')[-1])
    print_results(categories, suffix)


# example to run python3 favites_validation.py Favites_inputs/ data
if __name__ == "__main__":
    # the folder with FAVITES outputs
    input_folder = sys.argv[1]
    # data - for data praparation, tools - to run the tools on prepared data, analyze - to run performance analysis on tools' outputs
    command = sys.argv[2]
    # output suffix for different inputs(for example if we want to analyze separately tree_0.time.tre and RAxML_bestTree.raxmltree)
    suffix = 'time'
    # the tree file to run inside the error_free_files folder in FAVITES output
    tree_file_path = 'phylogenetic_trees/tree_0.time.tre'
    folders = [join(input_folder, f) for f in listdir(input_folder) if not isfile(join(input_folder, f))]
    if command == 'data':
        process_data(folders, suffix, tree_file_path)
    elif command == 'tools':
        run_tools(folders, suffix)
    elif command == 'analyze':
        run_analysis(folders, suffix)
