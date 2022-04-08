import os
import re
import sys
from os import listdir
from os.path import join, isfile

from ete3 import Tree

metrics = ['sensitivity', 'specificity', 'f1', 'tp']
tools = ['tnet', 'phylo']
output_folder = 'real_formatted_HCV'
sources = {'AA': 45, 'AC': 124, 'AI': 4, 'AJ': 199, 'AQ': 89, 'AW': 2, 'BA': 3, 'BB': 45, 'BC': 46, 'BJ': 28}


def process_folder_phyloscanner(folder, suffix, tree_file_name):
    f = open(join(folder, tree_file_name))
    content = f.read()
    content = re.sub(r'Branch\s[0-9]+', '', content)
    content = content.replace('Root', '')
    t = Tree(content)
    for node in t.traverse("postorder"):
        if node.name.startswith("A") or node.name.startswith('B'):
            m = re.search(suffix + '([0-9]+)_.+_[0-9]+_([0-9]+)_[0-9]+', node.name)
            host_id = m.group(1)
            sequence_id = m.group(2)
            node.name = "H{}H_{}".format(host_id, sequence_id)

    output = open(join(output_folder, 'formatted_best_tree_scanner_' + suffix + '.raxmltree'), 'w')
    output.write(t.write())
    output.write("\n")
    output.close()


def process_folder_tnet(folder, suffix, tree_file_name):
    f = open(join(folder, tree_file_name))
    content = f.read()
    content = re.sub(r'Branch\s[0-9]+', '', content)
    content = content.replace('Root', '')
    t = Tree(content)
    for node in t.traverse("postorder"):
        if node.name.startswith("A") or node.name.startswith('B'):
            m = re.search(suffix + '([0-9]+)_.+_[0-9]+_([0-9]+)_[0-9]+', node.name)
            host_id = m.group(1)
            sequence_id = m.group(2)
            node.name = "{}_{}".format(host_id, sequence_id)

    output = open(join(output_folder, 'formatted_best_tree_tnet_' + suffix + '.raxmltree'), 'w')
    output.write(t.write())
    output.write("\n")
    output.close()


def run_phyloscanner(folder, suffix, sample):
    input_path = join(folder, 'formatted_best_tree_scanner_' + sample + '.raxmltree')
    outfolder = join('out_' + suffix + '/', sample)
    mkdir_command = 'mkdir -p {}'
    run_command = '~/phyloscanner/phyloscanner_analyse_trees.R {} network s,0 -od {}  --overwrite --tipRegex "(H[0-9]+H)_([0-9]+)" --verbose 2 --allClassifications --collapsedTrees'
    os.system(mkdir_command.format(outfolder))
    # print(run_command.format(input_path, outfolder))
    os.system(run_command.format(input_path, outfolder))


def run_tnet(folder, suffix, sample):
    input_path = join(folder, 'formatted_best_tree_tnet_' + sample + '.raxmltree')
    mkdir_command = 'mkdir -p {}'
    outfolder = 'tnet_out_' + suffix
    output_name = outfolder + '/' + sample + '.txt'
    run_command = '~/projects/tnet_python/tnet.py {} {}'
    os.system(mkdir_command.format(outfolder))
    # print(run_command.format(input_path, output_name))
    os.system(run_command.format(input_path, output_name))


def read_ground_truth_network(folder, sample_name, tree_file_name, source):
    f = open(join(folder, tree_file_name))
    content = f.read()
    content = re.sub(r'Branch\s[0-9]+', '', content)
    content = content.replace('Root', '')
    t = Tree(content)
    edges = []
    hosts = set()
    for node in t.traverse("postorder"):
        if node.name.startswith("A") or node.name.startswith('B'):
            m = re.search(sample_name + '([0-9]+)_.+_[0-9]+_([0-9]+)_[0-9]+', node.name)
            host_id = int(m.group(1))
            if host_id != source:
                hosts.add(host_id)
    for id in hosts:
        edges.append((source, id))
    return edges


def read_tnet_network(folder, suffix, sample):
    outfolder = 'tnet_out_' + suffix
    path = outfolder + '/' + sample + '.txt'
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


def read_phyloscanner_network(folder, suffix, sample):
    path = join('out_' + suffix + '/', sample, 'network_classification.csv')
    try:
        lines = open(path).readlines()
    except FileNotFoundError:
        print("Didn't find " + path)
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
    tp = len([e for e in ground_truth if e in inferred])
    sensitivity = len([e for e in ground_truth if e in inferred]) / len(ground_truth)
    if len(inferred) == 0:
        specificity = 0
    else:
        specificity = len([e for e in inferred if e in ground_truth]) / len(inferred)
    if sensitivity < 0.00001 or specificity < 0.00001:
        f1 = 0
    else:
        f1 = 2 * (sensitivity * specificity) / (sensitivity + specificity)
    return sensitivity, specificity, f1, tp


# prints statistics to console and save detailed statistics by categories in results.txt
def print_results(results, suffix):
    for category in results:
        n_samples = len(results[category][tools[0]]['sensitivity'])
        print("Category {}. Total samples {}".format(category, n_samples))
        total_ground_truth_edges = sum(results[category]['total_gt'])

        for tool in tools:
            values = [sum(results[category][tool][metric]) / n_samples for metric in metrics]
            total_tp = sum(results[category][tool]['tp'])
            total_inferred = sum(results[category][tool]['total_inferred'])
            t_sensitivity = total_tp/total_ground_truth_edges
            t_specificity = total_tp/total_inferred
            print("total_gt {} total_inferred {} total_tp {} ".format(total_ground_truth_edges, total_inferred, total_tp))
            print(tool + " ", end='')
            print(",".join(["{} = {:.3f}".format(metric, value) for (metric, value) in zip(metrics, values)]), end='')
            print(', t_sens = {:.3f}, t_spec = {:.3f}'.format(t_sensitivity, t_specificity))

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


def process_data(folder, files):
    for file in files:
        process_folder_phyloscanner(folder, file[-2:], file)
        process_folder_tnet(folder, file[-2:], file)


def run_tools(folder, files, suffix):
    for file in files:
        run_phyloscanner(folder, suffix, file[-2:])
        run_tnet(folder, suffix, file[-2:])


# reads the outputs for all tools and calculate metrics with regard to the ground truth
def run_analysis(input_folder, output_folder, files, suffix):
    categories = {}
    for file in files:
        sample_name = file[-2:]
        category = 'real'
        ground_truth = read_ground_truth_network(input_folder, sample_name, file, sources[sample_name])
        tnet = read_tnet_network(output_folder, suffix, sample_name)
        phylo = read_phyloscanner_network(output_folder, suffix, sample_name)
        reconstructed_networks = {'tnet': tnet, 'phylo': phylo}
        if len(ground_truth) == 0:
            continue
        if category not in categories:
            categories[category] = {tool: {metric: [] for metric in metrics} for tool in tools}
            for tool in tools:
                categories[category][tool] ={metric: [] for metric in metrics}
                categories[category][tool]['total_inferred'] = []
            categories[category]['sample_names'] = []
            categories[category]['total_gt'] = []
        for tool in tools:
            values = calculate_stat(ground_truth, reconstructed_networks[tool])
            categories[category][tool]['total_inferred'].append(len(reconstructed_networks[tool]))
            for metric, value in zip(metrics, values):
                categories[category][tool][metric].append(value)
        categories[category]['sample_names'].append(sample_name)
        categories[category]['total_gt'].append(len(ground_truth))

    print_results(categories, suffix)


if __name__ == "__main__":
    input_folder = sys.argv[1]
    command = sys.argv[2]
    suffix = 'real_HCV'
    files = [f for f in listdir(input_folder) if isfile(join(input_folder, f))]
    os.system("mkdir " + output_folder)
    if command == 'data':
        process_data(input_folder, files)
    elif command == 'tools':
        run_tools(output_folder, files, suffix)
    elif command == 'analyze':
        run_analysis(input_folder, output_folder, files, suffix)
