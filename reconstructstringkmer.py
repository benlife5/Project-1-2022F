import collections
import matplotlib.pyplot as plt
import numpy as np
import utils
from copy import deepcopy
from tqdm import tqdm

sequence_reads, qualities = utils.read_fastq('TeleTubby.fastq')
# find all unique 7-mer nodes
k7mernodes = []
edges = {}
reverse_edges = {}
for kmer in sequence_reads:
    prefix = kmer[:7]
    suffix = kmer[1:]
    if prefix not in k7mernodes:
        k7mernodes.append(prefix)
    if suffix not in k7mernodes:
        k7mernodes.append(suffix)
    if prefix not in edges:
        edges[prefix] = []
    if suffix not in reverse_edges:
        reverse_edges[suffix] = []
    # populate the adjacency list
    edges[prefix].append(suffix)
    reverse_edges[suffix].append(prefix)
# print("Hello")

start, stop = "", ""
for node in k7mernodes:
    if node not in edges:
        stop = node
    if node not in reverse_edges:
        start = node

path = []
edges_copy = deepcopy(edges)
# want to walk each edge once, we created one edge for each read
current_edge = start
current_path = []
branched_index = 0
while len(path) < len(sequence_reads):
    current_path.append(current_edge)

    if current_edge != stop and len(edges_copy[current_edge]) > 0:
        current_edge = edges_copy[current_edge].pop(0)
    else:
        path = path[:branched_index - 1] + current_path + path[branched_index:]
        current_path = []
        new_start_edge = ""
        for edge in edges_copy:
            if len(edges_copy[edge]) > 0 and edge in path:
                # print(edge)
                new_start_edge = edge
                break
        if new_start_edge == "":
            break
        # print(edges_copy[new_start_edge], new_start_edge)
        branched_index = path.index(new_start_edge)
        current_edge = new_start_edge

final_path = path[0]
for val in path[1:]:
    final_path += val[-1]
print(final_path)
print(len(final_path))