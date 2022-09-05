import collections
import matplotlib.pyplot as plt
import numpy as np
import utils
from copy import deepcopy


# from tqdm import tqdm

# sequence_reads, qualities = utils.read_fastq('TeleTubby.fastq')
# # find all unique 7-mer nodes
# k7mernodes = []
# edges = {}
# reverse_edges = {}
# for kmer in sequence_reads:
#     prefix = kmer[:7]
#     suffix = kmer[1:]
#     if prefix not in k7mernodes:
#         k7mernodes.append(prefix)
#     if suffix not in k7mernodes:
#         k7mernodes.append(suffix)
#     if prefix not in edges:
#         edges[prefix] = []
#     if suffix not in reverse_edges:
#         reverse_edges[suffix] = []
#     # populate the adjacency list
#     edges[prefix].append(suffix)
#     reverse_edges[suffix].append(prefix)
# # print("Hello")
#
# start, stop = "", ""
# for node in k7mernodes:
#     if node not in edges:
#         stop = node
#     if node not in reverse_edges:
#         start = node
#
# path = []
# edges_copy = deepcopy(edges)
# # want to walk each edge once, we created one edge for each read
# current_edge = start
# current_path = []
# branched_index = 0
# while len(path) < len(sequence_reads):
#     current_path.append(current_edge)
#
#     if current_edge != stop and len(edges_copy[current_edge]) > 0:
#         current_edge = edges_copy[current_edge].pop(0)
#     else:
#         path = path[:branched_index - 1] + current_path + path[branched_index:]
#         current_path = []
#         new_start_edge = ""
#         for edge in edges_copy:
#             if len(edges_copy[edge]) > 0 and edge in path:
#                 # print(edge)
#                 new_start_edge = edge
#                 break
#         if new_start_edge == "":
#             break
#         # print(edges_copy[new_start_edge], new_start_edge)
#         branched_index = path.index(new_start_edge)
#         current_edge = new_start_edge
#
# final_path = path[0]
# for val in path[1:]:
#     final_path += val[-1]
# print(final_path)
# print(len(final_path))


def breakkmer(k, read):
    rtn_list = []
    last_index = len(read) - 1
    i = k - 1
    j = 0
    while i + k <= last_index:
        rtn_list.append(read[j:i])
        j += k
        i += k
    return rtn_list


sequence_reads_real, qualities_real = utils.read_fastq('ABS2-LN-R1_cleaned_paired.fastq.gz')
kmer_list = []
qual_list = []

# break each kmer down into a uniform, smaller size and add to final list of values
for index, qual_read in enumerate(qualities_real):
    temp1 = breakkmer(10, sequence_reads_real[index])
    temp2 = breakkmer(10, qual_read)
    for val in range(len(temp2)):
        qual_list.append(temp2[val])
        kmer_list.append(temp1[val])

print(range(len(qualities_real)))
print(len(sequence_reads_real[0]))
best_sequence = ""
best_score = 0

avg_len = 0

for val in kmer_list:
    avg_len += len(val)
avg_len = avg_len / len(kmer_list)
print("avg_len is ", avg_len)

print(len(kmer_list), "is length before")

# todo
for index, qual in enumerate(kmer_list):
    score = 0
    for char in qual_list[index]:
        score += (ord(char) - 33)
    score_avg = score / len(qual_list[index])
    if score_avg <= 30:
        kmer_list.pop(index)
        qual_list.pop(index)
    #     print(score, "for this sequence", qual)
    # if score > best_score:
    #     best_score = score
    #     best_sequence = qual


# print(score, "average score is ", score / len(kmer_list))
# print(len(kmer_list), "is length after")


def build_graph(sequence_reads):
    k7mernodes = []
    edges = {}
    reverse_edges = {}
    for kmer in sequence_reads:
        prefix = kmer[:-1]
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
    return k7mernodes, edges, reverse_edges


nodes_tt, edges_tt, reverse_edges_tt = build_graph(kmer_list)


def assemble_sequence_algo(nodes, edges, reverse_edges):
    start, stop = "", ""
    for node in nodes:
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
        #         print(current_edge)
        current_path.append(current_edge)

        if current_edge in edges_copy and len(edges_copy[current_edge]) > 0:
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

    #     print('path', path)
    final_path = path[0]
    for val in path[1:]:
        final_path += val[-1]
    return final_path


path = assemble_sequence_algo(nodes_tt, edges_tt, reverse_edges_tt)

print(path)
