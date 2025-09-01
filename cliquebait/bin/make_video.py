#!/usr/bin/env python3

import os
import sys
import pandas
import argparse
import os
import sys
from sys import stderr
from cliquebait.utils import parse_fastani_output, checkm_parser
from cliquebait.guidetrees import available_guidetrees, get_guidetree_class
from cliquebait.clustering.cliquebaitClustering import cliqueblocksClustering
import cliquebait
from cliquebait import default_gap_size, default_strain_cutoff, default_lower_cutoff, default_denoising_cutoff, default_guidetree, default_min_size
from plotnine import * 


anis = parse_fastani_output("/data/moritz/0079_pelaginet/fastani.all_v_all.tsv")
checkm_genomes = checkm_parser("/data/moritz/0079_pelaginet/checkm.txt", 70, 5)

anis.filter_genomes(checkm_genomes)

clustering = cliqueblocksClustering(default_guidetree, anis, size_cutoff=6)
clustering.cluster_simple()


# position map
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import networkx
import matplotlib

plt.clf()
matplotlib.rcParams['figure.figsize'] = 400,400

edges = { frozenset((g1,g2))  for cluster in clustering.final_clusters for g1 in cluster for g2 in cluster if g1 != g2}
genomes = list(frozenset.union(*edges))
graph = networkx.Graph(list(edges))

pos = networkx.spring_layout(graph, method = "energy", gravity = 20)

data = clustering.get_clusters_stats()
classes = {gg : dd['id'] for dd in  data for gg in dd['genomes']}
communities =  [classes.get(g, 0) for g in graph]
numCommunities = max(communities)+1
# color map from http://colorbrewer2.org/
cmapLight = colors.ListedColormap(['#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f', '#cab2d6'], 'indexed', numCommunities)
cmapDark = colors.ListedColormap(['#1f78b4', '#33a02c', '#e31a1c', '#ff7f00', '#6a3d9a'], 'indexed', numCommunities)


from tqdm import tqdm
for i in tqdm(list(range(100, 89, -1))):
    plt.clf()
    edges = {k for k,v in clustering.ani_dictionary.items() if v > i and all([g in genomes for g in k])}
    graph = networkx.Graph(list(edges))
    communities =  [classes.get(g, 0) for g in graph]
    nodeCollection = networkx.draw_networkx_nodes(graph,
                pos = pos,
                node_color = communities,
                cmap = cmapLight
        )
    networkx.draw_networkx_edges(graph, pos)
    plt.axis('off')
    plt.savefig(f"network_{i:03d}.png")
