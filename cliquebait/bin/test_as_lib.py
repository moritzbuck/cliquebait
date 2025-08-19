#!/usr/bin/env python3

import os
import sys

import argparse
import os
import sys
from sys import stderr
from cliquebait.utils import parse_fastani_output
from cliquebait.guidetrees import available_guidetrees, get_guidetree_class
from cliquebait.clustering.cliquebaitClustering import cliqueblocksClustering
import cliquebait
from cliquebait import default_gap_size, default_strain_cutoff, default_lower_cutoff, default_denoising_cutoff, default_guidetree, default_min_size

description_text = "TO DO"

anis = parse_fastani_output("/data/moritz/0079_pelaginet/fastani_bacillus.all_v_all.tsv")

clustering = cliqueblocksClustering(default_guidetree, anis, gap_size=default_gap_size, strain_cutoff=default_strain_cutoff, bottom_cutoff=default_lower_cutoff, denoising_cutoff=default_denoising_cutoff, size_cutoff=default_min_size)
clustering.cluster_simple()


clustering.guide_tree.draw_dendrogram(clusters=clustering.final_clusters)

clustered = frozenset.union(*clustering.final_clusters)

unclustered_nodes = [n for n in clustering.guide_tree.iterate_nodes() if len(clustered.intersection(clustering.guide_tree.get_node_info(n)['genomes'])) == 0]

key_node = [(n, len(clustering.guide_tree.get_node_info(n)['genomes'])) for n in unclustered_nodes]

top_node = max(key_node, key =  lambda x : x[1])[0]

missed_cluster = [frozenset(clustering.guide_tree.get_node_info(top_node)['genomes'])]

clustering.guide_tree.draw_dendrogram(clusters=missed_cluster)

tree_stats = {v['id'] : v for v in clustering.guide_tree.get_nodes_info()}
test_genomes = list(missed_cluster[0])[350]

node = [v for v in tree_stats.values() if test_genomes in v['genomes'] and len(v['genomes']) == 1][0]
