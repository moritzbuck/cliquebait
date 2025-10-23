#!/usr/bin/env python3

import os
import sys
import pandas
import argparse
import os
import sys
from sys import stderr

import tqdm
from cliquebait.utils import parse_fastani_output, checkm_parser
from cliquebait.guidetrees import available_guidetrees, get_guidetree_class
from cliquebait.clustering.cliquebaitClustering import cliqueblocksClustering
import cliquebait
from cliquebait import default_gap_size, default_strain_cutoff, default_lower_cutoff, default_denoising_cutoff, default_guidetree, default_min_size
from plotnine import * 
import matplotlib.pyplot as plt
from plotnine import ggplot, aes, geom_histogram, facet_wrap, geom_vline
from tqdm import tqdm
import matplotlib
import networkx
from matplotlib import colors, cm

with open("/home/moritz/data/0064_bis/metadata/bac120_taxonomy_r226.tsv")  as handle:
    tax = {l.split()[0] : l.strip().split("\t")[1] for l in tqdm(handle.readlines())}

genuses = {k[3:] : v.split(";")[-2] for k,v in tax.items() }
species = {k[3:] : v.split(";")[-1] for k,v in tax.items() }


def get_ani_data(ani_file, checkm = False):
    bac_anis = parse_fastani_output(ani_file)
    if checkm:
        checkm_genomes = checkm_parser(checkm, 70, 5)
        bac_anis.filter_genomes(checkm_genomes)

    ani_df = pandas.DataFrame.from_records([ { f"g{i}" : "_".join(kk.split("/")[-1].split("_")[:2]) for i,kk in enumerate(k) } | {'ani' : v} for k,v in bac_anis.store.items()])
    ani_df['genus0'] = ani_df['g0'].map(genuses)
    ani_df['species0'] = ani_df['g0'].map(species)
    ani_df['genus1'] = ani_df['g1'].map(genuses)
    ani_df['species1'] = ani_df['g1'].map(species)
    ani_df['same_genus'] = ani_df['genus0'] == ani_df['genus1']
    ani_df['same_species'] = ani_df['species0'] == ani_df['species1']
    return ani_df

def simple_ani_graph(anis, cutoff = 95, pos = None, cols = False, file = None):
    plt.clf()
    matplotlib.rcParams['figure.figsize'] = 400,400

    edges = { frozenset((dd[1]['g0'], dd[1]['g1']))  for dd in anis.iterrows() if dd[1]['ani'] > cutoff }

    genomes = list(frozenset.union(*edges))
    graph = networkx.Graph(list(edges))

    if pos is None : 
        pos = networkx.spring_layout(graph, method = "energy", gravity = 3)

    if cols :
        spec = [species[g] for g in graph]
        spec2index = {s:i for i,s in enumerate(set(spec))}
        cmapLight  = colors.LinearSegmentedColormap.from_list('custom pastel', ['#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f', '#cab2d6'], N=max(spec2index.values())+1)

        node_cols = [cmapLight(spec2index[s]) for s in spec]
    else :
        node_cols = ["#000000" for _ in graph]

    networkx.draw_networkx_nodes(graph, pos, node_color=node_cols, node_size=50)
    networkx.draw_networkx_edges(graph, pos)
    plt.axis('off')
    fig = plt.gcf()
    fig.set_size_inches(30,30)
    if file is None :
        fig.show()
    else :
        fig.savefig(file)
    return pos

nano_ani = get_ani_data("/data/moritz/0079_pelaginet/fastani_nano.all_v_all.tsv")
poly_ani = get_ani_data("/data/moritz/0079_pelaginet/fastani_poly.all_v_all.tsv")
pelagic_ani = get_ani_data("/data/moritz/0079_pelaginet/fastani.all_v_all.tsv")
pelagic_ani_filter = get_ani_data("/data/moritz/0079_pelaginet/fastani.all_v_all.tsv", checkm="/data/moritz/0079_pelaginet/checkm.txt")


g = (ggplot(data=nano_ani.loc[nano_ani.same_genus], mapping=aes(x='ani'))+geom_histogram(bins = 100) +geom_vline(xintercept = 95, color = "red")+xlim(80, 100)).draw()
g.show()
g.savefig("ani_gap_nano.pdf", format='pdf')
pp = simple_ani_graph(nano_ani, file = "nano_graph_no_col.pdf")
_ = simple_ani_graph(nano_ani, cols = True, pos = pp, file = "nano_graph_col.pdf")

g = (ggplot(data=poly_ani.loc[poly_ani.same_genus], mapping=aes(x='ani'))+geom_histogram(bins = 100) +geom_vline(xintercept = 95, color = "red")+xlim(80, 100)).draw()
g.show()
g.savefig("ani_gap_poly.pdf", format='pdf')
_ = simple_ani_graph(poly_ani, cols = True, file = "poly_graph_col.pdf")

g = (ggplot(data=pelagic_ani.loc[pelagic_ani.same_genus], mapping=aes(x='ani'))+geom_histogram(bins = 100) +geom_vline(xintercept = 95, color = "red")+xlim(80, 100)).draw()
g.show()
g.savefig("ani_gap_pelagic.pdf", format='pdf')
pos = simple_ani_graph(pelagic_ani, cols = True, file = "pelagic_graph_col.pdf")
_ = simple_ani_graph(pelagic_ani_filter, pos = _, cols = True, file = "pelagic_graph_col.pdf")

pelagic_gs = pelagic_ani_filter.genus0.value_counts()[0:9].index
g = (ggplot(data=pelagic_ani_filter.loc[(pelagic_ani_filter.genus0.isin(pelagic_gs) & pelagic_ani_filter.same_genus)], mapping=aes(x='ani'))+geom_histogram(bins = 100) +geom_vline(xintercept = 95, color = "red")+
     xlim(80, 100)+facet_wrap("~genus0", scales = "free_y")+theme(figure_size=(11, 8), dpi=100)
     ).draw()
g.show()
g.savefig("ani_gap_pelagic_by_genus.pdf", format='pdf')

pelagic_gs = pelagic_ani_filter.species0.value_counts()[0:9].index
zz = [ [a for a in ab if a in pelagic_gs][:1]  for ab in zip(pelagic_ani_filter.species0,pelagic_ani_filter.species1)]
pelagic_ani_filter['helper'] = [a[0] if a else "None" for a in zz]
g = (ggplot(data=pelagic_ani_filter.loc[pelagic_ani_filter.helper != "None"], mapping=aes(x='ani', fill = "same_species"))+geom_histogram(bins = 100) +
     geom_vline(xintercept = 95, color = "red")+xlim(80, 100)+facet_wrap("~helper", scales = "free_y")+theme(figure_size=(11, 8), dpi=100)
     ).draw()
g.show()
g.savefig("ani_gap_pelagic_by_species.pdf", format='pdf')



pelaAs = [g for g,v in genuses.items() if v == "g__Pelagibacter_A" ]
pelagAs_ani = pelagic_ani_filter.loc[(pelagic_ani_filter.genus0 == "g__Pelagibacter_A") & (pelagic_ani_filter.same_genus)]
pp = simple_ani_graph(pelagAs_ani, cols=True, cutoff = 85)


anis = parse_fastani_output("/data/moritz/0079_pelaginet/fastani.all_v_all.tsv")
checkm_genomes = checkm_parser("/data/moritz/0079_pelaginet/checkm.txt", 70, 5)
anis.filter_genomes(checkm_genomes)
clustering = cliqueblocksClustering(default_guidetree, anis, size_cutoff=0, gap_size = 1.5, strain_cutoff=99)
clustering.cluster_simple()

genome2clust = { "_".join(g.split("/")[-1].split("_")[:2]) : i+1 for i,gg in enumerate(clustering.final_clusters) for g in gg}
pelagic_ani_filter['clust0'] = pelagic_ani_filter['g0'].map(genome2clust).fillna(0).astype(int)
pelagic_ani_filter['clust1'] = pelagic_ani_filter['g1'].map(genome2clust).fillna(0).astype(int)
pelagic_ani_filter['same_clust'] = pelagic_ani_filter['clust0'] == pelagic_ani_filter['clust1']

pelagic_gs = pelagic_ani_filter.clust0.value_counts()[0:12].index
zz = [ [a for a in ab if a in pelagic_gs][:1]  for ab in zip(pelagic_ani_filter.clust0,pelagic_ani_filter.clust1)]
pelagic_ani_filter['helper'] = [a[0] if a else "None" for a in zz]
g = (ggplot(data=pelagic_ani_filter.loc[pelagic_ani_filter.helper!= "None"], mapping=aes(x='ani', fill = "same_clust"))+geom_histogram(bins = 30) +
     geom_vline(xintercept = 95, color = "red")+facet_wrap("~helper", scales="free_y" )+xlim(85, 100)+theme(figure_size=(14, 8), dpi=100)+theme_bw()
     ).draw()
g.show()
g.set_size_inches(14,8)
g.savefig("ani_gap_pelagic_by_cliquebait.pdf", format='pdf')


# position map
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import networkx
import matplotlib

def make_video():
    plt.clf()
    matplotlib.rcParams['figure.figsize'] = 700,700

    inclique_edges = { frozenset((g1,g2))  for cluster in clustering.final_clusters for g1 in cluster for g2 in cluster if g1 != g2 and len(cluster) > 4}
    inclusters = frozenset.union(*inclique_edges)

    edges = { k  for k, v in clustering.ani_dictionary.items() if v > 92 if any([g in inclusters for g in k])}

    genomes = list(frozenset.union(*edges))
    self_edges = [(g,g) for g in genomes]

    graph = networkx.Graph(list(edges.union(inclique_edges))+ self_edges)


    pos = networkx.spring_layout(graph, method = "energy", gravity = 20, seed = 123)

    data = clustering.get_clusters_stats()
    classes = {gg : dd['id'] for dd in  data for gg in dd['genomes']}
    communities =  [classes.get(g, 0) for g in graph]
    numCommunities = max(communities)+1
    # color map from http://colorbrewer2.org/
    cmapLight  = colors.LinearSegmentedColormap.from_list('custom pastel', ['#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f', '#cab2d6'], N=numCommunities)
    cmapDark  = colors.LinearSegmentedColormap.from_list('custom dark', ['#1f78b4', '#33a02c', '#e31a1c', '#ff7f00', '#6a3d9a'], N=numCommunities)


    from tqdm import tqdm
    for i in tqdm(list(range(100, 89, -1))):
        plt.clf()
        fig = plt.gcf()
        fig.set_size_inches(100,100)
        edges = {k for k,v in clustering.ani_dictionary.items() if v > i and all([g in genomes for g in k])}
        graph = networkx.Graph(list(edges)+ self_edges)
        communities =  [classes.get(g, 0) for g in graph]
        edgelist = [e for e in graph.edges if e not in networkx.selfloop_edges(graph)]
        edge_cols = ['#888888' if frozenset(e) not in inclique_edges else cmapLight(classes[e[0]]) for e in edgelist ]
        edge_alpha = [0.1 if frozenset(e) not in inclique_edges else 1 for e in edgelist ]
        networkx.draw_networkx_edges(graph, pos, edge_color=edge_cols, edgelist=edgelist, width=5, alpha=edge_alpha)
        networkx.draw_networkx_nodes(graph, pos, alpha=0.2, node_color=[cmapDark(classes.get(g, "#000000")) for g in graph], node_size=100)    
        plt.axis('off')
        fig.savefig(f"network_{i:03d}.png")

clustering2= cliqueblocksClustering(default_guidetree, anis, size_cutoff=7)
clustering2.cluster_simple()
clustering2.draw_histograms(file="clustering2_histograms.pdf")