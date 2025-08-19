import argparse
import os
import sys
from sys import stderr
from cliquebait.utils import parse_fastani_output
from cliquebait.guidetrees import available_guidetrees, get_guidetree_class
from cliquebait.clustering.cliquebaitClustering import cliqueblocksClustering
import cliquebait
import json


description_text = "TO DO"


def main(**arg): 
    cliquebait.set_verbose(arg['verbose'])
    output = arg['output'][0] if arg['output'] else sys.stdout
    force = arg['force']

    if output != sys.stdout and os.path.exists(output) and not arg['force']:
        print(f"Output file {output} already exists. Use --force to overwrite.", file=stderr)
        sys.exit(1)

    if cliquebait.get_verbose() > 1:
        print("SuperVerbose mode is enabled", file=stderr)
        print(f"Arguments: {arg}", file=stderr)

    gap_size = arg['gap_size'][0]
    strain_cutoff = arg['strain_cutoff'][0]

    bottom_cutoff = arg['lower_cutoff'][0]
    denoising_cutoff = arg['denoising_cutoff'][0]
    size_cutoff = arg['min_size'][0]

    anis = parse_fastani_output(arg['similarities'][0])
    guide_tree_type = arg['guide_tree'][0]

    clustering = cliqueblocksClustering(guide_tree_type, anis, gap_size=gap_size, strain_cutoff=strain_cutoff, bottom_cutoff=bottom_cutoff, denoising_cutoff=denoising_cutoff, size_cutoff=size_cutoff)  
    clustering.cluster_simple()
#    clustering.draw_network(clusters=clustering.final_clusters)
    dendro = clustering.guide_tree.draw_dendrogram(clusters=clustering.final_clusters, file = "dendrogram.pdf")
    clustering.draw_network(clusters=clustering.final_clusters, file = "karate.pdf")

    with open(output, 'w') as f:
        json.dump(clustering.get_clusters_stats(), f, indent=4)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(prog = "cliquebait", description=description_text, epilog = "Let's do this")
    parser.add_argument('--output', '-o', nargs = 1, help = "send output to this file")
    parser.add_argument('--force', '-f', action='store_true', help = "force execution (overwritting existing files)")
    parser.add_argument('--similarities', '-I', nargs = 1, type=str, help = "file containing similarities between MAGs, only accepts output of fastANI right now", required=True)
    parser.add_argument('--gap_size', '-g', nargs = 1, type=float, default = [cliquebait.default_gap_size] , help = f"size of the gap used to cut the guide tree, default : {cliquebait.default_gap_size}")
    parser.add_argument('--min_size', '-m', nargs = 1, type=int, default = [cliquebait.default_min_size] , help = f"minimum size of returned cluster, default : {cliquebait.default_min_size}")
    parser.add_argument('--guide_tree', '-G', nargs = 1, type=str, default = [cliquebait.default_guidetree] , help = f"type of guide tree used, default : {cliquebait.default_guidetree}", choices=cliquebait.guidetrees.available_guidetrees)
    parser.add_argument('--strain_cutoff', '-s', nargs = 1, type=float, default = [cliquebait.default_strain_cutoff], help = f"similarities above that threshold are ignored in the clustering (considered intra-species clusters, e.g. strain), default {cliquebait.default_strain_cutoff}")
    parser.add_argument('--lower_cutoff', '-l', nargs = 1, type=float, default = [cliquebait.default_lower_cutoff], help = f"similarities above that threshold are ignored in the clustering (considered intra-species clusters, e.g. strain), default {cliquebait.default_lower_cutoff}")
    parser.add_argument('--denoising_cutoff', '-d', nargs = 1, type=float, default = [cliquebait.default_denoising_cutoff], help = f"Similarities of a branch are ordered in decreasing order above that threshold are ignored in the clustering (considered intra-species clusters, e.g. strain), default {cliquebait.default_strain_cutoff}")
    parser.add_argument('--version','-V', action="store_true", help = "get version and exit")
    parser.add_argument('--verbose','-v', nargs="?", const=1, default=1, help = "Set the verbosity level (0-3), default 1. 0 is silent, 1 (or no value) is normal, above is debug mode")

    args = parser.parse_args()
    
    if args.version:
        print(f"{os.path.basename(__file__)} Version {cliquebait.__version__}")
        sys.exit(1)

    main(**vars(args))
