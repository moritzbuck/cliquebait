from cliqueblocks import get_verbose
from statistics import mean
from cliqueblocks.guidetrees import get_guidetree_class

class cliqueblocksClustering:

    def __init__(self, guide_tree_class, ani_dictionary, gap_size = 2.5, strain_cutoff = 97.5, size_cutoff = 5, bottom_cutoff = 89, denoising_cutoff = 0.5):
        self.guidetree_class = guide_tree_class
        self.verbose = get_verbose()
        self.ani_dictionary = ani_dictionary
        self.guide_tree = get_guidetree_class(self.guide_tree_class)
        self.guide_tree.compute_tree(self.ani_dictionary)
        self.clusters = None
        self.gap_size = gap_size
        self.strain_cutoff = strain_cutoff
        self.size_cutoff = size_cutoff
        self.bottom_cutoff = bottom_cutoff
        self.denoising_cutoff = denoising_cutoff
        

    def _get_cluster_info(self, cluster):
        if len(cluster) == 1:
            return { 'clique_ani'      : 100.0,
                    'stab_clique_ani' : 100.0,
                    'nb_genomes'      : 1,
                    'mean_ani'        : 100,
                    'cliqueness'      : 1.0
                    }
            
        zanis = self.anis.getset(cluster))
        cliqueness = sum([a != self.anis.default for a in zanis])/len(zanis)
        clique_ani = 100.0 if len(zanis) == 0 else min(zanis)

        sub_anis = sorted([a for a in zanis if a < self.strain_cutoff], reverse=True)    
        if len(sub_anis) == 0:
            stab_clique_ani = clique_ani
        else :
            stab_clique_ani = sub_anis[int(len(sub_anis)*denoising_cutoff)]

        nb_genomes = len(cluster)
        mean_ani = mean([a for a in zanis if a ])
        return { 'clique_ani'      : clique_ani,
                'stab_clique_ani' : stab_clique_ani,
                'nb_genomes'      : nb_genomes,
                'mean_ani'        : mean_ani,
                'cliqueness'      : cliqueness,
            }

    def leaf2cluster(self, leaf, stats, gap_size, strain_cutoff = 99, bottom_cutoff = 90):
        current_stats = stats[leaf]
        parent_id = current_stats['parent_id']
        if parent_id is None :
            return current_stats['genomes']

        parent_stats = stats[parent_id]

        if current_stats['is_leaf'] and parent_stats['stab_clique_ani'] > bottom_cutoff:
            return leaf2cluster(parent_id, stats, gap_size)
        
        if parent_stats['stab_clique_ani'] > strain_cutoff:
            return leaf2cluster(parent_id, stats, gap_size)

        if (current_stats['stab_clique_ani'] - parent_stats['stab_clique_ani']) < gap_size:
            return leaf2cluster(parent_id, stats, gap_size)

        return current_stats['genomes'] 

    def _get_guidetree_stats(self, tree):
        nstats = {node_stats['node_id'] : node_stats for node_stats in tree.get_nodes_info()  }
        for k,v in nstats.items():
            v.update(self.cluster_info(v['genomes']))     
        return nstats
    
    def cluster_simple(self):
        ori_set = self.ani_dictionary.genomes
        to_cluster = ori_set.copy()
        final_clusters = []
        while to_cluster:
            if to_cluster == ori_set:
                tree = self.guide_tree
            else:
                tree = get_guidetree_class(self.guide_tree_class)
                tree.compute_tree(self.ani_dictionary, subset=to_cluster)


            ids_ = [k for k,v in  nstats.items() if v['is_leaf']]
            clusters = []
            for i in ids_ :
                clusters += [frozenset(leaf2cluster(i, nstats,  strain_cutoff = strain_cutoff, gap_size = gap_size, bottom_cutoff = bottom_cutoff))]
            print(f"{len(set(clusters))} unique clusters obtained from {len(ids_)} genomes")
                
            print("derep and merge clusters")

            print(f"{len(set(clusters))} unique clusters left after removing clusters {size_cutoff} or smaller")
            
            clusters = list(set(clusters) - {frozenset(ori_set)})

            to_rm = set()
            for i,c1 in enumerate(clusters):
                for j,c2 in enumerate(clusters):
                    if j > i and len(c2.intersection(c1)) > 0 :
                        if len(c1) > len(c2):
                            to_rm.add(c1)
                        else :
                            to_rm.add(c2)


            final_clusters +=  [c for c in clusters if c not in to_rm and c not in final_clusters]
            
            if to_rm:
                torm = list(frozenset(ori_set) - frozenset().union(*final_clusters))
                print( f"{len(torm)} removed as they were engulfing smalled clusters" )
                to_cluster = torm
            else :
                to_cluster = []

        final_clusters = list(set([c for c in set(final_clusters) if len(c) > size_cutoff]))    
            
        print(f"Final cluster count {len(final_clusters)} accounting for {len(frozenset.union(*final_clusters))} genomes (e.g. {100*len(frozenset.union(*final_clusters))/len(ori_set)}% of the genomes")
        return final_clusters
