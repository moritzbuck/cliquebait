from cliqueblocks import get_verbose
from statistics import mean


def _get_cluster_stats(cluster, anis, genomes, parent_id = None, top_cutoff = 99.9 , quantile_filter = 0.95 ):
    if cluster.is_leaf():
        return { 'clique_ani'      : 100.0,
                 'stab_clique_ani' : 100.0,
                 'nb_genomes'      : 1,
                 'mean_ani'        : 100,
                 'cliqueness'      : 1.0,
                 'linkage_dist'    : 0.0,
                 'cluster_id'      : cluster.id,
                 'parent_id'       : parent_id,
                 'is_leaf'         : cluster.is_leaf(),
                 'genomes'         : [genomes[cluster.id]],
        }
        
    leaves = [genomes[n.id] for n in get_leaves(cluster)]
    zanis = anis.getset(leaves)
    cliqueness = sum([a != anis.default for a in zanis])/len(zanis)
    clique_ani = 100.0 if len(zanis) == 0 else min(zanis)

    sub_anis = sorted([a for a in zanis if a < top_cutoff], reverse=True)    
    if len(sub_anis) == 0:
        stab_clique_ani = clique_ani
    else :
        stab_clique_ani = sub_anis[int(len(sub_anis)*quantile_filter)]

    nb_genomes = len(leaves)
    mean_ani = mean([a for a in zanis if a ])
    return { 'clique_ani'      : clique_ani,
             'stab_clique_ani' : stab_clique_ani,
             'nb_genomes'      : nb_genomes,
             'mean_ani'        : mean_ani,
             'cliqueness'      : cliqueness,
             'linkage_dist'    : cluster.dist,
             'cluster_id'      : cluster.id,
             'parent_id'       : parent_id,
             'is_leaf'         : cluster.is_leaf(),
             'genomes'         : leaves 
        }


def get_leaves(cluster):
    if cluster.is_leaf() :
        return [cluster]
    else :
        return get_leaves(cluster.left) + get_leaves(cluster.right)

def iterate_nodes(cluster, parent = None , cutoff = 0.01):
    if cluster.is_leaf() or cluster.dist < cutoff:
        return [(cluster, None if not parent else parent.id)]
    else :
        return iterate_nodes(cluster.left, cluster) + iterate_nodes(cluster.right, cluster) + [(cluster, None if not parent else parent.id)]

def leaf2cluster(leaf, stats, gap_size, strain_cutoff = 99, bottom_cutoff = 90):
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

def cluster_simple(ani_dictionary, guide_tree, gap_size = 2.5, strain_cutoff = 97.5, size_cutoff = 5, bottom_cutoff = 89, denoising_cutoff = 0.5):
    ori_set = {g for gg in ani_dictionary.keys() for g in gg}
    to_cluster = ori_set.copy()
    final_clusters = []
    while to_cluster:
        tree, genomes = guide_tree.compute_tree(ani_dictionary, subset=to_cluster)

        nstats = {c[0].id :  _get_cluster_stats(c[0], ani_dictionary, genomes, parent_id = c[1]) for c in  iterate_nodes(tree, cutoff = denoising_cutoff)}


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
