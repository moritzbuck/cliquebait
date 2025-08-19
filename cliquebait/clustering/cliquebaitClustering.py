from cliquebait import get_verbose
from statistics import mean
from cliquebait.guidetrees import get_guidetree_class

class cliqueblocksClustering:

    def __init__(self, guidetree_class, ani_dictionary, gap_size = 2.5, strain_cutoff = 97.5, size_cutoff = 5, bottom_cutoff = 89, denoising_cutoff = 0.5):
        self.guidetree_class = guidetree_class
        self.verbose = get_verbose()
        self.ani_dictionary = ani_dictionary
        self.guide_tree = get_guidetree_class(self.guidetree_class)
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
            
        zanis = self.ani_dictionary.getset(cluster)
        cliqueness = sum([a != self.ani_dictionary.default for a in zanis])/len(zanis)
        clique_ani = 100.0 if len(zanis) == 0 else min(zanis)

        sub_anis = sorted([a for a in zanis if a < self.strain_cutoff], reverse=True)
        if len(sub_anis) == 0:
            stab_clique_ani = clique_ani
        else :
            stab_clique_ani = sub_anis[int(len(sub_anis)*self.denoising_cutoff)]

        nb_genomes = len(cluster)
        mean_ani = mean([a for a in zanis if a ])
        return { 'clique_ani'      : clique_ani,
                'stab_clique_ani' : stab_clique_ani,
                'nb_genomes'      : nb_genomes,
                'mean_ani'        : mean_ani,
                'cliqueness'      : cliqueness,
            }

    def _leaf2cluster(self, leaf, stats):
        current_stats = stats[leaf]
        parent_id = current_stats['parent_id']
        if parent_id is None :
            return current_stats['genomes']

        parent_stats = stats[parent_id]

        if current_stats['is_leaf'] and parent_stats['stab_clique_ani'] > self.bottom_cutoff:
            return self._leaf2cluster(parent_id, stats)
        
        if parent_stats['stab_clique_ani'] > self.strain_cutoff:
            return self._leaf2cluster(parent_id, stats)

        if (current_stats['stab_clique_ani'] - parent_stats['stab_clique_ani']) < self.gap_size:
            return self._leaf2cluster(parent_id, stats)

        return current_stats['genomes'] 

    def _get_guidetree_stats(self, tree):
        nstats = {node_stats['id'] : node_stats for node_stats in tree.get_nodes_info()  }
        for k,v in nstats.items():
            v.update(self._get_cluster_info(v['genomes']))     
        return nstats
    
    def cluster_simple(self):
        ori_set = self.ani_dictionary.genomes
        to_cluster = ori_set.copy()
        final_clusters = []
        while to_cluster:
            if to_cluster == ori_set:
                tree = self.guide_tree
            else:
                tree = get_guidetree_class(self.guidetree_class)
                tree.compute_tree(self.ani_dictionary, subset=to_cluster)
            
            nstats = self._get_guidetree_stats(tree)

            ids_ = [k for k,v in  nstats.items() if v['is_leaf']]
            clusters = []
            for i in ids_ :
                clusters += [frozenset(self._leaf2cluster(i, nstats))]
            print(f"{len(set(clusters))} unique clusters obtained from {len(ids_)} genomes")
                
            print("derep and merge clusters")

            print(f"{len(set(clusters))} unique clusters left after removing clusters {self.size_cutoff} or smaller")
            
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

        final_clusters = list(set([c for c in set(final_clusters) if len(c) > self.size_cutoff]))    
            
        print(f"Final cluster count {len(final_clusters)} accounting for {len(frozenset.union(*final_clusters))} genomes (e.g. {100*len(frozenset.union(*final_clusters))/len(ori_set)}% of the genomes")
        self.final_clusters = final_clusters


    def draw_network(self, clusters = None, subset = None, gravity = 3, cutoff = 0.95, layout = None, file = None):
        # position map
        import matplotlib.pyplot as plt
        from matplotlib import colors, cm
        import networkx

        plt.clf()

        if subset:
            genomes = subset
        else :
            genomes = self.ani_dictionary.genomes

        edges = { frozenset(k)  for k,v in self.ani_dictionary.items() if all([kk in genomes for kk in k]) and v > cutoff}
        graph = networkx.Graph(list(edges))

        if layout : 
            pos = layout
        else:
            pos = networkx.spring_layout(graph, method = "energy", gravity = gravity)
        if clusters : 
            genome2cluster = {gg : i+1 for i,g in  enumerate(clusters) for gg in g}
            clustering = {g : genome2cluster.get(g, 0) for g in genomes}
        else :
            clustering = {gg : 0 for gg in genomes}

        # community ids
        
        communities =  [clustering.get(g, 0) for g in graph]
        numCommunities = max(communities)+1
        # color map from http://colorbrewer2.org/
        cmapLight = colors.ListedColormap(['#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f', '#cab2d6'], 'indexed', numCommunities)
        cmapDark = colors.ListedColormap(['#1f78b4', '#33a02c', '#e31a1c', '#ff7f00', '#6a3d9a'], 'indexed', numCommunities)

        # edges
        networkx.draw_networkx_edges(graph, pos)

        # nodes
        nodeCollection = networkx.draw_networkx_nodes(graph,
                pos = pos,
                node_color = communities,
                cmap = cmapLight
        )
        # set node border color to the darker shade
        darkColors = [cmapDark(v) for v in communities]
        nodeCollection.set_edgecolor(darkColors)

        # Print node labels separately instead
        for n in graph.nodes():
                plt.annotate(n,
                        xy = pos[n],
                        textcoords = 'offset points',
                        horizontalalignment = 'center',
                        verticalalignment = 'center',
                        xytext = [0, 2],
                        color = cmapDark(clustering.get(n,0))
                )

        plt.axis('off')
        if file:
            plt.savefig(file, format='pdf')
        else:   
            plt.show()

    def get_clusters_stats(self):

        cinfor = [ self._get_cluster_info(cluster) for cluster in self.final_clusters]
        for i, c in enumerate(cinfor):
            c['id'] = i+1
            c['genomes'] = list(self.final_clusters[i])
            c['anis'] = self.ani_dictionary.getset(c['genomes'])
        
        return cinfor