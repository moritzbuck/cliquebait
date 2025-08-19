from scipy.cluster.hierarchy import linkage, to_tree
from scipy.spatial.distance import pdist, squareform
from numpy import array

available_guidetrees = [
    "complete_hierarchical",
    "average_hierarchical",
    "single_hierarchical",
    "ward_hierarchical",
    "centroid_hierarchical",
    "hierarchical"
]


def get_guidetree_class(method):
    if method not in available_guidetrees:
        raise ValueError(f"Unknown guide tree method: {method}, available methods are: {available_guidetrees}")
    class_ = method.split("_")[-1]
    arguments_ = method.split("_")[:-1]
    if class_ == "hierarchical":
        return Hierarchicalnodeing(method=arguments_[0] if arguments_ else "complete")
    return Hierarchicalnodeing()

class Hierarchicalnodeing:
    def __init__(self, method="complete"):
        self.method = method
        self.genomes = None
        self.tree = None
        self.links = None
        self.node2parent = {}


    def compute_tree(self, ani_dictionary : dict, subset : set = None):
        # Implement the tree computation logic here
        if not subset:
            self.genomes = ani_dictionary.genomes
        else:
            self.genomes = subset

        self.genomes = list(self.genomes)
        square = [[(100 - ani_dictionary[frozenset((g1,g2)) ])/100 for g2 in self.genomes ] for g1 in self.genomes]
        panis = pdist(array(square, dtype = float))
        self.links = linkage(panis, method=self.method)
        self.tree = to_tree(self.links)
        self.node2parent = { node.id : parent for node, parent in self.iterate_nodes_parents()}

    def iterate_nodes_parents(self):
        return self._iterate_node_parents(self.tree)

    def _iterate_node_parents(self, node, parent = None):
        if node.is_leaf():
            return [(node, None if not parent else parent.id)]
        else :
            return self._iterate_node_parents(node.left, node) + self._iterate_node_parents(node.right, node) + [(node, None if not parent else parent.id)]

    def iterate_nodes(self):
        return self._iterate_nodes(self.tree)

    def get_nodes_info(self):
        return [self.get_node_info(node) for node in self.iterate_nodes()]

    def _iterate_nodes(self, node):
        if node.is_leaf():
            return [node]
        else :
            return self._iterate_nodes(node.left) + self._iterate_nodes(node.right) + [node]

    def get_leaves(self, node = None):
        if node is None:
            node = self.tree
        return self._get_leaves(node)

    def _get_leaves(self, node):
        if node.is_leaf() :
            return [node]
        else :
            return self._get_leaves(node.left) + self._get_leaves(node.right)

    
    def get_node_info(self, node):
        if node.is_leaf():
            return {
                'id': node.id,
                'is_leaf': True,
                'genomes': [self.genomes[node.id] ],
                'parent_id' : self.node2parent.get(node.id, None)
            }
        else:
            return {
                'id': node.id,
                'is_leaf': False,
                'genomes': [self.genomes[leaf.id] for leaf in self.get_leaves(node)],
                'parent_id' : self.node2parent.get(node.id, None)
            }
        
    def find_nodes(self, clstr):
        return {node.id for node in self.iterate_nodes()  if clstr.intersection(self.get_node_info(node)['genomes'])}

    def draw_dendrogram(self, clusters, file = None):
        from matplotlib import pyplot as plt
        from matplotlib import colormaps
        from scipy.cluster.hierarchy import dendrogram

        clusters = sorted(clusters, key = len, reverse= False)
        node2cluster = {i : clstr for clstr in clusters for i in self.find_nodes(clstr)}
        cmapLight = colormaps['gist_rainbow'].resampled(len(clusters))

        colorMap = { cluster : '#%02x%02x%02x' %  tuple([int(f*255) for f in  cmapLight(i)[:3]]) for i, cluster in enumerate(clusters)}
        plt.figure(figsize=(10, 7))
        dendrogram(self.links, color_threshold = None, link_color_func = lambda x : "#000000" if node2cluster.get(x, None) == None else colorMap[node2cluster[x]])
        plt.title('Dendrogram')
        plt.xlabel('Sample Index')
        plt.ylabel('Distance')
        if file:
            plt.savefig(file, format='pdf')
        else:
            plt.show()
        