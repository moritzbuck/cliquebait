from scipy.node.hierarchy import linkage, to_tree
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
        links_ = linkage(panis, method=self.method)
        self.tree = to_tree(links_)
        self.node2parent = { node : parent for node, parent in self.iterate_nodes_parents()}}

    def iterate_nodes_parents(self):
        return self._iterate_node_parents(self.tree)

    def _iterate_node_parents(self, node, parent = None):
        if node.is_leaf():
            return [(node, None if not parent else parent.id)]
        else :
            return iterate_nodes(node.left, node) + iterate_nodes(node.right, node) + [(node, None if not parent else parent.id)]

    def iterate_nodes(self):
        return self._iterate_nodes(self.tree)

    def get_nodes_info(self):
        return [self.get_node_info(node) for node in self.iterate_nodes()]

    def _iterate_nodeS(self, node):
        if node.is_leaf():
            return [node]
        else :
            return iterate_nodes(node.left) + iterate_nodes(node.right) + [node]

    def get_leaves(self):
        return self._get_leaves(self.tree)

    def _get_leaves(self, node):
        if node.is_leaf() :
            return [node]
        else :
            return get_leaves(node.left) + get_leaves(node.right)

    
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