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
        return HierarchicalClustering(method=arguments_[0] if arguments_ else "complete")
    return HierarchicalClustering()

class HierarchicalClustering:
    def __init__(self, method="complete"):
        self.method = method


    def compute_tree(self, ani_dictionary : dict, subset : set = None):
        # Implement the tree computation logic here
        if not subset:
            subset = ani_dictionary.genomes
    
        subset = list(subset)
        square = [[(100 - ani_dictionary[frozenset((g1,g2)) ])/100 for g2 in subset ] for g1 in subset]
        panis = pdist(array(square, dtype = float))
        links_ = linkage(panis, method=self.method)
        tree = to_tree(links_)
        return (tree, subset)
