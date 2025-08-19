from statistics import mean
from cliquebait import get_verbose
import os


def parse_fastani_output(fastani_file: str, na_val = None) -> dict:
    """
    Parse the output of fastANI and return a dictionary with relevant information.
    Or any tab-separated text file with three columns where the two first columns are the name of the genomes
    and the third column is the ANI value (between 0 and 100). 
    """

    if not os.path.isfile(fastani_file):
        raise FileNotFoundError(f"File not found: {fastani_file}")

    with open(fastani_file) as handle:
        anis = [l.strip().split() for l in handle]
    anis = {(v[0], v[1]) : float(v[2]) for v in anis}

    dreped_anis = {frozenset((k[0],k[1])) :  v if (k[1], k[0]) not in anis else mean([v, anis[(k[1], k[0])]]) for k, v in anis.items() if k[0] != k[1]}

    return ani_dict(dreped_anis, na_val=na_val)

from collections.abc import MutableMapping

class ani_dict(MutableMapping):
    def __init__(self, dictionary = dict(), na_val=0.0):
        self.store = dict()
        self.genomes = set()
        for k, v in dictionary.items():
            self[k] = v
            for kk in k:
                self.genomes.add(kk)
        self.default = na_val if na_val is not None else min(self.store.values()) 


    def __missing__(self, key):
        return self.default

    def __getitem__(self, key):
        if key not in self.store:
            return self.default
        return self.store[key]

    def __setitem__(self, key, value):
        assert value >= 0.0 and value <= 100.0, "ANI values must be between 0 and 100"
        assert isinstance(key, frozenset) and len(key) == 2, f"Key must be a frozenset of two genome names not {key}"
        self.store[key] = value
        for kk in key:
            self.genomes.add(kk)

    def __delitem__(self, key):
        del self.store[key]

    def __iter__(self):
        return iter(self.store)
    
    def __len__(self):
        return len(self.store)
    
    def getset(self, set_):
        return [self.get(frozenset((g1,g2))) for i,g1 in enumerate(set_) for j, g2 in enumerate(set_) if j > i]
    


