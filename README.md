# DISCLAIMER

THIS IS STILL VERY ALPHA but looks like it works

# CliqueBait

CliqueBait identifies bacterial species, pseudo-species, or metagenomic Operational Taxonomic Units (mOTUs). Not sure yet "which"...
It uses pairwise ANIs from whatever set of genomes you wanna use (as computed by `fastANI`). These are used to build as guide-tree (by default hierarchical cluster with complete linkage).
This tree  will be parsed from the leaves up, and cut whenever the miniumn ANI of the underlying clique drops by more than a predetermined value (2.5%). the underlying "clique-ANI" gets computed by taking all genomes of that branch of the guide-tree and collecting the ANIs of the clique formed by theses (e.g. a network were all are connected to all), if no value is present for a link it is set to the minimum of all ANIs (for now, other options might be implemented in the future if necessary). All ANIs above a "strain-cutoff" will be ignored (maybe clone-cutoff is a better name, default is 99.9) to remove the effect of near clonal genomes from the data (mainly for a future implementation of model fitting I want to do). Also the bottom-tail end of the distribution will be removed to denoise a bit the effect of incomplete genomes (still better to remove them before hand though I guess), default is bottom 5% of ANIs are removed. Whatever is the minimum of the leftover ANIs is the "clique-ANI", basically the defacto cutoff for that clique.

The parsing is repeated for every genome to obtain a node/cluster for each leaf. Some leaves will be in multiple-clusters, so the process will be repeated removing the non-ambiguous clusterings, until all resolved.


### Command-line Usage

Basic usage:

```bash
python clique_bait.py -I fastANI.out -o clustering.json
```

For a full list of options, run:

```bash
python clique_bait.py --help
```

## Output
right now is mainly a `json`-file, with a dictionary for each cluster containing:

* clique_ani : minimum of all ANIs of the clique formed by that cluster
* stab_clique_ani : minimum after denoising (removing "clonal"-ANIs and bottom 5%)
* mean_ani : mean of all ANIs of that cluster
* cliqueness : fraction of all edges of the clique that are defined in the original dataset
* id : id of the cluster
* genomes : list of genomes in that cluster
* anis : all anis of that cluster

Also it generated an ugly netowkr figure (to improve a lot), and a somwhat useable dendrogram which needs to be made more readable too.

## dependencies:

* scipy
* numpy
* matplotlib, networkx, pandas and plotnine if you want the plots
* and probably more to come .....

## License

GNU AFFERO GENERAL PUBLIC LICENSE v3