library(data.table)
library(ggplot2)
full = fread("/home/moritz/data/0079_pelaginet/cliquebait_nf/cliquebait_full.csv")
summ = fread("/home/moritz/data/0079_pelaginet/cliquebait_nf/cliquebait_summary.csv")

kept = summ[stab_clique_ani < 95 & nb_genomes > 50]$V1
ggplot(full[cluster_id %in% kept], aes(x = ani, group=cluster_id, coll=family))+geom_density()+facet_wrap(~phylum)


qc = function(tax){
gmms = lapply(1:5, function(x) Mclust(full[cluster_id == "g__Campylobacter_B_01"]$ani, G=x))

logliks = sapply(gmms, function(x) x$loglik)
logliks[2:5]/logliks[1:4]
bics = sapply(gmms, function(x) x$bic)

mean_diffs = sapply(gmms, function(x) max(x$parameters$mean)-min(x$parameters$mean))
}


 ggplot(full[cluster_id %in% kept], aes(x = ani, y=after_stat(ndensity), group=cluster_id, fill = phylum))+
    geom_histogram(bins =100)+xlim(90,100)+geom_vline(xintercept = 95, col = "red")+
    facet_wrap(~cluster_id, scale = "free")