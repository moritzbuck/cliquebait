library(data.table)
library(ggplot2)
library(hues)

#full = fread("/home/moritz/data/0079_pelaginet/cliquebait_nf/cliquebait_full.csv")
#summ = fread("/home/moritz/data/0079_pelaginet/cliquebait_nf/cliquebait_summary.csv")

full = fread("cliquebait_full.csv")
full[, cluster_id := sub("g__", "", cluster_id)]

summ = fread("cliquebait_summary.csv")
summ[, V1 := sub("g__", "", V1)]

kept = summ[stab_clique_ani < 95 & nb_genomes > 40]$V1
motu_kept = summ[stab_clique_ani > 95 & nb_genomes > 40]$V1

ggplot(full[cluster_id %in% kept], aes(x = ani, group=cluster_id, coll=family))+geom_density()+facet_wrap(~phylum)+geom_vline(xintercept = 95, col = "red")+theme(text=element_text(size=30))

ggplot(full[cluster_id %in% kept], 
       aes(x = ani, group = cluster_id, coll = family)) +
  geom_density() +
  geom_vline(xintercept = 95, col = "red")

ggplot(summ, aes(x=mean_ani, y=after_stat(ndensity), col = stab_clique_ani < 95)) +
  geom_freqpoly()

ggplot(summ, aes(x=mean_ani, col = stab_clique_ani < 95)) +
  geom_freqpoly()

qc = function(tax){
gmms = lapply(1:5, function(x) Mclust(full[cluster_id == "g__Campylobacter_B_01"]$ani, G=x))

logliks = sapply(gmms, function(x) x$loglik)
logliks[2:5]/logliks[1:4]
bics = sapply(gmms, function(x) x$bic)

mean_diffs = sapply(gmms, function(x) max(x$parameters$mean)-min(x$parameters$mean))
}

ggplot(full[cluster_id %in% kept],
       aes(x = ani, y=after_stat(ndensity), group=cluster_id, fill = class))+
  geom_histogram(bins =100)+xlim(90,100)+scale_fill_manual(values = iwanthue(9))+
  geom_vline(xintercept = 95, col = "red", size = 5)+guides(fill = FALSE)+
      facet_wrap(~cluster_id, scale = "free", ncol=6)+theme(text = element_text(size = 30), strip.text=element_text(size=50), axis.label.text = element_text(size = 40))+
      xlab("ANI (average nucleotide identity)")+ylab("Density")

  ggsave("weird_clades.pdf", width = 48, height = 28)

  ggplot(full[cluster_id %in% sample(motu_kept, 24)],
       aes(x = ani, y=after_stat(ndensity), group=cluster_id, fill = class))+
  geom_histogram(bins =100)+xlim(90,100)+scale_fill_manual(values = iwanthue(9))+
  geom_vline(xintercept = 95, col = "red", size = 5)+guides(fill = FALSE)+
      facet_wrap(~cluster_id, scale = "free", ncol=6)+theme(text = element_text(size = 30), strip.text=element_text(size=50), axis.label.text = element_text(size = 40))+
      xlab("ANI (average nucleotide identity)")+ylab("Density")

  ggsave("normal_clades.pdf", width = 48, height = 28)

summ[, clade_type := "normal"]
summ[stab_clique_ani < 95 , clade_type := "weird"]
ggplot(summ, aes(x = clade_type, y = core_length/mean_est_genome_size))+geom_violin()+geom_boxplot(width = 0.20)+theme(text = element_text(size = 50))+ylab("core fraction")
ggplot(summ, aes(x = clade_type, y = mean_est_genome_size))+geom_violin()+geom_boxplot(width = 0.20)+theme(text = element_text(size = 50))+ylab("genome size (nb of gene-clusters)")
ggsave("genome_size.pdf", width = 10, height = 14)

ggplot(summ, aes(x = stab_clique_ani, y = core_length/mean_est_genome_size))+geom_point(size = 5)+theme(text = element_text(size = 50))+ylab("Fraction of genome that is core")+xlab("clique ANI")
ggsave("core_fract.pdf", width = 20, height = 14)
