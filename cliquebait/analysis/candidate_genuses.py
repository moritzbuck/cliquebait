import pandas

data = pandas.read_csv("/home/moritz/data/0064_bis/metadata/bac120_taxonomy_r226.tsv", sep = "\t", index_col=0, header=None )
data['genus'] = data[1].apply(lambda x : x.split(";")[-2])

with open("cliquebait/analysis/genuses.txt", "w") as out:
    out.write("\n".join([k for k,v in  data.genus.value_counts().items() if 100 < v  < 2000 ]))