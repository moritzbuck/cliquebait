import json
import os 
from tqdm import tqdm 
import pandas 
cliquebait_out = '/home/moritz/data/0079_pelaginet/cliquebait_nf/cliquebaits/'  
gtdb_tax = "/home/moritz/data/dbs/gtdb/release226/taxonomy/bac120_taxonomy_r226_reps.tsv"

tax_strs = pandas.read_csv(gtdb_tax, sep='\t', header=None, index_col=0)[1]
levels = ['phylum', 'class', 'order', 'family']
tax_db = pandas.DataFrame.from_dict({tstr.split(";")[-2] : { l : ";".join(tstr.split(";")[1:(i+2)]) for i,l in enumerate(levels)} for tstr in tax_strs}, orient='index')

data = []
for fo in tqdm(os.listdir(cliquebait_out)):
        if fo.endswith('.json'):
            with open(os.path.join(cliquebait_out, fo), 'r') as f:
                data += json.load(f)

df = pandas.DataFrame(data)
df.index = [f"{g[0].split('/')[0]}_{id_:02}" for id_ , g in zip(df.id, df.genomes)]
df['genus'] = df.genomes.apply(lambda x: x[0].split("/")[0])
del df['genomes']
df =  df.join(tax_db, on = "genus")
df_full = pandas.DataFrame([dict([("cluster_id", k)] + list(v.items()) + [("ani", vv)] ) for k, v in zip(df.index, data) for vv in v['anis']])
del df_full['genomes']
del df_full['anis']
del df_full['cliqueness']
df_full['genus'] = list(df.loc[df_full.cluster_id].genus)
del df['anis']
df_full =  df_full.join(tax_db, on = "genus")
df.to_csv('/home/moritz/data/0079_pelaginet/cliquebait_nf/cliquebait_summary.csv')
df_full.to_csv('/home/moritz/data/0079_pelaginet/cliquebait_nf/cliquebait_full.csv')
