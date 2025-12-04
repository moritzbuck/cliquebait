import json
import os 
from tqdm import tqdm 
import pandas 
cliquebait_out = '/home/moritz/data/0079_pelaginet/cliquebait_nf/cliquebaits/'  

motupan_out = '/home/moritz/data/0079_pelaginet/cliquebait_nf/motupan/'  

gtdb_tax = "/home/moritz/data/dbs/gtdb/release226/taxonomy/bac120_taxonomy_r226_reps.tsv"
gtdb_ar_tax = "/home/moritz/data/dbs/gtdb/release226/taxonomy/ar53_taxonomy_r226_reps.tsv"

tax_strs = pandas.concat([pandas.read_csv(gtdb_tax, sep='\t', header=None, index_col=0)[1], pandas.read_csv(gtdb_ar_tax, sep='\t', header=None, index_col=0)[1]])
levels = ['domain','phylum', 'class', 'order', 'family']
tax_db = pandas.DataFrame.from_dict({tstr.split(";")[-2] : { l : ";".join(tstr.split(";")[:(i+1)]) for i,l in enumerate(levels)} for tstr in tax_strs}, orient='index')

data = []
for fo in tqdm(os.listdir(cliquebait_out)):
        if fo.endswith('.json'):
            with open(os.path.join(cliquebait_out, fo), 'r') as f:
                dd = json.load(f)
                for d in dd:
                    d['genus'] = fo.split('_baits')[0]
                data += dd

def fix_name(gg):
    ss = gg.split("_")
    l = "_".join(ss[:-1])
    r = int(ss[-1])
    return f"{l}_{r:02}"

data_motupan = {}
for fo in tqdm(os.listdir(motupan_out)):
        if fo.endswith('.csv'):
            with open(os.path.join(motupan_out, fo), 'r') as f:
                head = []
                gc_table = []
                for l in f:
                    if l.startswith("#"):
                        head += [l] 
                    else:
                        gc_table += [l]
            hh = gc_table[0].strip().split("\t")
            gcs = [ dict(zip(hh, gc.strip().split("\t")) )for gc in gc_table[1:-1]]
            dd = {l.lstrip("#").split("=")[0] : l.strip().split("=")[1].split(";")[0] for l in head if "posterior" not in l and not l.startswith("#mOTUlizer:mOTUpan") and l !=  '#\n'}
            dd['nb_gcs'] = len(gcs)
            dd['nb_singletons'] = sum([gc['genome_occurences']  == '1' for gc in gcs])
            data_motupan[fix_name(dd['run_name'])] = {k : float(v) for k,v in  dd.items() if k != "run_name"}


df = pandas.DataFrame(data)
df.index = [f"{g}_{id_:02}" for id_ , g in zip(df.id, df.genus)]
#df['genus'] = df.genomes.apply(lambda x: x[0].split("/")[0])
del df['genomes']
df =  df.join(tax_db, on = "genus")
df_full = pandas.DataFrame([dict([("cluster_id", k)] + list(v.items()) + [("ani", vv)] ) for k, v in zip(df.index, data) for vv in v['anis']])
del df_full['genomes']
del df_full['anis']
del df_full['cliqueness']
df_full['genus'] = list(df.loc[df_full.cluster_id].genus)
del df['anis']
df_full =  df_full.join(tax_db, on = "genus")

df = pandas.merge(df, pandas.DataFrame.from_dict(data_motupan, orient = "index"),  left_index=True, right_index=True, how = "outer")
df.to_csv('/home/moritz/data/0079_pelaginet/cliquebait_nf/cliquebait_summary.csv')
df_full.to_csv('/home/moritz/data/0079_pelaginet/cliquebait_nf/cliquebait_full.csv')
