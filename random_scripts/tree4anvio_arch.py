import os, shutil
from ete3 import Tree
import pandas
from random import choice
from tqdm import tqdm

tree = Tree('representative_tree/archaea_raw.tree')

gtdb_md = pandas.read_csv("/home/moritz/dbs/gtdb.realease89/ar122_metadata_r89.tsv", sep="\t", index_col=0)
gtdbid2fam = {i :  ";".join(t.split(";")[:5]) for i,t in zip(gtdb_md.index, gtdb_md.gtdb_taxonomy)}
fam2gtdbid = { l : [] for l in set(gtdbid2fam.values())}

for k,v in gtdbid2fam.items():
    fam2gtdbid[v] += [k]

fam2rnd = {k : choice(v) for k,v in fam2gtdbid.items()}
with open("representative_tree/archaea_random_fam_reps.txt", "w") as handle:
    handle.writelines([k + "\t" + v + "\n" for k,v in fam2rnd.items()])

rnds = set(fam2rnd.values())
nodes = list(tree.iter_leaves())

tree.prune([n for n in tree.iter_leaves() if "mOTU" in n.name or n.name in rnds], preserve_branch_length=True)

sub_tree.write(outfile="representative_tree/archaea_pruned.tree")
motu_md = pandas.read_csv("metadata/mOTU_stats.csv", index_col=0)

gtdbid2fam.update(motu_md.consensus_tax.to_dict())

ptree = tree.copy()
for l in ptree.iter_leaves():
    l.name = gtdbid2fam[l.name]

ptree.write(outfile="representative_tree/archaea_pretty.tree")

motu_md = pandas.read_csv("metadata/mOTU_stats.csv", index_col=0)
tt = set(list(tree.iter_leaf_names())).difference(set(motu_md.index))
gtdb_md = pandas.read_csv("/home/moritz/data/gtdb/ar122_metadata_r89.tsv", index_col=0, sep= "\t")
mag_md = pandas.read_csv("metadata/master_table.csv", index_col=0)

tt2 = mag_md.loc[motu_md.representative_MAGs]
tt2.index = tt2.mOTU

motu_md = motu_md.join(tt2)
del motu_md['mean_good_ANIs']
del motu_md['mean_decent_ANIs']
del motu_md['mashmd5']
del motu_md['gtdbtk_classification']
del motu_md['sourmash_classification']
del motu_md['mOTU']

gtdb_md = gtdb_md.loc[tt, ['gtdb_taxonomy', 'genome_size', 'protein_count', 'gc_percentage', 'coding_density', 'checkm_contamination', 'checkm_completeness']]
gtdb_md.columns = ['consensus_tax','length','nb_proteins','GC','coding_density','contamination','completeness']
gtdb_md['genbank_accession'] = [c[3:] for c in gtdb_md.index]
gtdb_md['GC'] = gtdb_md.GC/100
gtdb_md['coding_density'] = gtdb_md.coding_density/100

anvio_table = pandas.concat([motu_md, gtdb_md])

del anvio_table['nb_good_genomes']
del anvio_table['nb_decent_genomes']
del anvio_table['nb_contigs']
del anvio_table['nb_proteins']
del anvio_table['strain_heterogeneity']
del anvio_table['path']
del anvio_table['accession']
del anvio_table['taxonomy']

anvio_table['database'] = ["stratfreshdb"  if "anoxic" in a else "gtdb_r89" for a in anvio_table.index]

anvio_table['novelty'] = ""

anvio_table.loc[["mOTU" in t.split(";")[-1] for t in anvio_table.consensus_tax], "novelty"] = "new species"
anvio_table.loc[["unclassified" in t.split(";")[-2] for t in anvio_table.consensus_tax], "novelty"] = "new genus"
anvio_table.loc[["unclassified" in t.split(";")[-3] for t in anvio_table.consensus_tax], "novelty"] = "new family"
anvio_table.loc[["unclassified" in t.split(";")[-4] for t in anvio_table.consensus_tax], "novelty"] = "new order (or higher)"

tax = pandas.DataFrame.from_dict({k : {kk : vv[3:] for kk,vv in zip(["Phylum", "Class", "Order", "Family"],v.split(";")[1:5]) } for k,v in anvio_table.consensus_tax.items()}, orient="index")
anvio_table = anvio_table.join(tax)
anvio_table= anvio_table.loc[[t.startswith("d__Archaea") for t in anvio_table.consensus_tax]]

anvio_table.to_csv("anvio/ar_metadata.csv", sep="\t", index_label="species")
