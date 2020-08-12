from Bio import Phylo
import pandas
from tqdm import tqdm
from ete3 import Tree
tree = list(Phylo.parse('data/bac120_r89_unroot.pplacer.tree', 'newick'))[0]

for c in tree.find_clades():
    if not c.is_terminal():
        c.name = None

full_md = pandas.read_csv("data/bac120_metadata_r89.tsv", index_col=0, sep="\t")
mOTUs = {k : [] for k in set(full_md.gtdb_genome_representative)}
for v, k in tqdm(full_md.gtdb_genome_representative.items()):
    mOTUs[k] += [v]

Phylo.write(tree, "data/gtdbtk.tree", "newick")
ete_tree = Tree("data/gtdbtk.tree")

subset = {k : v  for k,v in mOTUs.items() if len(v) > 2}
to_rm = [c for c in list(ete_tree.iter_leaves()) if c.name and c.name in subset]
ete_tree.prune(to_rm)

ete_tree.write(outfile = "data/sub_gtdbtk.tree", "newick")

sub_md = full_md.loc[[l.name for l in tree.find_clades() if l.name]]

anvi_md = pandas.DataFrame.from_dict(dict(sub_md.gtdb_taxonomy.apply(lambda s : {ss.split("__")[0] : ss.split("__")[1] for ss in  s.split(";")})), orient = "index")
anvi_md['nb_genomes'] = [len(mOTUs[k]) for k in anvi_md.index]
anvi_md['has_isolate'] = [any([vv == 'none' for vv in full_md.ncbi_genome_category[mOTUs[k]]]) for k in anvi_md.index]

anvi_md.to_csv("data/full_gtdb_md.tsv", sep="\t", index_label= "ID")

del anvi_md['s']
del anvi_md['g']

anvi_md.loc[[l.name for l in to_rm]].to_csv("data/gtdb_md.tsv", sep="\t", index_label= "ID")


"""
library(ggplot2)
library(data.table)
library(ggrepel)
dd = fread("data/full_gtdb_md.tsv")
dd[,rank := rank(nb_genomes, ties.method="random")]
ggplot(dd[nb_genomes > 10] , aes(y=nb_genomes, x=rank, label=g))+geom_point(size = 5)+scale_y_log10()+theme_minimal()+xlab("rank")+ylab("number of genomes/[MS]AGs per species")+ theme(text = element_text(size=20))+geom_label_repel(data = dd[nb_genomes > 2000], col="black", size=5)
ggsave("zoom_expo.svg")

ggplot(dd , aes(y=nb_genomes, x=rank, label=g))+geom_point(size = 3)+theme_minimal()+xlab("rank")+ylab("number of genomes/[MS]AGs per species")+ theme(text = element_text(size=20))
ggsave("expo.svg")
"""
