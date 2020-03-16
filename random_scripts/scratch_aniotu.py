simil = 95
pairs = "fastani_pairs.csv"
magstats = "magstats.csv"
otus = "ani_otus.csv"
import igraph
from numpy import  logical_and,logical_or
import pandas

magstats = pandas.read_csv(magstats, index_col=0)
decents = magstats.index[logical_and( magstats.completeness > 40, magstats.contamination <5)]
subs = set(magstats.index ).difference(decents)
with open(pairs) as handle:
    mag_sim = {tuple(l.split()[0:2]) : float(l[:-1].split()[2]) for l in tqdm(handle) if l[0] != "q"  and float(l[:-1].split()[2]) > simil}

def make_ani_otus(sim_dict, simil) :
    valid_pair = lambda k : mag_sim.get((k[1],k[0])) and k[1] in decents and k[0] in decents
    good_pairs = [k for k in tqdm(sim_dict) if valid_pair(k)]
    species_graph = igraph.Graph()
    vertexDeict = { v : i for i,v in enumerate(set([x for k in good_pairs for x in k]))}
    rev_vertexDeict = { v : i for i,v in vertexDeict.items()}
    species_graph.add_vertices(len(vertexDeict))
    species_graph.add_edges([(vertexDeict[k[0]], vertexDeict[k[1]]) for k in good_pairs])
    genome_clusters = [[rev_vertexDeict[cc] for cc in c ] for c in species_graph.components(mode=igraph.STRONG)]
    return genome_clusters

genome_clusters = make_ani_otus(mag_sim, simil)


left_pairs = {k : v for k, v in mag_sim.items() if k[0] not in decents and k[1] in decents and v > 95}
subs = {l[0] : (None,0) for l in left_pairs}
for p,ani in tqdm(left_pairs.items()):
    if subs[p[0]][1] < ani:
        subs[p[0]] = (p[1], ani)

for k, v in subs.items():
    for g in genome_clusters:
        if v[0] in g :
            g += [k]

with open(otus, "w") as handle:
    handle.writelines([ "aniOTU_{id}\t".format(id = i) + ";".join(gs) + "\n" for i, gs in enumerate(genome_clusters)])



plot_ordination(physeq, ord, type='samples', color='when', label= day, title='PCA of the samples from the MiSeq SOP') +theme_minimal()



bacteroidetes <- subset_taxa(physeq, Family  %in% c('Erysipelotrichaceae'))
plot_tree(bacteroidetes, ladderize='left', size='abundance',
          color='when', label.tips='Genus')


@numba.jit(nopython=True, parallel=True)
def get_inters(tt):
     return [(k,l, len(v.intersection(w))/len(v)) for k,v in tt for l,w  in tt if len(v) > 100 and len(w) > 100]

@numba.jit(nopython=True)
def get_simi(l1, l2):
    v = len(l1.intersection(l2))/len(l2)
    return v
