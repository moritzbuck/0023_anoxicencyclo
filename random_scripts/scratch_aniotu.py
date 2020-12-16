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




header = ["query_name" , "ortholog" , "evalue" , "score" , "taxo" , "name" , "GO-terms" , "EC" , "KO" , "Pathway" , "Module" , "Reaction" , "rclass" , "BRITE" , "KEGG_TC" , "CAZy " , "BiGG" , "tax_scope" , "OG" , "deprecated" , "COG_category" , "description"]

all_motus = set(mag_md.mOTU)
all_motus.remove(None)

def find_rep(motu, method = "checkm", tolerance = 5, max_contamin = 5):
        if method == "checkm":
            magdat = mag_md.loc[mag_md.mOTU==motu,["completeness","contamination"]].to_dict(orient = "index")
            magdat = {k : v for k,v in magdat.items() if v['contamination'] < max_contamin}
            most_complete = max(magdat.values(), key = lambda a : a['completeness'])['completeness']
            cands = {k : v for k,v in magdat.items() if v['completeness'] > (most_complete - tolerance)}
            winner = min(cands.items(), key = lambda a : a[1]['contamination'])
            return winner[0]

motu2rep = {m : find_rep(m) for m in all_motus}

rep_dir = "/home/moritz/data/data_submit/representative_genomes"
os.makedirs(rep_dir, exist_ok=True)

for k,v in tqdm(motu2rep.items()):
    patty = mag_md.loc[v]['path']
#    shutil.copy(pjoin(rep_dir, "..",patty ), rep_dir)
    call("tar -C" + rep_dir + " --strip-components 3 -xzf " + pjoin(raw_folder, patty) + " " + pjoin(patty[:-7], os.path.basename(patty[:-7]) + ".fna.gz"), shell = True)
    patty2 = pjoin(rep_dir, os.path.basename(patty)).split(".")[0]
    shutil.move(patty2 + ".fna.gz", pjoin(rep_dir, k + ".fna.gz"))

with open("/home/moritz/data/data_submit/metadata/anoxic_mOTUs.txt") as handle:
    motu_full = json.load(handle)
bin2motu = {vv['name'] : k for k, v in motu_full.items() for vv in v['genomes']}

mag_md.mOTU = [bin2motu.get(i) for i in mag_md.index]
mag2tax = mag_md.taxonomy.to_dict()



from itertools import takewhile

common_prefix = lambda stri : [c[0] for c in takewhile(lambda x: all(x[0] == y for y in x), zip(*stri))]

def consensus_tax(motu):
    taxos = [mag2tax[bin['name']]  for bin in motu_full[motu]['genomes'] if bin['checkm_contamin'] <5 and bin['checkm_complet'] > 40  ]
    hits = [[hh for hh in h.split(";") if len(hh) > 3] for h in taxos if h != 'not_classified' and h != "disagree"]
    if len(hits) > 0:
        if all([t == taxos[0] for t in taxos ]) :
            out = hits[0]
        else :
            out = __consensus_tax(hits)
    else :
        out = ['d__unclassidied_organism']
    prefix = ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    if len(out) < 7:
        suffix = [p + "unclassified_" + out[-1][3:] if p != "s__" else p + motu for p in prefix[len(out):]]
        return ";".join(out + suffix)
    return ";".join(out)

def __consensus_tax(hits) :
    common_core = common_prefix(hits)
    left = [h for h in hits if h != common_core]
    if len(left) == 0:
        return common_core
    if all([h == left[0] for h in  left]):
        return left[0]
    longest = len(max(left, key = len))-1
    return __consensus_tax([l[:longest] for l in left])

motu2tax = {k : consensus_tax(k) for k in motu_full}

def copy_mag(bin_id, taxonomy, motu):

    ipath = pjoin("/home/moritz/data/data_submit/bins",bin_id.split("_")[0], bin_id)
    call("tar xf {}".format(ipath + ".tar.gz"), shell = True)
    ipath2 = ipath.replace("/home/moritz/data/data_submit/","")
    opath = "/home/moritz/data/0039_mOTUlizer_play/root/" + taxonomy.replace(";","/")
    os.makedirs(pjoin(opath, bin_id), exist_ok = True)
    call("unpigz  {}".format(ipath2 + "/" + bin_id + ".emapper.gz"), shell = True)

    shutil.copy(pjoin(ipath2, bin_id + ".faa.gz"), pjoin(opath, bin_id))
    shutil.copy(pjoin(ipath2, bin_id + ".fna.gz"), pjoin(opath, bin_id))
    shutil.copy(pjoin(ipath2, bin_id + ".ffn.gz"), pjoin(opath, bin_id))
    shutil.copy(pjoin(ipath2, bin_id + ".gff.gz"), pjoin(opath, bin_id))
    shutil.copy(pjoin(ipath2, bin_id + ".emapper"), pjoin(opath, bin_id))

    shutil.rmtree("bins")
    while not opath.endswith("_play"):
        if os.path.exists(pjoin(opath, os.path.basename(opath) + ".gids")):
            with open(pjoin(opath, os.path.basename(opath) + ".gids")) as handle:
                lines = handle.readlines()
                if bin_id + "\n" not in lines:
                    lines += [bin_id + "\n"]
            os.remove(pjoin(opath, os.path.basename(opath) + ".gids"))
        else :
            lines = [bin_id + "\n"]
        with open(pjoin(opath, os.path.basename(opath) + ".gids"), "w") as handle:
            handle.writelines(lines)
        to_rm = [f for f in os.listdir(opath) if f.startswith( os.path.basename(opath)) and not f.endswith(".gids")]
        for f in to_rm :
            os.remove(pjoin(opath, f))
        opath = os.path.dirname(opath)

def remove_mag(bin_id, taxonomy):
    opath = "/home/moritz/data/0039_mOTUlizer_play/root/" + taxonomy.replace(";","/")

    if os.path.exists(pjoin(opath,bin_id)):
        shutil.rmtree(pjoin(opath, bin_id))

        while not opath.endswith("_play"):
            with open(pjoin(opath, os.path.basename(opath) + ".gids")) as handle:
                lines = [l for l in  handle.readlines() if l != bin_id + "\n"]
            with open(pjoin(opath, os.path.basename(opath) + ".gids"), "w") as handle:
                handle.writelines(lines)
            to_rm = [f for f in os.listdir(opath) if f.startswith( os.path.basename(opath)) and not f.endswith(".gids")]
            for f in to_rm :
                os.remove(pjoin(opath, f))
            opath = os.path.dirname(opath)

for k, v in tqdm(tt2[["taxonomy", "mOTU"]].iterrows()):
    if v['mOTU'] in to_fix:
        remove_mag(k, broken_motu2tax[v['mOTU']])

found = 0
for k, v in tqdm(tt2[["taxonomy", "mOTU"]].iterrows()):
    if v['mOTU'] in to_fix:
        copy_mag(k, motu2tax[v['mOTU']])
        found += 1
    elif "s__" in v['taxonomy']:
        copy_mag(k, v['taxonomy'])
        found += 1

with open("checkm_file.txt", "w") as handle:
    handle.write("Bin Id\tCompleteness\tContamination\n")
    for k, v in tqdm(tt2[["taxonomy", "mOTU", "completeness", "contamination"]].iterrows()):
        if v['mOTU'] or "s__" in v['taxonomy']  :
            handle.write("{}\t{}\t{}\n".format(k, v['completeness'], v['contamination']))
