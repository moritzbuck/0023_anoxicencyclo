import json
from tqdm import tqdm

def parse_fastqc_file(lib):
    with open('fastqcs/{}/fastqc_data.txt'.format(lib)) as handle:
        lines = "".join(handle.readlines())
        by_modules = [l.split("\n") for l in lines.split(">>END_MODULE\n")]
        by_modules = [[l for l in m if not l.startswith("##") and l != ''] for m in by_modules]
        by_modules = {m[0].split("\t")[0][2:] : [{ k : ll for k,ll in zip(m[1][1:].split("\t"),l.split("\t")) } for l in m[2:] ]  for m in by_modules[:-1]}
        return by_modules

libs = [l for l in os.listdir('fastqcs/') if l.endswith("_fastqc")]
data_by_lib = { l : parse_fastqc_file(l) for l in tqdm(libs)}
with open("fastqc_data_by_library.json", "w") as handle:
    json.dump(data_by_lib, handle, indent=2, sort_keys=True)

samples = {l.split("_")[0] : {'fwd' : [], 'rev' : []} for l in libs}
for l in libs:
    samples[l.split("_")[0]]['fwd' if "_R1_" in l else 'rev'] += [l]

def join_q_hists(sample):
    fwds = samples[sample]['fwd']
    revs = samples[sample]['rev']

    fwds_counts = [ll for l in fwds for ll in data_by_lib[l]['Per sequence quality scores']]
    revs_counts = [ll for l in revs for ll in data_by_lib[l]['Per sequence quality scores']]

    fwd_full = {float(k['Quality']) : { 'count' : 0, 'direction' : "forward-reads" }  for k in fwds_counts}
    for ll in fwds_counts:
        fwd_full[float(ll['Quality'])]['count'] += float(ll['Count'])

    rev_full = {float(k['Quality']) : { 'count' : 0, 'direction' : "reverse-reads" }  for k in revs_counts}
    for ll in revs_counts:
        rev_full[float(ll['Quality'])]['count'] += float(ll['Count'])
    nb_reads = sum([v['count'] for v in rev_full.values()])

    for k,v in fwd_full.items():
        fwd_full[k]['proportion'] = float(v['count']/nb_reads)
    for k,v in rev_full.items():
        rev_full[k]['proportion'] = float(v['count']/nb_reads)

    out_data = ["{},{},{},{},{}\n".format(sample.replace('-','_') if sample.startswith('Loc') else  sample, k, v['count'], v['proportion'] , v['direction']) for  k, v in fwd_full.items()]
    out_data += ["{},{},{},{},{}\n".format(sample.replace('-','_') if sample.startswith('Loc') else  sample, k, v['count'], v['proportion'] , v['direction']) for  k, v in rev_full.items()]

    return out_data

q_hist = ["sample,q_score,count,proportion,direction\n"] + [ll for sample in samples for ll in join_q_hists(sample)]
with open("q_hist.csv", "w") as handle:
    handle.writelines(q_hist)
