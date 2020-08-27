import sys
from subprocess import Popen, PIPE, call
import os
from os.path import join as pjoin
from tqdm import tqdm
import xml.etree.ElementTree as ET
import re
import pandas
from numpy import logical_and

user = "Webin-41995"
pwd = "775632"

template_folder = "data"
submission_file = pjoin(template_folder, 'submission.xml')

raw_folder = "/home/moritz/data/data_submit/"
raw_reads = pjoin(raw_folder, "reads/mg/")
sample_sheet = pjoin(raw_folder, "metadata/samples_contextual_working_copy.csv")
mag_sheet = pjoin(raw_folder, "metadata/master_table.csv")
url = "https://www.ebi.ac.uk/ena/submit/drop-box/submit/"
curl_line = 'curl -u {user}:{pwd} -F "SUBMISSION=@{submission_file}" -F "{submission_type}=@{data_file}" "{url}" > {out_file}'

# PARSING SAMPLE SHEET

with open(pjoin(raw_folder, "metadata", "older_libraries.csv")) as handle:
    older_libraries = {l.strip().split(",")[0] : l.strip().split(",")[2] for l in handle.readlines()}

with open(sample_sheet) as handle:
    trash = handle.readline()
    colnames = handle.readline().rstrip().split(",")[1:]
    units = {k : v for k,v in  zip(colnames, handle.readline().rstrip().split(",")[1:])}
    sample_data = {v.split(",")[0] : {c : None if vv in ['NA', '', '-'] else vv for c,vv in zip(colnames,v.rstrip().split(",")[1:])} for v in handle.readlines()}

sample_data = {k : {kk : vv for kk, vv in v.items() if vv} for k, v in sample_data.items()}

units['geographic location (longitude)'] = "DD"
units['geographic location (latitude)'] = "DD"


def fix_coordinates(v):
    allowed = set([ str(i) for i in range(10)] + [" ", "N", "E", "W", "S", ";","."])
    tt = "".join([vv if vv in allowed else " " for vv in v ])
    if "N;" not in tt:
        tt = tt.replace("N", "N;")
    coord = tt.split(";")
    coord = [v.strip() for v in coord]
    coord = [re.split("\ +", v) for v in coord]

    to_num = lambda v : ([float(vv) for vv in v[:-1]], v[-1])
    fill = lambda v : v if len(v) == 3 else v + [0]*(3-len(v))
    to_dec = lambda v : v[0] + v[1]/60 + v[2]/3600
    coord = [to_num(v) for v in coord]
    coord = [(fill(n),d) for n,d in coord]
    coord = [(to_dec(n), d) for n,d in coord]
    coord = [n * (1 if d in ["N","E"] else -1) for n,d in coord]
    return coord
convert_date = lambda d : "-".join(d.split("/")[::-1])

common_md = { 'project name' : 'microbes in stratified freshwaters',
              'ENA-CHECKLIST' : 'ERC000024',
              'sequencing method' : "NovaSeq",
              'investigation type' : 'metagenome',
              'collection data' : 'metagenome',
              'water environmental package' : 'water',
}

for k,v in  sample_data.items():
    coords = fix_coordinates(v['coordinates'])
    v['geographic location (longitude)'] = coords[1]
    v['geographic location (latitude)'] = coords[0]
    v['collection date'] = convert_date(v['collection date'])
    v.update(common_md)
    del v['coordinates']




###############
#MAKE PROJECT #
###############

projet_subm = curl_line.format(user = user, pwd = pwd,
                               submission_file = submission_file,
                               submission_type = "PROJECT",
                               url = url,
                               data_file = pjoin(template_folder, "project.xml"),
                               out_file = pjoin(template_folder, "project.out.xml")
                               )

#call(projet_subm, shell = True)

with open(pjoin(template_folder,"project.out.xml")) as handle:
    root = ET.parse(handle).getroot()

data = dict([(child.tag, child.attrib) for child in root])
project_id = data['PROJECT']['accession']
project_sub_id = data['SUBMISSION']['accession']

################
# MAKE SAMPLES #
################

def make_attributes(key, val, unit = None):
    proto = """    <SAMPLE_ATTRIBUTE>
        <TAG>{key}</TAG>
        <VALUE>{val}</VALUE>{unit_stuff}
    </SAMPLE_ATTRIBUTE>"""

    unit_stuff = "\n        <UNITS>{unit}</UNITS>".format(unit = unit)
    return proto.format(key = key, val = val, unit_stuff = unit_stuff if unit else "")

def make_sample(sample_key, att_dict):
    blacklist = []

    header = """<SAMPLE alias="{name}" center_name="">
    <TITLE>{title}</TITLE>
    <SAMPLE_NAME>
      <TAXON_ID>449393</TAXON_ID>
      <SCIENTIFIC_NAME>freshwater metagenome</SCIENTIFIC_NAME>
      <COMMON_NAME></COMMON_NAME>
    </SAMPLE_NAME>
    <SAMPLE_ATTRIBUTES>
    """
    #1851193	pond metagenome
    title = "{type} {name} sample {sample}".format(type = att_dict['environment (feature)'], name = att_dict['site name'], sample = sample_key)
    formated_header = header.format(name = sample_key, title = title)
    attributes = "".join([make_attributes(k,v, units.get(k, None)) for k,v in att_dict.items() if k not in blacklist])
    return formated_header + attributes + "</SAMPLE_ATTRIBUTES></SAMPLE>\n"


sample_set_xml = """<?xml version="1.0" encoding="UTF-8"?>
<SAMPLE_SET>
{sample}
</SAMPLE_SET>
""".format(sample = "".join([make_sample(k,v) for k,v in sample_data.items()]))

with open(pjoin(template_folder, "all_samples.xml"), "w") as handle:
    handle.writelines(sample_set_xml)

# REGISTERING SAMPLES

all_samples_subm = curl_line.format(user = user, pwd = pwd,
                               submission_file = submission_file,
                               submission_type = "SAMPLE",
                               url = url,
                               data_file = pjoin(template_folder, "all_samples.xml"),
                               out_file = pjoin(template_folder, "all_samples.out.xml")
                   )

#fix_subm = curl_line.format(user = user, pwd = pwd,
#                               submission_file = pjoin(template_folder, 'fix.xml'),
#                               submission_type = "SAMPLE",
#                               url = url,
#                               data_file = pjoin(template_folder, "fix_samples.xml"),
#                               out_file = pjoin(template_folder, "fix_samples.out.xml")
#                   )

#call(all_samples_subm, shell = True)

with open(pjoin(template_folder,"all_samples.out.xml")) as handle:
    root = ET.parse(handle).getroot()

data = [child.attrib for child in root if child.tag == "SAMPLE"]
for dd in data:
    sample_data[dd['alias']]['sample_accession'] = dd['accession']

read_manifest = """STUDY {study_id}
SAMPLE {sample_id}
NAME  microbes in stratified freshwaters - {library_name}
INSTRUMENT  Illumina NovaSeq 6000
LIBRARY_NAME  {library_name}
LIBRARY_SOURCE  METAGENOMIC
LIBRARY_SELECTION RANDOM
LIBRARY_STRATEGY  WGS
FASTQ {fwd}
FASTQ {rev}
"""

title = lambda sample_key, att_dict : "{type} {name} sample {sample}".format(type = att_dict['environment (feature)'], name = att_dict['site name'], sample = sample_key)

for k,v in tqdm(sample_data.items()):
    if not 'exp_id' in v:
        print("Doing", k)
        instanced_manifest = read_manifest.format(study_id = project_id, sample_id = v['sample_accession'], library_name = title(k,v), fwd = pjoin(raw_reads, k + "_R1.fastq.gz"), rev = pjoin(raw_reads, k + "_R2.fastq.gz"))
        with open("data/manifests/{}_manifest.txt".format(k), "w") as handle:
            handle.writelines(instanced_manifest)

samples_to_submit = [pjoin("data/manifests/" + f) for f in os.listdir("data/manifests/") if f.endswith("_manifest.txt")]
parallel_script = "parallel -j{threads} ./helper.sh {{}} ::: {samples}".format(threads = 16, samples = " ".join(samples_to_submit))
call(parallel_script, shell=True)

######################
# PREPING ASSEMBLIES #
######################

raw_assemblies_dir = pjoin(raw_folder, "assemblies")
coassemblies_dir = pjoin(raw_folder, "metadata/coassemblies")
coass_list = [c[:-4] for c in os.listdir(coassemblies_dir)]

single_sample_assemblies = [f[:-8] for f in os.listdir(raw_assemblies_dir) if f[:-8] not in coass_list and (f[:-8] in sample_data or "Loc" in f)  ]

###################################
# SUBMIT SINGLE SAMPLE ASSEMBLIES #
###################################

ass_manifest = """STUDY {study_id}
SAMPLE {sample_id}
ASSEMBLYNAME assembled MG of {sample_name}
ASSEMBLY_TYPE primary metagenome
COVERAGE 1
PROGRAM megahit
PLATFORM  Illumina NovaSeq 6000
MOLECULETYPE genomic DNA
FASTA {seq}
"""

def fix_loc_samples(s):
    if not s.startswith("Loc0") :
        return s
    else :
        return k.replace("_","-")


for k,v in tqdm(sample_data.items()):
    if not 'exp_id' in v and os.path.exists(pjoin(raw_assemblies_dir, fix_loc_samples(k) + ".fna.bz2")):
        print("Doing", fix_loc_samples(k))
        instanced_manifest = ass_manifest.format(study_id = project_id, sample_id = v['sample_accession'], sample_name = k, seq = pjoin(raw_assemblies_dir, fix_loc_samples(k) + ".fna.bz2"))
        with open("data/manifests/{}_assmanifest.txt".format(k), "w") as handle:
            handle.writelines(instanced_manifest)

asses_to_submit = [pjoin("data/manifests/" + f) for f in os.listdir("data/manifests/") if f.endswith("_assmanifest.txt")]

i = 0
while len(asses_to_submit) >0  and i < 20:
    assparallel_script = "parallel -j{threads} ./asshelper.sh {{}} ::: {samples}".format(threads = 20, samples = " ".join(asses_to_submit))
    call(assparallel_script, shell=True)
    asses_to_submit = [pjoin("data/manifests/" + f) for f in os.listdir("data/manifests/") if f.endswith("_assmanifest.txt")]
    print("retry",i)
    i += 1

#########################
# PREPING CO-ASSEMBLIES #
#########################

coassemblies = [f[:-8] for f in os.listdir(raw_assemblies_dir) if f[:-8] in coass_list  if not f.startswith("SA")]

#######################################
# REGISTERING SAMPLES 4 CO-ASSEMBLIES #
#######################################

keys = {dd for d in sample_data.values() for dd in d}

coass2sample = {c : [] for c in coassemblies}
for c in coass2sample:
    with open(pjoin(coassemblies_dir, c + ".txt")) as handle:
        coass2sample[c] = [s.strip() for s in handle]
from datetime import datetime
coass2sample['Loclat'] = [l.replace("-","_") for l in coass2sample['Loclat']]

coass2md = {k : {} for k in coass2sample }
for k,v in coass2sample.items():
    for key in keys:
        values = list(set([sample_data[vv][key] for vv in v if vv[1] != "-" if key in sample_data[vv]]))
        if len(values) == 1 and values[0]:
            coass2md[k][key] = values[0]
        if key == "collection date" and len(values) > 1:
            dates = [datetime.strptime(d, "%Y-%M-%d") for d in values]
            dates = sorted(dates)
            coass2md[k][key] = str(dates[0]).split()[0] + "/" + str(dates[-1]).split()[0]
        if key in ["geographic location (latitude)", "geographic location (longitude)"] and len(values) > 1:
            coass2md[k][key] = 'not provided'

    coass2md[k].update(common_md)
    coass2md[k]['sample derived from'] = ",".join([sample_data[vv]['sample_accession'] if (vv in sample_data) else older_libraries[vv]for vv in v])
    del coass2md[k]['ENA-CHECKLIST']

def make_coassample(sample_key, att_dict):
    blacklist = []

    header = """<SAMPLE alias="{name}" center_name="">
    <TITLE>{title}</TITLE>
    <SAMPLE_NAME>
      <TAXON_ID>449393</TAXON_ID>
      <SCIENTIFIC_NAME>freshwater metagenome</SCIENTIFIC_NAME>
      <COMMON_NAME></COMMON_NAME>
    </SAMPLE_NAME>
    <SAMPLE_ATTRIBUTES>
    """
    #1851193	pond metagenome
    title = "mock sample for coassembly {name} sample of samples: {sample}".format(name = sample_key, sample = att_dict['sample derived from'])
    formated_header = header.format(name = sample_key, title = title)
    attributes = "".join([make_attributes(k,v, units.get(k, None)) for k,v in att_dict.items() if k not in blacklist])
    return formated_header + attributes + "</SAMPLE_ATTRIBUTES></SAMPLE>\n"

coasssample_set_xml = """<?xml version="1.0" encoding="UTF-8"?>
<SAMPLE_SET>
{sample}
</SAMPLE_SET>
""".format(sample = "".join([make_coassample(k,v) for k,v in coass2md.items()]))

with open(pjoin(template_folder, "all_coass_samples.xml"), "w") as handle:
    handle.writelines(coasssample_set_xml)

# REGISTERING SAMPLES

all_samples_subm = curl_line.format(user = user, pwd = pwd,
                               submission_file = submission_file,
                               submission_type = "SAMPLE",
                               url = url,
                               data_file = pjoin(template_folder, "all_coass_samples.xml"),
                               out_file = pjoin(template_folder, "all_coass_samples.out.xml")
                               )

call(all_samples_subm, shell = True)

with open(pjoin(template_folder,"all_coass_samples.out.xml")) as handle:
    root = ET.parse(handle).getroot()

coass_accessions = {}
data = [child.attrib for child in root if child.tag == "SAMPLE"]
for dd in data:
    coass_accessions[dd['alias']] = dd['accession']


#########################
# SUBMIT CO-ASSEMBLIES #
########################
ass_manifest = """STUDY {study_id}
SAMPLE {sample_id}
ASSEMBLYNAME coassembled MG --- {sample_name}
ASSEMBLY_TYPE primary metagenome
COVERAGE 1
PROGRAM megahit
PLATFORM  Illumina NovaSeq 6000
MOLECULETYPE genomic DNA
FASTA {seq}
"""

for k,v in tqdm(coass2md.items()):
        instanced_manifest = ass_manifest.format(study_id = project_id, sample_id = coass_accessions[k], sample_name = k, seq = pjoin(raw_assemblies_dir, k + ".fna.bz2"))
        with open("data/manifests/{}_assmanifest.txt".format(k), "w") as handle:
            handle.writelines(instanced_manifest)

asses_to_submit = [pjoin("data/manifests/" + f) for f in os.listdir("data/manifests/") if f.endswith("_assmanifest.txt")]

i = 0
while len(asses_to_submit) >0  and i < 20:
    i += 1
    print("try",i, ",", len(asses_to_submit), "coassemblies left to submit")
    assparallel_script = "parallel -j{threads} ./asshelper.sh {{}} ::: {samples}".format(threads = 20, samples = " ".join(asses_to_submit))
    call(assparallel_script, shell=True)
    asses_to_submit = [pjoin("data/manifests/" + f) for f in os.listdir("data/manifests/") if f.endswith("_assmanifest.txt")]

################
# PREPING MAGS #
################

mag_md = pandas.read_csv(mag_sheet, index_col=0)

general_md =   { 'sequencing method' : 'Illumina NovaSeq',
                 'assembly software' : 'megahit',
                 'completeness software' : 'checkm',
                 'binning software' : 'metabat',
                 'assembly quality' : 'Many fragments with little to no review of assembly other than reporting of standard assembly statistics',
                 'investigation type' : 'metagenome-assembled genome',
                 'binning parameters' : 'min_len:1500,min_bin_size:10000,maxP:93,minS:50',
                 'taxonomic identity marker' : 'gtdbtk followed by sourmash-lca, see publication',
                 'isolation_source' : 'freshwater',
                 'metagenomic source' : '449393'
                 }
units['completeness score'] = "%"
units['contamination score'] = "%"

mag_md = mag_md.loc[logical_and(mag_md.completeness > 40 , mag_md.contamination < 5) ]
mag_md = mag_md.loc[[not i.startswith("SA") for i in mag_md.index]]

gtdb_md = pandas.concat([pandas.read_csv("/home/moritz/data/gtdb/ar122_metadata_r89.tsv", sep="\t", index_col=0), pandas.read_csv("/home/moritz/data/gtdb/bac120_metadata_r89.tsv", sep="\t", index_col=0)])
gtdb_md = gtdb_md.loc[gtdb_md.gtdb_representative == 't', ['gtdb_taxonomy','ncbi_taxonomy']]
gtdb2ncbi = dict({tuple(v) for k,v in gtdb_md.iterrows()})

from itertools import takewhile

common_prefix = lambda stri : [c[0] for c in takewhile(lambda x: all(x[0] == y for y in x), zip(*stri))]

found_taxo = {}
found_taxo['d__Bacteria;p__Cyanobacteria;c__Cyanobacteriia;o__Cyanobacteriales;f__Nostocaceae;g__Nodularia'] = ('Nodularia', 159191)
found_taxo['not_classified'] = ('root', 1)
found_taxo['d__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Mycobacteriales;f__Mycobacteriaceae;g__Gordonia'] = ('Gordonia', 2053)
found_taxo['d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Chromatiales;f__Chromatiaceae;g__Lamprocystis'] = ('Lamprocystis', 53452)
found_taxo['d__Bacteria;p__Planctomycetota;c__Planctomycetes;o__Planctomycetales;f__Planctomycetaceae;g__Fuerstia'] = ('Fuerstia', 1936111)

def convert2ncbi(tax) :
    tax = ";".join(tax.split(";")[0:6])
    if tax in found_taxo:
        return found_taxo[tax]
    hits = [v for k, v in tqdm(gtdb2ncbi.items()) if k.startswith(tax)]
    hits = [[hh for hh in h.split(";") if len(hh) > 3] for h in hits if h != 'none']
    common_core = common_prefix(hits)[0:6]
    if len(common_core) == 0:
        return convert2ncbi(";".join(tax.split(";")[:-1]))

    trans = {}
    i = 1
    while len(trans) == 0:
        level, taxo = common_core[-i].split("__")
        trans = ncbi.get_name_translator([taxo])
        i += 1

    taxid = trans[taxo]
    if len(taxid) > 1:
        taxid = [tid for tid,rank in ncbi.get_rank(taxid).items() if rank[0] == level]
    assert len(taxid) == 1
    found_taxo[tax] = (taxo, taxid[0])
    return taxo, taxid[0]

mag2md = {}
for k,v in tqdm(mag_md.iterrows()):
    #print("doing",k)
    mag2md[k] = {}
    mag2md[k].update(general_md)
    mag2md[k]['completeness score'] = v['completeness'] if v['completeness'] < 100 else 99.9
    mag2md[k]['contamination score'] = v['contamination']
    mag2md[k]['gtdb_taxonomy'] = v['taxonomy']
    mag2md[k]['scientific_name'], mag2md[k]['tax_id'] = convert2ncbi(v['taxonomy'])
    mag2md[k]['metagenomic OTU'] = v['mOTU']
    ass = k.split("_")[0]
    if ass.startswith("Loc"):
        ass = ass.replace("-", "_")
    if ass[1] == "-":
        ass = ass.replace("-","")
    if ass in sample_data:
        ass_md = sample_data[ass]
    else :
        ass_md = coass2md[ass]
        ass_md['sample_accession'] = coass_accessions[ass]
    mag2md[k].update(ass_md)
    mag2md[k]['sample derived from'] = ass_md['sample_accession']
    mag2md[k]['ENA-CHECKLIST'] = 'ERC000047'
    if 'sample_accession' in mag2md[k]:
        del mag2md[k]['sample_accession']
    if 'Run' in mag2md[k]:
        del mag2md[k]['Run']
    if 'Lake_code' in mag2md[k]:
            del mag2md[k]['Lake_code']

get_parent = lambda taxid : max({k : v for k, v in ncbi.get_lineage_translator(ncbi.get_lineage(taxid)).items() if k != taxid}.items(), key=lambda l : len(l[1]) )[0]
tax2uncul = {}
for k, v in set(found_taxo.values()):
    options = ncbi.get_name_translator(["uncultured " + k + " bacterium", "uncultured " + k + " archaeon"])
    if k == 'root':
        tax2uncul[k] = ('uncultured prokaryote', 198431)
        continue
    if not options:
        options = ncbi.get_name_translator(["uncultured " + k + " sp."])
    if not options:
        options = ncbi.get_name_translator(["uncultured " + k + " cyanobacterium"])
    if not options:
        options = ncbi.get_name_translator([k + " bacterium", k + " archaeon"])
    if not options:
            mama = ncbi.get_taxid_translator([get_parent(v)])
            mama = list(mama.values())[0]
            options = ncbi.get_name_translator(["uncultured " + mama + " bacterium", "uncultured " + mama + " archaeon", "uncultured " + mama + " cyanobacterium"])
    if not options:
            options = ncbi.get_name_translator([k + " sp."])
    if not options:
        tax2uncul[k] = (k,v)
        continue
    assert len(options) == 1
    taxo = list(options.keys())[0]
    tax_id = list(options.values())[0]
    assert len(tax_id) == 1
    tax2uncul[k] = (taxo, tax_id[0])
tax2uncul['Alsobacter'] = ('uncultured Rhizobiales bacterium', 208549)
tax2uncul['Candidatus Paracaedibacteraceae'] = ('uncultured Alphaproteobacteria bacterium', 91750)
tax2uncul['Mageeibacillus'] = ('uncultured Clostridiales bacterium', 172733)
tax2uncul['Candidatus Riflemargulisbacteria'] = ('uncultured Candidatus Marinamargulisbacteria bacterium', 447829)



#######################
# REGISTERING SAMPLES #
#######################

def make_mag_sample(sample_key, att_dict):
    blacklist = ['scientific_name', 'tax_id']

    header = """<SAMPLE alias="{name}" center_name="">
    <TITLE>{title}</TITLE>
    <SAMPLE_NAME>
      <TAXON_ID>{tax_id}</TAXON_ID>
      <SCIENTIFIC_NAME>{tax_name}</SCIENTIFIC_NAME>
      <COMMON_NAME></COMMON_NAME>
    </SAMPLE_NAME>
    <SAMPLE_ATTRIBUTES>
    """

    title = "MAG {name} from assembly {sample}".format(name = sample_key.replace("_megahit_metabat", ""), sample = att_dict['sample derived from'])
    formated_header = header.format(name = sample_key.replace("_megahit_metabat", ""), title = title, tax_id = tax2uncul[att_dict['scientific_name']][1], tax_name = tax2uncul[att_dict['scientific_name']][0])
    attributes = "".join([make_attributes(k,v, units.get(k, None)) for k,v in att_dict.items() if k not in blacklist])
    return formated_header + attributes + "</SAMPLE_ATTRIBUTES></SAMPLE>\n"

for i in list(range(1000, len(mag2md), 1000)):
    mag_sample_set_xml = """<?xml version="1.0" encoding="UTF-8"?>
    <SAMPLE_SET>
    {sample}
    </SAMPLE_SET>
    """.format(sample = "".join([make_mag_sample(k,v) for k,v in list(mag2md.items())[i:(i+1000)]]))


    with open(pjoin(template_folder, "mag_samples_{i}.xml").format(i=i), "w") as handle:
        handle.writelines(mag_sample_set_xml)

    mags_subm = curl_line.format(user = user, pwd = pwd,
                                   submission_file = submission_file,
                                   submission_type = "SAMPLE",
                                   url = url,
                                   data_file = pjoin(template_folder, "mag_samples_{i}.xml").format(i=i),
                                   out_file = pjoin(template_folder, "mag_samples_{i}.out.xml").format(i=i)
                                   )
    if not os.path.exists(pjoin(template_folder, "mag_samples_{i}.out.xml").format(i=i)):
        call(mags_subm, shell = True)

bin_access = {}
for i in list(range(0, len(mag2md), 1000)):
    with open(pjoin(template_folder, "mag_samples_{i}.out.xml").format(i=i)) as handle:
        root = ET.parse(handle).getroot()
    coass_accessions = {}
    data = [child.attrib for child in root if child.tag == "SAMPLE"]
    if not all(["accession" in dd for dd in data]):
        data = [[t.text for t in child if t.tag == 'ERROR'] for child in root if child.tag ==  "MESSAGES"]
        bin_access.update({ l.split('\"')[1] : l.split('\"')[-2] for l in data[0]})
    else:
        bin_access.update({dd['alias'] : dd['accession'] for dd in data})

for k, v in bin_access.items():
    mag2md[k.replace("_", "_megahit_metabat_")]['accession'] = v

###############
# SUBMIT BINS #
###############

bin_manifest = """STUDY {study_id}
SAMPLE {sample_id}
ASSEMBLYNAME freshwater MAG --- {bin_name}
ASSEMBLY_TYPE Metagenome-Assembled Genome (MAG)
COVERAGE 1
PROGRAM megahit
PLATFORM  Illumina NovaSeq 6000
MOLECULETYPE genomic DNA
FASTA {seq}
"""

for k,v in tqdm(mag2md.items()):
    ass = k.split("_")[0]
    if ass.startswith("Loc"):
        ass = ass.replace("-", "_")
    if ass[1] == "-":
        ass = ass.replace("-","")
    instanced_manifest = bin_manifest.format(study_id = project_id, sample_id = v['accession'], bin_name = k.replace("_megahit_metabat", ""), seq = pjoin("bins", ass, k,  k + ".fna.gz"))
    with open("data/manifests/bin_manifests/{}_binmanifest.txt".format(k), "w") as handle:
        handle.writelines(instanced_manifest)

bins_to_submit = [pjoin("data/manifests/bin_manifests/" + f) for f in os.listdir("data/manifests/bin_manifests/") if f.endswith("_binmanifest.txt")]

i = 0
while len(bins_to_submit) >0  and i < 150:
    i += 1
    print("try",i, ",", len(bins_to_submit), "bins left to submit")
    for j in list(range(0, len(bins_to_submit), 500)):
        binparallel_script = "parallel -j{threads} ./binhelper.sh {{}} ::: {samples}".format(threads = 20, samples = " ".join(bins_to_submit[j:(j+500)]))
        call(binparallel_script, shell=True)
    bins_to_submit = [pjoin("data/manifests/bin_manifests/" + f) for f in os.listdir("data/manifests/bin_manifests/") if f.endswith("_binmanifest.txt")]

for j in bins_to_submit:
    binparallel_script = "./binhelper.sh {}".format(j)
    call(binparallel_script, shell=True)

########################
# SUBMITTING SAG READS #
########################

with open("../../../data/data_submit/metadata/sags_to_db_200716.txt") as handle:
    sag2sample = handle.readlines()

sag2sample = {l.split()[0] : l[:-1].split()[1] for l in sag2sample[1:]}
sag2sample = {k : v.replace("_","-").replace("Day1", "Day2") for k,v in sag2sample.items()}

with open("/home/moritz/data/data_submit/metadata/sag_vs_reps.txt") as handle:
    anis = [(l.split()[0][:-4].split("/")[-1], l.split()[1][:-7].split("/")[-1], float(l.split()[2])) for l in handle.readlines()[1:]]
max_anis = {vv : max([vvv for vvv in anis if vvv[0] == vv], key = lambda x:x[2]) for vv in {v[0] for v in anis}}
sag2motu = {k.split("_") : v[1] if v[2] > 95 else None  for k,v in max_anis.items()}

with open("/home/moritz/data/data_submit/temp/checkm.txt") as handle:
    lines = [l.strip().split() for l in handle if l[0] not in ["-", "["]]
    sag_checkm = {l[0].split(".")[0] : {'completeness' : float(l[12]), 'contamination' : float(l[13]), 'strain_heterogeneity' : float(l[14])} for l in lines}

with open("/home/moritz/data/data_submit/temp/gtdbtk_classify/gtdbtk.bac120.summary.tsv") as handle:
    lines = [l.strip().split() for l in handle][1:]
    sag_gtdb = {l[0] : {'gtdbtk_taxonomy' : l[1]} for l in lines}

with open("/home/moritz/data/data_submit/temp/sourmash.tax") as handle:
    lines = [l.strip().split(",") for l in handle][1:]
    sag_sourmash = {l[0][:-4] : {'sourmash_taxonomy' : ";".join([pre if val =="" else val for pre, val in zip(prefix, l[2:]) ]) } for l in lines}


def get_genome_stats(sag):
    stats = {}
    sag_folder = "/home/moritz/data/data_submit/sags"
    fna = [str(s.seq) for s in  SeqIO.parse(pjoin(sag_folder,sag, sag + ".fna"), "fasta")]
    faa = [str(s.seq) for s in  SeqIO.parse(pjoin(sag_folder,sag, sag + ".faa"), "fasta")]
    stats['length'] = sum([len(l) for l in fna])
    stats['GC'] = sum([l.count("G") + l.count("C") for l in fna])/stats['length']
    stats['coding_density'] = 3*sum([len(l) for l in faa])/stats['length']
    stats['nb_contigs'] = len(fna)
    stats['mOTU'] = sag2motu.get(sag,sag2motu.get(sag[:-4]))
    stats['nb_proteins'] = len(faa)
    stats['path'] = pjoin("sags",sag + ".tar.gz")
    stats['sample'] = sag2sample[sag]
    stats.update(sag_checkm[sag])
    stats['sample derived from'] = sample_data[stats['sample']]['sample_accession']
    tt = sag_gtdb.get(sag,sag_gtdb.get(sag[:-4]))
    if tt :
        stats.update(tt)
    else :
        stats['gtdbtk_taxonomy'] = ";".join(prefix)
    tt = sag_sourmash.get(sag,sag_sourmash.get(sag[:-4]))
    if tt:
        stats.update(tt)
    stats['taxonomy'] = max([stats['gtdbtk_taxonomy'], stats['sourmash_taxonomy']], key = len)
    stats['taxonomy'] = ";".join([v for v in stats['taxonomy'].split(";") if len(v) > 3])
    return stats

sag2stats = {s : get_genome_stats(s) for s in tqdm(sag2sample)}

general_sag_md =   { 'sequencing method' : 'Illumina HiSeq',
                 'assembly software' : 'SPAdes',
                 'completeness software' : 'checkm',
                 'assembly quality' : 'Many fragments with little to no review of assembly other than reporting of standard assembly statistics',
                 'investigation type' : 'metagenome-assembled genome',
                 'binning parameters' : 'min_len:1500,min_bin_size:10000,maxP:93,minS:50',
                 'taxonomic identity marker' : 'gtdbtk followed by sourmash-lca, see publication',
                 'isolation_source' : 'freshwater',
                 'metagenomic source' : '449393',
                 'ENA-CHECKLIST' : 'ERC000048',
                 'sorting technology' : 'flow cytometric cell sorting',
                 'WGA amplification approach' : 'mda based',
                 'single cell or viral particle lysis approach' : 'enzymatic'
                 }

sag2md = {}
for k,v in tqdm(sag2stats.items()):
    #print("doing",k)
    sag2md[k] = {}
    sag2md[k].update(general_md)
    sag2md[k]['completeness score'] = v['completeness'] if v['completeness'] < 100 else 99.9
    sag2md[k]['contamination score'] = v['contamination']
    sag2md[k]['scientific_name'], sag2md[k]['tax_id'] = convert2ncbi(v['taxonomy']) if v['taxonomy'] != '' else convert2ncbi('not_classified')
    sag2md[k]['gtdb_taxonomy'] = v['taxonomy']
    sag2md[k]['metagenomic OTU'] = v['mOTU']
    ass_md = sample_data[v['sample']]
    sag2md[k].update(ass_md)
    sag2md[k]['sample derived from'] = ass_md['sample_accession']
    if 'sample_accession' in sag2md[k]:
        del sag2md[k]['sample_accession']
    if 'Run' in sag2md[k]:
        del sag2md[k]['Run']
    if 'Lake_code' in sag2md[k]:
            del sag2md[k]['Lake_code']


def make_sag_sample(sample_key, att_dict):
    blacklist = ['scientific_name', 'tax_id']

    header = """<SAMPLE alias="{name}" center_name="">
    <TITLE>{title}</TITLE>
    <SAMPLE_NAME>
      <TAXON_ID>{tax_id}</TAXON_ID>
      <SCIENTIFIC_NAME>{tax_name}</SCIENTIFIC_NAME>
      <COMMON_NAME></COMMON_NAME>
    </SAMPLE_NAME>
    <SAMPLE_ATTRIBUTES>
    """

    title = "SAG {name} from sample {sample}".format(name = sample_key, sample = att_dict['sample derived from'])
    formated_header = header.format(name = sample_key, title = title, tax_id = tax2uncul[att_dict['scientific_name']][1], tax_name = tax2uncul[att_dict['scientific_name']][0])
    attributes = "".join([make_attributes(k,v, units.get(k, None)) for k,v in att_dict.items() if k not in blacklist])
    return formated_header + attributes + "</SAMPLE_ATTRIBUTES></SAMPLE>\n"

sag_sample_set_xml = """<?xml version="1.0" encoding="UTF-8"?>
<SAMPLE_SET>
{sample}
</SAMPLE_SET>
""".format(sample = "".join([make_sag_sample(k,v) for k,v in sag2md.items()]))


with open(pjoin(template_folder, "sag_samples.xml"), "w") as handle:
    handle.writelines(sag_sample_set_xml)

mags_subm = curl_line.format(user = user, pwd = pwd,
                               submission_file = submission_file,
                               submission_type = "SAMPLE",
                               url = url,
                               data_file = pjoin(template_folder, "sag_samples.xml").format(i=i),
                               out_file = pjoin(template_folder, "sag_samples.out.xml").format(i=i)
                               )
call(mags_subm, shell = True)

sag_access = {}
with open(pjoin(template_folder, "sag_samples.out.xml").format(i=i)) as handle:
    root = ET.parse(handle).getroot()
    coass_accessions = {}
    data = [child.attrib for child in root if child.tag == "SAMPLE"]
    if not all(["accession" in dd for dd in data]):
        data = [[t.text for t in child if t.tag == 'ERROR'] for child in root if child.tag ==  "MESSAGES"]
        sag_access.update({ l.split('\"')[1] : l.split('\"')[-2] for l in data[0]})
    else:
        sag_access.update({dd['alias'] : dd['accession'] for dd in data})

for k, v in sag_access.items():
    sag2md[k]['accession'] = v





sag_read_manifest = """STUDY {study_id}
SAMPLE {sample_id}
NAME  microbes in stratified freshwaters - {library_name}
INSTRUMENT  HiSeq X Ten
LIBRARY_NAME  {library_name}
LIBRARY_SOURCE  GENOMIC SINGLE CELL
LIBRARY_SELECTION RANDOM
LIBRARY_STRATEGY  WGA
FASTQ {fwd}
FASTQ {rev}
"""
sag_raw_folder = "/home/moritz/data/data_submit/reads/sags/"
sag_title = lambda sample_key, att_dict : "SAG {name} from sample {sample}".format(name = sample_key, sample = att_dict['sample derived from'])

for k,v in tqdm(sag2md.items()):
    if not 'exp_id' in v:
        print("Doing", k)
        instanced_manifest = sag_read_manifest.format(study_id = project_id, sample_id = v['accession'], library_name = sag_title(k,v), fwd = pjoin(sag_raw_folder, k + "_R1.fastq.gz"), rev = pjoin(sag_raw_folder, k + "_R2.fastq.gz"))
        with open("data/manifests/{}_manifest.txt".format(k), "w") as handle:
            handle.writelines(instanced_manifest)

samples_to_submit = [pjoin("data/manifests/" + f) for f in os.listdir("data/manifests/") if f.endswith("_manifest.txt")]

for f in samples_to_submit:
    call("./helper.sh {sample}".format(sample = f), shell=True)

sag_manifest = """STUDY {study_id}
SAMPLE {sample_id}
ASSEMBLYNAME freshwater SAG --- {sag_name}
ASSEMBLY_TYPE Environmental Single-Cell Amplified Genome (SAG)
COVERAGE 1
PROGRAM SPAdes
PLATFORM  HiSeq X Ten
MOLECULETYPE genomic DNA
FASTA {seq}
"""

sag_folder = "/home/moritz/data/data_submit/sags/"

for k,v in tqdm(sag2md.items()):
    if not 'exp_id' in v:
        print("Doing", k)
        instanced_manifest = sag_manifest.format(study_id = project_id, sample_id = v['accession'], sag_name = k, seq = pjoin(sag_folder,k, k + ".fasta.gz"))
        with open("data/manifests/{}_sagmanifest.txt".format(k), "w") as handle:
            handle.writelines(instanced_manifest)

sags_to_submit = [pjoin("data/manifests/" + f) for f in os.listdir("data/manifests/") if f.endswith("_sagmanifest.txt")]

for f in sags_to_submit:
    call("./sag_helper.sh {sample}".format(sample = f), shell=True)



def get_assembly_ids(sample_id):
    try :
        sample_id = Entrez.read(Entrez.esearch(db="biosample", term = sample_id, idtype="acc"))['IdList'][0]
        sample_accession = Entrez.read(Entrez.esummary(db="biosample", id = sample_id))['DocumentSummarySet']['DocumentSummary'][0]['Accession']
        assembly_ids = Entrez.read(Entrez.esearch(db="assembly", term = sample_accession, idtype="acc"))['IdList'][0]
        assembly_accession = Entrez.read(Entrez.esummary(db="assembly", id = assembly_ids))['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession']
    except:
        return None
    return assembly_accession

def get_sample_readinfo(id):
    sample_ids = Entrez.read(Entrez.esearch(db="sra", term = id, idtype="acc"))['IdList']
    out = {}
    for sample_id in sample_ids:
        xml = Entrez.read(Entrez.esummary(db="sra", id = sample_id))[0]['ExpXml']
        root = ET.XML( "<WTF>" + xml + "</WTF>")
        stats = [ll for ll in [l for l in root if l.tag =="Summary"][0] if ll.tag == "Statistics"][0].attrib
        sra_id = [l for l in root if l.tag =="Submitter"][0].attrib['acc']
        title = [ll for ll in [l for l in root if l.tag =="Summary"][0] if ll.tag == "Title"][0].text
        stats['title'] = title
        stats['sample_id'] = id
        out.update({sra_id : stats})
    return out

def get_MG_ids(sample_id):
    try :
        sample_ids = Entrez.read(Entrez.esearch(db="biosample", term = sample_id, idtype="acc", retmax=100))['IdList']
        assemblies = []
        nassembly = []
        for sample_id in tqdm(sample_ids):
            sample_accession = Entrez.read(Entrez.esummary(db="biosample", id = sample_id))['DocumentSummarySet']['DocumentSummary'][0]['Accession']
            assembly_ids = [ e for e in Entrez.read(Entrez.esearch(db="assembly", term = sample_accession))['IdList'] if len(e) > 0]
            if len(assembly_ids) > 0 :
                assemblies += [Entrez.read(Entrez.esummary(db="assembly", id = assembly_ids))['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession']]
            else :
                nassembly += [sample_accession]
    except:
        return None
    return assembly_accession


def motu_stats(motu) :
    tax = consensus_tax(motu)
    mean_anis = mean([l[2] for l in motu_full[motu]['ANIs']])
    mean_good_anis = mean([l[2] for l in motu_full[motu]['ANIs'] if is_good(l[0]) and is_good(l[1])])
    mean_decent_anis = mean([l[2] for l in motu_full[motu]['ANIs'] if is_decent(l[0]) and is_decent(l[1])])
    nb_genomes = len(motu_full[motu]['genomes'])
    nb_goods = sum([is_good(g['name'])  for g in motu_full[motu]['genomes']])
    nb_decents = sum([is_decent(g['name'])  for g in motu_full[motu]['genomes']])
    return {'mean_ANI' : mean_anis, 'consensus_tax' : tax, 'mean_good_ANIs' : mean_good_anis, 'mean_decent_ANIs' :mean_decent_anis, 'nb_genomes' : nb_genomes, 'nb_good_genomes' : nb_goods, 'nb_decent_genomes' : nb_decents}


sag_folder = "/home/moritz/data/data_submit/assemblies"
ass2md = { ass : get_assembly_stats(ass) for ass in tqdm([o for o in os.listdir(sag_folder) if o.endswith(".fna.bz2")])}

def get_assembly_stats(ass):
    stats = {}
    with bz2.open(pjoin(sag_folder, ass ), "rt") as handle:
        fna = [str(s.seq) for s in  SeqIO.parse(handle, "fasta")]

    stats['length'] = sum([len(l) for l in fna])
    stats['GC'] = sum([l.count("G") + l.count("C") for l in fna])/stats['length']
    stats['nb_contigs'] = len(fna)
    stats['path'] = pjoin("assembly",ass + ".tar.gz")

    sample = ass[:-8]
    sample = sample.replace("-","_") if sample.startswith("Loc") else sample

    if sample in sample2sra:
        stats['type'] = "single_sample_assembly"
        stats['libraries'] = sample2sra[sample]
    else :
        with open("metadata/coassemblies/" + sample + ".txt") as handle:
            samples = [l[:-1].replace("-","_") if l[:-1].startswith("Loc") else l[:-1]  for l in handle]
        stats['type'] = "coassembly_assembly"
        stats['libraries'] = ";".join([sample2sra[l] for l in samples])
    return stats
