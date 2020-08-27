import sys
from subprocess import Popen, PIPE, call
import os
from os.path import join as pjoin
from tqdm import tqdm
import xml.etree.ElementTree as ET
import re
import pandas
from numpy import logical_and

user = "$USER_NAME"
pwd = "$PASSWD"

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




#PROJ ID

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

    header = """<SAMPLE alias="{name}" center_name="">
    <TITLE>{title}</TITLE>
    <SAMPLE_NAME>
      <TAXON_ID>{tax_id}</TAXON_ID>
      <SCIENTIFIC_NAME>{tax}</SCIENTIFIC_NAME>
      <COMMON_NAME></COMMON_NAME>
    </SAMPLE_NAME>
    <SAMPLE_ATTRIBUTES>
    """
    #1851193	pond metagenome
    #449393 freshwater metagenome
    blacklist = ['ncbi_taxid', 'taxonomy']
    title = "{type} {name} sample {sample}".format(type = att_dict['environment (feature)'], name = att_dict['site name'], sample = sample_key)
    formated_header = header.format(name = sample_key, title = title, tax_id = att_dict['ncbi_taxid'], tax = att_dict['taxonomy'])
    attributes = "".join([make_attributes(k,v, units.get(k, None)) for k,v in att_dict.items() if k not in blacklist])
    return formated_header + attributes + "</SAMPLE_ATTRIBUTES></SAMPLE>\n"


sample_set_xml = """<?xml version="1.0" encoding="UTF-8"?>
<SAMPLE_SET>
{sample}
</SAMPLE_SET>
""".format(sample = "".join([make_sample(k,v) for k,v in sample_data.items()]))

with open(pjoin(template_folder, "fix_all_libs.xml"), "w") as handle:
    handle.writelines(sample_set_xml)


fix_subm = curl_line.format(user = user, pwd = pwd,
                               submission_file = pjoin(template_folder, 'fix.xml'),
                               submission_type = "SAMPLE",
                               url = url,
                               data_file = pjoin(template_folder, "fix_all_libs.xml"),
                               out_file = pjoin(template_folder, "fix_all_libs.out.xml")
                   )

call(fix_subm, shell = True)

with open(pjoin(template_folder,"all_samples.out.xml")) as handle:
    root = ET.parse(handle).getroot()

data = [child.attrib for child in root if child.tag == "SAMPLE"]
for dd in data:
    sample_data[dd['alias']]['sample_accession'] = dd['accession']

######################
# PREPING ASSEMBLIES #
######################

raw_assemblies_dir = pjoin(raw_folder, "assemblies")
coassemblies_dir = pjoin(raw_folder, "metadata/coassemblies")
coass_list = [c[:-4] for c in os.listdir(coassemblies_dir)]

single_sample_assemblies = [f[:-8] for f in os.listdir(raw_assemblies_dir) if f[:-8] not in coass_list and (f[:-8] in sample_data or "Loc" in f)  ]

def fix_loc_samples(s):
    if not s.startswith("Loc0") :
        return s
    else :
        return k.replace("_","-")

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

    header = """<SAMPLE alias="{name}" center_name="">
    <TITLE>{title}</TITLE>
    <SAMPLE_NAME>
      <TAXON_ID>{tax_id}</TAXON_ID>
      <SCIENTIFIC_NAME>{tax}</SCIENTIFIC_NAME>
      <COMMON_NAME></COMMON_NAME>
    </SAMPLE_NAME>
    <SAMPLE_ATTRIBUTES>
    """
    blacklist = ['ncbi_taxid', 'taxonomy']

    title = "mock sample for coassembly {name} sample of samples: {sample}".format(name = sample_key, sample = att_dict['sample derived from'])
    formated_header = header.format(name = sample_key, title = title, tax_id = att_dict['ncbi_taxid'], tax = att_dict['taxonomy'])
    attributes = "".join([make_attributes(k,v, units.get(k, None)) for k,v in att_dict.items() if k not in blacklist])
    return formated_header + attributes + "</SAMPLE_ATTRIBUTES></SAMPLE>\n"

coasssample_set_xml = """<?xml version="1.0" encoding="UTF-8"?>
<SAMPLE_SET>
{sample}
</SAMPLE_SET>
""".format(sample = "".join([make_coassample(k,v) for k,v in coass2md.items()]))

with open(pjoin(template_folder, "fix_all_coass_samples.xml"), "w") as handle:
    handle.writelines(coasssample_set_xml)

# REGISTERING SAMPLES

all_samples_subm = curl_line.format(user = user, pwd = pwd,
                               submission_file = pjoin(template_folder, 'fix.xml'),
                               submission_type = "SAMPLE",
                               url = url,
                               data_file = pjoin(template_folder, "fix_all_coass_samples.xml"),
                               out_file = pjoin(template_folder, "fix_all_coass_samples.out.xml")
                               )

call(all_samples_subm, shell = True)

with open(pjoin(template_folder,"all_coass_samples.out.xml")) as handle:
    root = ET.parse(handle).getroot()

coass_accessions = {}
data = [child.attrib for child in root if child.tag == "SAMPLE"]
for dd in data:
    coass_accessions[dd['alias']] = dd['accession']

################
# PREPING MAGS #
################
from ete3 import NCBITaxa
ncbi = NCBITaxa()
mag_md = pandas.read_csv(mag_sheet, index_col=0)

general_md =   { 'sequencing method' : 'Illumina NovaSeq',
                 'assembly software' : 'megahit',
                 'completeness software' : 'checkm',
                 'binning software' : 'metabat',
                 'assembly quality' : 'Many fragments with little to no review of assembly other than reporting of standard assembly statistics',
                 'investigation type' : 'metagenome-assembled genome',
                 'binning parameters' : 'min_len:1500,min_bin_size:10000,maxP:93,minS:50',
                 'taxonomic identity marker' : 'gtdbtk followed by sourmash-lca, see publication',
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
    hits = [v for k, v in gtdb2ncbi.items() if k.startswith(tax)]
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
    mag2md[k]['isolation_source'] = " ".join(ass_md['taxonomy'].split()[:-1])
    mag2md[k]['metagenomic source'] = ass_md['ncbi_taxid']
    del mag2md[k]['ncbi_taxid']
    del mag2md[k]['taxonomy']
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


    with open(pjoin(template_folder, "fix_mag_samples_{i}.xml").format(i=i), "w") as handle:
        handle.writelines(mag_sample_set_xml)

    mags_subm = curl_line.format(user = user, pwd = pwd,
                                   submission_file = pjoin(template_folder, 'fix.xml'),
                                   submission_type = "SAMPLE",
                                   url = url,
                                   data_file = pjoin(template_folder, "fix_mag_samples_{i}.xml").format(i=i),
                                   out_file = pjoin(template_folder, "fix_mag_samples_{i}.out.xml").format(i=i)
                                   )
    if not os.path.exists(pjoin(template_folder, "fix_mag_samples_{i}.out.xml").format(i=i)):
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

with open("coords2.csv", "w") as handle:
    tt = {(v['site name'],str(v['geographic location (longitude)']), str(v['geographic location (latitude)']))  for k,v in sample_data.items()}
    handle.writelines([",".join(t) + "\n" for t in tt])
