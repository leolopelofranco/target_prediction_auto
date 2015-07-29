# This script simply converts GO IDs of the mRNA targets to their respective Gene symbols

import sys
import re

from Bio import Entrez

# *Always* tell NCBI who you are
Entrez.email = "email@gmail.com"

d = open('dic.txt', 'r')
dic = d.read()
diction = dic.split('\n')

l = open('viral_miRNA_vs_hiv_linked_proteins_greater_1.sif', 'r')
lnc_dic = l.read()
to_be_converted = lnc_dic.split('\n')

diction_ref = []
convert_this = []

for i in diction:
    a = i.split()
    if len(a)>0:
        x ={}
        x['gi_number'] = a[0]
        x['symbol'] = a[1]
        diction_ref.append(x)

for i in to_be_converted:
    line = {}
    a = i.split()
    if len(a)>0:
        line['mir'] = a[0]
        line['hits'] = a[1]
        line['gi_number'] = a[2]
        convert_this.append(line)

no_symbol =  []

for i in convert_this:
    ctr = 0
    for z in diction_ref:
        if i['gi_number'] == z['gi_number']:
            ctr = 1
            print i['mir'] + '\t' + i['hits'] + '\t' + z['symbol']
    if ctr == 0:
        miRNA = {}
        miRNA['family'] = i['mir']
        miRNA['hits'] = i['hits']
        miRNA['target'] = i['gi_number']
        no_symbol.append(miRNA)
    ctr = 0

######


all_proteins = no_symbol

all_lnc_proteins = []


a = []
b = []


for i in all_proteins:
    a.append(i['target'])

id_list = list(set(a))

def chunks(l, n):
    n = max(1, n)
    return [l[i:i + n] for i in range(0, len(l), n)]


def converter(ids):

    chunky = chunks(ids, 400)
    all_translators = []
    translator = []
    exceptions = []
    all_exceptions = []

    ctr = 1
    for i in chunky:
        handle = Entrez.elink(dbfrom="nuccore", db="gene", id=i)
        result = Entrez.read(handle)
        handle.close()


        for each_record in result:
            try:
                mrna = {}
                print each_record
                mrna['gene_id'] = each_record["LinkSetDb"][0]["Link"][0]["Id"]
                mrna['mrna_id'] = each_record['IdList'][0]
                translator.append(mrna)

                """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
                submit the data to NCBI) and esummary to retrieve the information.
                Returns a list of dictionaries with the annotations."""


            except IndexError:
                exceptions.append(each_record['IdList'][0])

        all_translators.append(translator)
        all_exceptions.append(exceptions)
        translator = []
        exceptions = []


    flatten(all_translators, all_exceptions)

def flatten(a, b):
    gene_ids = []
    flattened = [item for sublist in a for item in sublist]
    flattened_exceptions = [item for sublist in b for item in sublist]

    for i in flattened:
        gene_ids.append(i['gene_id'])


    gene_ids = list(set(gene_ids))
    flattened_exceptions = list(set(flattened_exceptions))

    if len(flattened_exceptions) > 0:
        annotations_exceptions = retrieve_annotation(flattened_exceptions, 'nuccore')
        symbols_exceptions(annotations_exceptions)

    if len(gene_ids) > 0:
        annotations = retrieve_annotation(gene_ids, 'gene')
        symbols(annotations, flattened)

    gene_ids = []
    flattened_exceptions = []
    flattened = []
    annotations_exceptions = []
    annotations = []

def compare(s,f):
    for i in f:
        for a in s:
            if i['gene_id'] == a['gene_id']:
                print 'hl'
                #print i['mrna_id'] + '\t' + a['gene_symbol']

def compare_exceptions(a):
    for i in a:
        print 'hl'
        #print i['mrna_id'] + '\t' + i['symbol']

def symbols(annotation, flattened):
    symbols = []
    annotation = annotation['DocumentSummarySet']['DocumentSummary']

    for gene_data in annotation:
        g = {}
        g['gene_id'] = gene_data.attributes['uid']
        g['gene_symbol'] = gene_data["NomenclatureSymbol"]
        g['gene_name'] = gene_data["Description"]
        symbols.append(g)
    compare(symbols, flattened)
    symbols = []

def symbols_exceptions(annotation):
    symbols_exceptions= []

    for mrna_data in annotation:
        h = {}
        h['mrna_id'] = mrna_data["Id"]
        h['mrna_name'] = mrna_data["Title"]
        m = re.findall(r"\((.*?)\)", h['mrna_name'])

        if len(m)>1:
            h['symbol'] = m[1]
        else:
            h['symbol'] = m[0]

        symbols_exceptions.append(h)
    compare_exceptions(symbols_exceptions)
    symbols_exceptions = []

def retrieve_annotation(id_list, dbs):

    annotations = []

    """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
    submit the data to NCBI) and esummary to retrieve the information.
    Returns a list of dictionaries with the annotations."""

    request = Entrez.epost(dbs,id=",".join(id_list))

    try:
        result = Entrez.read(request)
    except RuntimeError as e:
        #FIXME: How generate NAs instead of causing an error with invalid IDs?
        print "An error occurred while retrieving the annotations."
        print "The error returned was %s" % e
        sys.exit(-1)

    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.esummary(db=dbs, webenv=webEnv, query_key =
            queryKey)
    annotations = Entrez.read(data)

    return annotations

converter(id_list)
