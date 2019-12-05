import matplotlib as mpl

import numpy as np
import pandas as pd

import requests

from io import StringIO
from PIL import Image
import re


def get_pw_info(pw_id):
    """Parse KEGG GET API output to get gene and compound information so far. Input is just \t or \n delimited
    raw string, output as a dictionary. A bit annoying that the API doesn't already serve json or xml.

    PARAMS
    ------
    pw_id: str; KEGG ID of pathway

    RETURNS
    -------
    out_dict: dictionary of output, where keys are API section headers, e.g. "genes", "compounds", and values
    are subdictionaries with key-value pairs of feature_ids and feature descriptions.
    """
    url = "http://rest.kegg.jp/get/"+pw_id
    r = requests.get(url)

    contents = str(r.content).split("\\n")
    gene_contents = []
    cpd_contents = []

    gene_reader_bool = 0
    cpd_reader_bool = 0

    for line in contents:
        # Set reader bools
        if "GENE" in line:
            gene_reader_bool = 1
        if "COMPOUND" in line:
            cpd_reader_bool = 1
            # Switch off gene_reader_bool
            if gene_reader_bool == 1:
                gene_reader_bool = 0
        if "REFERENCE" in line:
            # set all reader_bools to 0 (stop reading)
            gene_reader_bool = 0
            cpd_reader_bool = 0

        # Read contents
        if gene_reader_bool == 1:
            gene_contents.append(line)

        if cpd_reader_bool == 1:
            cpd_contents.append(line)

    # Wrangle some lists
    gene_contents[0] = gene_contents[0].replace("GENE", "    ")
    gene_contents = [re.sub('   +', '', x) for x in gene_contents]
    gene_contents = [x.split("  ") for x in gene_contents]

    cpd_contents[0] = cpd_contents[0].replace("COMPOUND", "        ")
    cpd_contents = [re.sub('   +', '', x) for x in cpd_contents]
    cpd_contents = [x.split("  ") for x in cpd_contents]

    gene_dict = {}
    for i in range(len(gene_contents)):
        gene_dict[gene_contents[i][0]] = gene_contents[i][1]

    cpd_dict = {}
    for i in range(len(cpd_contents)):
        cpd_dict[cpd_contents[i][0]] = cpd_contents[i][1]

    out_dict = {}
    out_dict["genes"] = gene_dict
    out_dict["compounds"] = cpd_dict

    return out_dict


def get_kegg_pathway_mapped_png(pw, gene_id_dict, fn_out, fn_out_format=".png", verbose=True):
    """Given a KEGG pathway ID and a dictionary of features with hexcode colours,
    Make POST call to KEGG and retrieve .png of data-mapped KEGG pathway.
    Note: this executes just fine even if certain input features that are passed through the constructed url
    do not appear in the specified pw.

    PARAMS
    ------
    pw: str; ID of KEGG Pathway, e.g. "hsa00010"
    gene_id_dict: dictionary of gene feature colours. Keys are any KEGG-recognized ID, e.g.
    EntrezIDs, KEGG compound IDs, values are hex colour codes, e.g. '#ff0000'.
    fn_out: str; output filename
    fn_out_format: str; file format of output filename
    verbose: bool; verbosity parameter.

    RETURNS
    -------
    Does not have an object returned; saves the Image.img output to file instead.
    """
    # input sanity check
    n_keys = len(gene_id_dict.keys())
    n_uq_keys = len(np.unique(list(gene_id_dict.keys())))
    if n_keys != n_uq_keys:
        print("WARNING: %s duplicate keys (features) detected in input dictionary" % (n_keys - n_uq_keys))

    # Form url for GET request
    multiquery_ls = ["&multi_query="]
    for k in list(gene_id_dict.keys()):
        new_str = str(k)+"+"+gene_id_dict[k].replace("#", "%23")
        multiquery_ls.append(new_str)
    multi_query_str = "%0d%0a".join(multiquery_ls)

    url = "https://www.kegg.jp/kegg-bin/show_pathway?map="+pw+multi_query_str
    if verbose:
        print("Mapping %s features to %s" % (len(gene_id_dict.keys()), pw), end="...")
    r = requests.get(url)
    if r.status_code == 200:
        if verbose:
            print(r.status_code)
        ms = re.finditer('src="(/tmp/mark_pathway[^"]*.png)"', r.text)
        m = list(ms).pop()

        # proc another GET request to pull the png
        img = Image.open(
            requests.get("http://www.kegg.jp%s" % m.group(1), stream=True).raw
        )

        #img_w, img_h = img.size
        img.save(fn_out+fn_out_format)
