#!/usr/bin/env python

#functions for preparing input data in detecting gene origins
#note the types of input files and their separators, input given by config.ini via master script

import sys
import os
import pandas as pd
from ete3 import NCBITaxa
ncbi = NCBITaxa()

def mkdir_and_cd(dir_name):
    try:
        os.mkdir(dir_name)
    except FileExistsError:
        pass
    os.chdir(dir_name)

def get_oglist(infile):
    with open(infile, 'r') as f:
        oglist = f.readlines()
        oglist = [x.rstrip("\n").split(':')[0] for x in oglist if x != ""]
    return oglist

def get_translation_abbr_taxids(infile):
    species_translation={}
    with open(infile, 'r') as f:
        for line in f:
            line = line.rstrip("\n")
            species_translation[line.split(",")[0]]=int(line.split(",")[1])
    return species_translation

def get_full_taxonomies(infile):
    """Collects the complete taxonomies of GTDB assemblies and stores them in a dictionary with assembly id as key and the taxonomy as a value in the type of a list with 7 items"""
    full_taxonomy={}
    with open(infile, 'r') as f:
        for line in f:
            line = line.rstrip("\n")
            assembly_id = line.split()[0]
            taxonomy = line.split()[1].split(";")
            full_taxonomy[assembly_id]=taxonomy
    return full_taxonomy

def get_full_taxonomies_all_levels(taxonomy):
    """Uses the full taxonomies to find the ancestors on all levels and uses the taxon name itself (without the rank indicator) as a key"""
    taxonomy_extended = {}
    for k, v in taxonomy.items():
        for i, l in enumerate(v):
            l = l.split("__")[1]
            new_v = v[0:i]
            taxonomy_extended[l] = new_v
    return taxonomy_extended

def get_metadata(infile):
    species_names={}
    phyla={}
    data_types = {}
    contig_collections = {}
    with open(infile, 'r') as f:
        for line in f:
            entries = line.rstrip("\n").split("\t")
            species_names[entries[0]] = entries[1]
            phyla[entries[0]] = entries[2]
            data_types[entries[0]] = entries[3]
            contig_collections[entries[0]] = entries[4]
    return species_names, phyla, data_types, contig_collections

def get_contaminants(infile):
    sequence_source = {}
    with open(infile, 'r') as f:
        for line in f:
            entries = line.rstrip("\n").split("\t")
            sequence_source[entries[0]] = entries[1]
    return sequence_source

def get_supergroups_taxids(inlist):
    supergroups = inlist.split(",")
    supergroups_dict = {ncbi.get_name_translator([s])[s][0]:s for s in supergroups}
    return supergroups_dict

def get_vertical_sequences(seq_results):
    seq_results_df = pd.read_csv(seq_results, index_col='sequence')
    vertical_sequences = seq_results_df[seq_results_df.sequence_origin == 'vertical'].index.tolist()
    del seq_results_df
    return vertical_sequences


def insert_lca_prediction_data(sequences, lca, dataset_prediction):
    """Append the lca feature to the dataset for origin prediction, using the binary pattern"""
    column = f"rhizarian_node_lca_{lca}"
    dataset_prediction.loc[sequences, column] = 1
    return dataset_prediction

