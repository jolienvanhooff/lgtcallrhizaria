#!/usr/bin/env python3

import os
import argparse
import numpy as np
import pandas as pd
from configparser import ConfigParser

from Bio import SeqIO
from data_preparations import mkdir_and_cd, get_full_taxonomies, get_full_taxonomies_all_levels, get_supergroups_taxids, get_translation_abbr_taxids, get_metadata, get_contaminants, get_vertical_sequences, insert_lca_prediction_data
from orthogroup_examinations import *
from gene_tree_operations import open_species_tree, get_phylum_domain_lineage, taxid_to_name, verify_contamination, project_events, create_event_table, visualise_species_tree
from ete3 import TreeStyle
ts = TreeStyle()

from sklearn.ensemble import GradientBoostingClassifier
import joblib

# Parse arguments
parser = argparse.ArgumentParser(description = "This script predicts the origins of rhizarian genes without a phylogeny")
parser.add_argument("-l", metavar = "list", type=str, help = "list of OGs")
parser.add_argument("-i", metavar = "inputdir", type=str, help = "directory containing the FASTA files of OGs that should be analysed")
parser.add_argument('-o', metavar = 'outputdir', type=str, default='output_predict', help='directory for output files (DEFAULT: output_predict)')
parser.add_argument('-c', metavar = 'config', type=str, default='config_files/predict_config.ini', help='path to configuration file (DEFAULT: config_files/predict_config.ini)')
args = parser.parse_args()


# Parse configuration settings from the config file in the current working directory
config_object = ConfigParser()
config_object.read(args.c)
settings = config_object["DEFAULT"]

# Prepare arguments and configurations
#oglist = get_oglist(args.l)
oglist = [l.rstrip("\n") for l in open(args.l, 'r').readlines()]
inputdir = os.path.abspath(args.i)
outputdir = args.o
mkdir_and_cd(outputdir)

# Extract parameters from the configuration file
taxonomy_prokaryotes = get_full_taxonomies(settings["taxonomy_prokaryotes"]) #dict
taxonomy_prokaryotes_extended = get_full_taxonomies_all_levels(taxonomy_prokaryotes) #dict
eukaryotic_supergroups = get_supergroups_taxids(settings["supergroups_eukaryotes"]) #dict
sar_species_translation = get_translation_abbr_taxids(settings["taxonomy_identifiers_sar"]) #dict
species_names, phyla, data_types, contig_collections = get_metadata(settings["data_rhizaria"]) #dict
monophyly_cutoff = int(settings["monophyly_cutoff"]) #maximum number of sequences for a single rhizarian origin to be assumed
newick_rhizaria = settings["newick_rhizaria"] #bifurcating species phylogeny
threads = int(settings["threads"]) # number of threads to be used (e.g. in diamond runs)
evaluecutoff = float(settings["evaluecutoff"]) 
interspersing_proportion = float(settings["interspersing_proportion"]) #maximum proportion of sequences not belonging to clade of interest and/or majority clade
diamond_path = settings["diamond_path"]
do_contamination_check = True if settings["contamination_check"] == 'True' else False
contamination_sequences = get_contaminants(settings["contamination_signal_diamond"])
do_genome_check = True if settings["genome_check"] == 'True' else False
vertical_sequences = get_vertical_sequences(settings["genome_check_information"])
mcl_inflation_factor=float(settings["mcl_inflation_factor"])
dataset_prediction = pd.read_csv(settings['dataset_prediction'], index_col='sequence')
classifier = joblib.load(settings['classifier'])
cluster_unification_strategy = settings['unification']
mean_probability_threshold = float(settings['mean_probability_threshold'])

# Load Rhizaria species tree as PhyloTree object with events as attributes
species_tree = open_species_tree(newick_rhizaria)

# Create or load the output dataframes
if os.path.isfile(f'results_per_og.csv') is False:
    results_per_og = pd.DataFrame(columns=["orthogroup", "in_fasta", "origin_inferred", "viruses", "acceptors"])
    results_per_acquisition = pd.DataFrame(columns=["orthogroup", "rhizarian_node", "rhizarian_node_sequences", "tree", "tree_parents", "tree_reconciled", "data_check", "contamination_check", "lca", "origin", "donor", "donor_phylum", "donor_domain", "inference", "median_distance_rhizarian_node", "median_distance_sister1", "median_distance_sister2", "duplication", "loss", "acceptor_proportion_node", "virus_proportion_og", "prediction_probability", "genome_check"])
    prediction_results_sequences = pd.DataFrame(columns=['sequence', 'orthogroup', 'rhizarian_node','rhizarian_node_origin_predicted_seq_level', 'rhizarian_node_origin_predicted_node_level'])
    prediction_results_sequences.set_index('sequence', inplace=True)
else:
    results_per_og = pd.read_csv("results_per_og.csv")
    results_per_acquisition = pd.read_csv("results_per_acquisition.csv")
    prediction_results_sequences = pd.read_csv("prediction_results_sequences.csv", index_col='sequence')

# Analysis per OrthoGroup+
for i, og in enumerate(oglist):
    # Print the so far collected results
    results_per_og.to_csv("results_per_og.csv", index=False)
    results_per_acquisition.to_csv("results_per_acquisition.csv", index=False)
    prediction_results_sequences.to_csv("prediction_results_sequences.csv", index=True)
    # Check if the current og has already been investigated - if so do not study it   
    if og in results_per_og.orthogroup.values:
        continue 
    print("\n\n"+og)
    fasta = f'{inputdir}/{og}.fa'
    results_per_og = results_per_og.append({"orthogroup" : og}, ignore_index=True)
    if not os.path.exists(fasta):
        print(f"no fasta file found")
        continue
    results_per_og.loc[i, ['in_fasta']] = fasta
    results_per_og.loc[i, ['origin_inferred']] = 0
    og_records = SeqIO.index(fasta, 'fasta')
    seq_species, seq_taxonomy, seq_group, group_members = obtain_og_contents(og_records, taxonomy_prokaryotes, sar_species_translation)

    # Remove viruses and store the original tree
    filename = os.path.basename(fasta)
    og_records, viruses, virus_proportion_og = delete_group(og_records, seq_group, 'viruses', filename)
    print(f"viruses: {viruses}")
    results_per_og.loc[i, ['viruses']] = viruses

    # Remove the unknown leaves
    og_records, _, _ = delete_group(og_records, seq_group, 'unknown', filename)

    rhizaria_sequences = group_members['rhizarians']
    if len(rhizaria_sequences) == 0:
        print(f"no rhizarian sequences in og+")
        continue
    rhizaria_proportion = len(rhizaria_sequences) / len(og_records)

    # If the vast majority of OG is Rhizarian, infer a Rhizaria origin 
    if (1 - rhizaria_proportion) <= interspersing_proportion:
        # Perform data or contamination check: either assess whether sequences from a single transcriptome species are contaminants ('contamination_check') or perform a regular 'data_check'
        data_check, contamination_check = verify_contamination(rhizaria_sequences, contamination_sequences, data_types, do_contamination_check)
        if (do_contamination_check == True and contamination_check == False) or (do_contamination_check == False and data_check == False):
            results_per_acquisition = results_per_acquisition.append({"orthogroup" : og, "rhizarian_node": "anc.1", "rhizarian_node_sequences": ';'.join(rhizaria_sequences), "data_check": data_check, "contamination_check": contamination_check}, ignore_index=True)
            continue
        lca = find_lca(rhizaria_sequences, species_tree, species_names)
        origin = 'invention'
        acceptors = False
        acceptor_proportion_node = 0
        if rhizaria_proportion < 1:
            acceptors = True, 
            acceptor_proportion_node = 1 - rhizaria_proportion
        print(f"acceptors: {acceptors}")
        # Report results
        results_per_acquisition = results_per_acquisition.append(
            {
            "orthogroup" : og, 
            "rhizarian_node" : "anc.1", 
            "rhizarian_node_sequences" : ";".join(rhizaria_sequences),
            "data_check": data_check,
            "contamination_check": contamination_check,
            "lca" : lca, 
            "origin": origin, 
            "inference": 'rhizaria_dominance',
            "acceptor_proportion_node": acceptor_proportion_node,
            "virus_proportion_og": virus_proportion_og
            }, 
            ignore_index=True)
        project_events(species_tree, lca, origin, None)
        results_per_og.loc[i, ['acceptors']] = acceptors
        results_per_og.loc[i, ['origin_inferred']] = 1
    
    # If there is a considerable number of non-Rhizaria sequences, use the HGT index
    else: 
        # Create clusters within the OG
        mkdir_and_cd(f"{og}_clusters")
        mcl_clusters = perform_mcl_clustering(og, og_records, threads, evaluecutoff, diamond_path=diamond_path, inflation_factor=mcl_inflation_factor)
        # If somehow mcl failed, continue to the next OG
        if mcl_clusters == None:
            os.chdir("../")
            continue
        mcl_clusters_rhizaria = retrieve_subclusters_rhizaria(mcl_clusters, rhizaria_sequences, og_records)
        # Continue only if there is at least one with Rhizaria
        if len(mcl_clusters_rhizaria) == 0:
            sys.exit(f"{og}: no subclusters with rhizarian sequences found!")
        origin_inferred = 0
        for subcluster, subcluster_sequences in mcl_clusters_rhizaria.items():
            rhizarian_node = subcluster 
            rhizarian_node_sequences = [sequence for sequence in subcluster_sequences if sequence in rhizaria_sequences]
            # Perform data or contamination check: either assess whether sequences from a single transcriptome species are contaminants ('contamination_check') or perform a regular 'data_check'
            data_check, contamination_check = verify_contamination(rhizarian_node_sequences, contamination_sequences, data_types, do_contamination_check)
            if (do_contamination_check == True and contamination_check == False) or (do_contamination_check == False and data_check == False):
                results_per_acquisition = results_per_acquisition.append({"orthogroup" : og, "rhizarian_node": rhizarian_node, "rhizarian_node_sequences": ';'.join(rhizarian_node_sequences), "data_check": data_check, "contamination_check": contamination_check}, ignore_index=True)
                continue
            subcluster_fasta = f"{og}_{rhizarian_node}_and_outgroup.fa"
            create_output_fasta(og_records, subcluster_sequences, subcluster_fasta)
            # Create a new recorddict for the subcluster
            record_dict_subcluster = SeqIO.index(subcluster_fasta, 'fasta')
            # Append the LCA to the prediction dataset 
            lca = find_lca(rhizarian_node_sequences, species_tree, species_names)
            dataset_prediction = insert_lca_prediction_data(rhizarian_node_sequences, lca, dataset_prediction)
            # Make predictions for the ancestral rhizarian node of each sequence
            dataset_prediction_subcluster = dataset_prediction[dataset_prediction.index.isin(rhizarian_node_sequences)]
            prediction_results_sequences, predictions_dict, probabilities_dict, classes = classify_sequence_origins(dataset_prediction_subcluster, classifier, prediction_results_sequences)
            # Integrate the results of the individual sequences of the cluster into a single origin 
            origin, donor, donor_phylum, donor_domain, genome_check = np.nan, np.nan, np.nan, np.nan, np.nan
            if cluster_unification_strategy == "majority":
                predictions_list = list(predictions_dict.keys())
                origin = max(set(predictions_list), key = predictions_list.count)
            elif cluster_unification_strategy == "mean_probability":
                origin = get_origin_mean_class_probabilities(classes, probabilities_dict, mean_probability_threshold)
            # Don't analyse / report this cluster any further if none of the origins meets the threshold
            if origin == "ND":
                results_per_acquisition = results_per_acquisition.append({"orthogroup" : og, "rhizarian_node": rhizarian_node, "rhizarian_node_sequences": ';'.join(rhizarian_node_sequences), "data_check": data_check, "contamination_check": contamination_check, "prediction_probability": False}, ignore_index=True)
                print(f"None of the origins meets the mean probability threshold")
                continue
            # Check the genome if predicted to be lateral
            if origin == "lateral" and do_genome_check == True and lca in data_types:
                if data_types[lca] == 'genome':
                    protein_contig = pd.read_csv(contig_collections[lca], sep="\t", index_col='protein_og')
                    genome_check = genome_native_check_sequences(rhizarian_node_sequences, vertical_sequences, protein_contig)
                    print(f"Genome check was applied: {rhizarian_node_sequences}")
                    if genome_check == False:
                        results_per_acquisition = results_per_acquisition.append({"orthogroup" : og, "rhizarian_node": rhizarian_node, "rhizarian_node_sequences": ';'.join(rhizarian_node_sequences), "data_check": data_check, "contamination_check": contamination_check, "prediction_probability": True, "genome_check": genome_check}, ignore_index=True)
                        # Eliminate the predicted origins of these sequences
                        prediction_results_sequences.drop(labels=rhizarian_node_sequences, axis=0, inplace=True)
                        print(f"Genome check was false: {rhizarian_node_sequences}")
                        continue
            # Get donor phylum and domain, and the corresponding names
            if pd.isna(donor) == False:
                donor_phylum, donor_domain = get_phylum_domain_lineage(taxonomy_prokaryotes_extended, eukaryotic_supergroups, donor, taxid=True)
                donor = taxid_to_name(donor)
                donor_phylum = taxid_to_name(donor_phylum)
                donor_domain = taxid_to_name(donor_domain)
            # Report results
            prediction_results_sequences = update_pred_sequences(prediction_results_sequences, rhizarian_node_sequences, og, rhizarian_node, origin)
            results_per_acquisition = results_per_acquisition.append(
                {
                "orthogroup" : og, 
                "rhizarian_node" : rhizarian_node, 
                "rhizarian_node_sequences" : ";".join(rhizarian_node_sequences),
                "data_check": data_check,
                "contamination_check": contamination_check,
                "lca" : lca, 
                "origin": origin, 
                "donor": donor, 
                "donor_phylum": donor_phylum,
                "donor_domain": donor_domain,
                "inference": 'sequence_classification',
                "virus_proportion_og": virus_proportion_og, 
                "prediction_probability": True,
                "genome_check": genome_check
                }, 
                ignore_index=True)
            project_events(species_tree, lca, origin, None)
            origin_inferred += 1
        results_per_og.loc[i, ['origin_inferred']] = origin_inferred
        os.chdir('../') # Return to general output dir


# Append prediction results to hgt_indices_results

# Print result tables
results_per_og.to_csv("results_per_og.csv", index=False)
results_per_acquisition.to_csv("results_per_acquisition.csv", index=False)
prediction_results_sequences.to_csv("prediction_results_sequences.csv", index=True)
dataset_prediction.to_csv("dataset_prediction.csv", index=True)

# Print the global results of and to the species tree
species_tree_results = create_event_table(species_tree)
species_tree_results.to_csv("species_tree_results.csv", index=False)
species_tree.write(features=[], outfile="species_tree_results.NHX", format_root_node=True)

visualise_species_tree(species_tree, ts)



