#!/usr/bin/env python3

# Functions for analysing the contents and distributions of orthogroups

import sys
import os
import re
import pandas as pd
import numpy as np
import subprocess
import collections
import operator

from Bio import SeqIO
from ete3 import NCBITaxa
ncbi = NCBITaxa()

from sklearn.ensemble import GradientBoostingClassifier

from gene_tree_operations import get_phylum_domain_lineage
from functions import create_dict_dataframe, protein_id_to_species_abbr


def obtain_og_contents(record_dict, taxonomy_prokaryotes, sar_species_translation):
    """For each sequence in a BioPython SeqIO index record dictionary, determine the corresponding species and its taxonomy and group"""
    seq_species = {}
    seq_taxonomy = {}
    seq_group = {}
    group_members = {"prokaryotes": [], "eukaryotes": [], "halvarians": [], "rhizarians": [], "viruses": [], "unknown": []} # dict of lists
    for seq in record_dict:
        sp = ""
        group = ""
        taxonomy = []
        if "__" in seq:
            sp = seq.split("__")[0]
            if sp in taxonomy_prokaryotes.keys():
                taxonomy = taxonomy_prokaryotes[sp]
                group = "prokaryotes"
                sp = taxonomy[-1].split('__')[1]
            else:
                try:
                    taxid = int(sp)
                    taxonomy = ncbi.get_lineage(taxid)
                    if 10239 in taxonomy:
                        group = "viruses"
                    elif 2759 in taxonomy:
                        group = "eukaryotes"
                except:
                    group = ""
        elif re.search("_OG\d{7}$", seq) != None:
            sp = protein_id_to_species_abbr(seq)
            if sp in sar_species_translation.keys():
                taxid = sar_species_translation[sp]
                try:
                    taxonomy = ncbi.get_lineage(taxid)
                    if 543769 in taxonomy:
                        group = "rhizarians"
                    else:
                        group = "halvarians"
                except:
                    group = ""
        if group == "":
            group = "unknown"
        seq_species[seq] = sp 
        seq_taxonomy[seq] = taxonomy 
        seq_group[seq] = group
        group_members[group].append(seq)
    return seq_species, seq_taxonomy, seq_group, group_members


def delete_group(record_dict, group_assignments, group_ignore, filename):
    """Remove sequences from a BioPython SeqIO index object based on their group identity"""
    sequences_to_keep = []
    sequences_to_remove = []
    ignore_present = False
    ignore_proportion = 0
    for seq in record_dict:
        group = group_assignments[seq]
        if group == group_ignore:
            sequences_to_remove.append(seq)
            ignore_present = True 
        else:
            sequences_to_keep.append(seq)
    if len(sequences_to_remove) == len(record_dict):
        sys.exit(f"found no remaining sequences after removal!")
    ignore_proportion = len(sequences_to_remove) / len(record_dict)
    conserve_ref_records = [record_dict[key] for key in sequences_to_keep]
    SeqIO.write(conserve_ref_records, filename, "fasta")
    record_dict = SeqIO.index(filename, 'fasta')
    os.remove(filename)
    return record_dict, ignore_present, ignore_proportion


def create_output_fasta(record_dict, sequence_identifiers, filename):
    sequence_records = [record_dict[s] for s in record_dict if s in sequence_identifiers]
    SeqIO.write(sequence_records, filename, "fasta")


def find_lca(sequence_ids, species_tree, species_names):
    """Finds the 'last common ancestor' of rhizarian sequences based on the species tree"""
    sp_names = []
    for name in sequence_ids:
        sp = protein_id_to_species_abbr(name)
        sp_names.append(sp)
    sp_nodes = set([species_tree&sp for sp in sp_names])
    if len(sp_nodes) > 1:
        lca = species_tree.get_common_ancestor(list(sp_nodes)).name
    else:
        lca = sp_names[0]
    return lca


def generate_fasta_subset(name, record_dict, seq_group, subject, native_group):
    query_file = f"{name}.subject.fa"
    query_records = []
    native_file = f"{name}.native.fa"
    native_records = []
    alien_file = f"{name}.alien.fa"
    alien_records = []
    for seq in record_dict:
        group = seq_group[seq]
        if group == subject:
            query_records.append(record_dict[seq])
        elif group == native_group:
            native_records.append(record_dict[seq])
        else:
            alien_records.append(record_dict[seq])
    SeqIO.write(query_records, query_file, 'fasta')
    SeqIO.write(native_records, native_file, 'fasta')
    SeqIO.write(alien_records, alien_file, 'fasta')
    return query_file, native_file, alien_file


def run_diamond(query_fasta, db_fasta, threads, evaluecutoff, diamond_path='diamond'):
    """Create a database and run diamond"""
    # Create database 
    db_name = os.path.basename(db_fasta).rstrip(".fa")
    cmd = f"{diamond_path} makedb --in {db_fasta} --db {db_name} -p {threads}" 
    subprocess.run(cmd, shell=True)
    # Run diamond
    query_name = os.path.basename(query_fasta).rstrip(".fa")
    diamond_output = f"{query_name}-vs-{db_name}.blastp.out"
    cmd = f"touch {diamond_output}"     # Create an empty file in case no results obtained from diamond blastp
    subprocess.run(cmd, shell=True)
    cmd = f"{diamond_path} blastp -q {query_fasta} -d {db_name} -p {threads} --evalue {evaluecutoff} --out {diamond_output} --ultra-sensitive --max-target-seqs 2000"
    subprocess.run(cmd, shell=True)
    return diamond_output


def perform_mcl_clustering(og, record_dict, threads, evaluecutoff, diamond_path='diamond', inflation_factor=2):
    """Perform MCL clustering to subdivide a group of sequences"""
    # First run blast
    inf_string = str(inflation_factor).replace(".", "")
    outfilename_re = f"out.{og}.mci.I{inf_string}"
    # First check if the cluster file already exists
    if os.path.exists(outfilename_re):
        print(f"cluster file exists")
        return outfilename_re
    diamond_output = f"{og}-{og}.blastp.out"
    if os.path.exists(diamond_output) == False:
        fasta_all = f"{og}.clean.fa" # clean fasta file: no viruses and unannotated sequences
        seq_all = [seq for seq in record_dict.keys()]
        create_output_fasta(record_dict, seq_all, fasta_all)
        diamond_output = run_diamond(fasta_all, fasta_all, threads, evaluecutoff, diamond_path=diamond_path)
    elif os.stat(diamond_output).st_size == 0:
        return None
    cmd = f"cut -f 1,2,11 {diamond_output} > {og}.abc"
    subprocess.run(cmd, shell=True)
    cmd = f"mcxload -abc {og}.abc --stream-mirror --stream-neg-log10 -stream-tf 'ceil(200)' -o {og}.mci -write-tab {og}.tab"
    subprocess.run(cmd, shell=True)
    cmd = f"mcl {og}.mci -I {inflation_factor} -use-tab {og}.tab -te {threads}"
    subprocess.run(cmd, shell=True)
    outfile = None
    for filename in os.listdir():
        if filename == outfilename_re:
            outfile = filename
        elif filename.endswith(('.abc', '.mci', '.tab')):
            os.remove(filename) 
    return outfile


def retrieve_subclusters_rhizaria(mcl_clusters, rhizaria_sequences, record_dict):
    """Return a dictionary with ancestral rhizarian nodes with rhizarian members"""
    cluster_counter = 0
    subclusters_selected = {}
    with open(mcl_clusters, 'r') as f:
        for line in f:
            members = line.rstrip("\n").split("\t")
            # Only select the members that are in the cleaned OG+
            members = [m for m in members if m in record_dict.keys()]
            if len(set(members).intersection(set(rhizaria_sequences))) > 0:
                cluster_counter += 1
                rhizarian_node = f"anc.{cluster_counter}"
                subclusters_selected[rhizarian_node] = members
    return subclusters_selected


def calculate_hgt_indices(sequence_ids, native_results, alien_results, remove_self=False):
    """Read diamond blastp results and use them to, for each query sequence among a list of supplemented sequences (sequence_ids) calculate a HGT index"""
    native_df = pd.read_csv(native_results, names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'], sep="\t")
    alien_df = pd.read_csv(alien_results, names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'], sep="\t")
    if remove_self == True:
        native_df = native_df[native_df.qseqid != native_df.sseqid]
        alien_df = alien_df[alien_df.qseqid != alien_df.sseqid]
    output_df = pd.DataFrame(columns=['query', 'best_native', 'best_native_hit', 'best_alien', 'best_alien_hit', 'hgt_index1', 'hgt_index2', 'hgt_index3', 'origin_predicted'])
    for query in sequence_ids:
        output_df = hgt_index(output_df, query, native_df, alien_df)
    return output_df
    

def hgt_index(output_df, query, native_df, alien_df):
    output_dict = {'query' : query}
    native_best_score = 1
    native_best_hit = ""
    alien_best_score = 1
    alien_best_hit = ""
    for i in range(2,-1,-1):
        if query in native_df.qseqid.values:
            native_hits = native_df[native_df.qseqid == query]
            try:
                native_best_score = native_hits.bitscore.sort_values(ascending=False).iloc[i]
            except IndexError:
                native_best_score = 1
            if i == 0:
                native_best_score_index = native_hits.bitscore.sort_values(ascending=False).index[0]
                native_best_hit = native_hits.loc[native_best_score_index, 'sseqid']
        if query in alien_df.qseqid.values:
            alien_hits = alien_df[alien_df.qseqid == query]
            try:
                alien_best_score = alien_hits.bitscore.sort_values(ascending=False).iloc[i]
            except IndexError:
                alien_best_score = 1
            if i == 0:
                alien_best_score_index = alien_hits.bitscore.sort_values(ascending=False).index[0]
                alien_best_hit = alien_hits.loc[alien_best_score_index, 'sseqid']
        hgt_index = alien_best_score / native_best_score
        output_dict[f'hgt_index{i+1}'] = hgt_index
    output_dict['best_native'] = native_best_score
    output_dict['best_native_hit'] = native_best_hit
    output_dict['best_alien'] = alien_best_score
    output_dict['best_alien_hit'] = alien_best_hit
    output_df = output_df.append(output_dict, ignore_index=True)   
    return output_df


def get_taxids(sequences, taxonomy_dict):
    taxids = []
    for s in sequences:
        taxids.append(taxonomy_dict[s][-1])
    return taxids


def get_lowest_level_eukaryotes(sequences, taxonomy_dict):
    taxids = get_taxids(sequences, taxonomy_dict)
    if len(set(taxids)) == 1:
        clade = taxids[0]
        rank = ncbi.get_rank([clade])[clade]
        return clade, rank
    try:
        tree = ncbi.get_topology(taxids)
        clade = tree.get_tree_root().taxid
        rank = tree.get_tree_root().rank
    except:
        clade = ""
        rank = ""
    return clade, rank


def get_lowest_level_prokaryotes(sequences, taxonomy_dict):
    clade = ""
    rank = ""
    all_taxonomies = []
    for s in sequences:
        all_taxonomies.extend(taxonomy_dict[s])
    all_taxonomies_set = set(all_taxonomies)
    hierarchy = collections.OrderedDict({'s__':'species', 'g__':'genus', 'f__':'family', 'o__':'order', 'c__':'class', 'p__':'phylum', 'd__':'domain'})
    for h, j in hierarchy.items():
        hits = [x for x in all_taxonomies_set if h in x]
        if len(hits) == 1:
            r, clade = hits[0].split('__')
            rank = hierarchy[f'{r}__']
            return clade, rank
    return clade, rank


def select_type_sequences(sequences, taxonomy_dict, clade, rank, prokaryotes=True):
    if prokaryotes == True:
        select_clade = rank[0]+"__"+clade
    else:
        select_clade = clade
    type_sequences = [s for s in sequences if select_clade in taxonomy_dict[s]]
    return type_sequences


def get_counts_taxonomic_level(sequences, taxonomy_dict, rank='p__'):
    taxon_counts = {}
    for s in sequences:
        # Find the phylum
        taxa = [t for t in taxonomy_dict[s] if rank in t]
        # In principle, there should be only one phylum
        if len(taxa) > 0:
            if taxa[0] not in taxon_counts:
                taxon_counts[taxa[0]] = 1
            else:
                taxon_counts[taxa[0]] += 1
        else:
            sys.exit(f"more than one phylum in taxonomy detected! {taxonomy_dict[s]}")
    return taxon_counts


def annotate_prokaryotes(prokaryotic_sequences, seq_taxonomy, sequence_number, interspersing_proportion):
    """Annotate a predominantly prokaryotic group of sequences"""
    clade, rank = get_lowest_level_prokaryotes(prokaryotic_sequences, seq_taxonomy)
    # If a clade, with a rank lower than 'domain' was found
    if clade != "" and rank != "domain":
        # Found a clade at <= phylum level - check if it meets the interspersing proportion criterium, taking into account all sequences in the set (not just prokaryotic)
        type_sequences = select_type_sequences(prokaryotic_sequences, seq_taxonomy, clade, rank)
        if 1 - (len(type_sequences) / sequence_number) <= interspersing_proportion:
            return clade, rank 
        # Reset; try later at phylum level with interspersing value
        else:
            clade, rank = "", ""
    # Reset; try later at phylum level with interspersing value
    else: 
        clade, rank = "", ""
    # If the lowest level was either not found, or only found at the domain level, try at phylum level by taking into account interspersing sequences
    if clade == "":
        phylum_counts = get_counts_taxonomic_level(prokaryotic_sequences, seq_taxonomy)
        # Find the phylum that occurs most often
        phylum_frequent = [l[0] for l in sorted(phylum_counts.items(), key=operator.itemgetter(1), reverse=True)][0]
        p_sequences = [s for s in prokaryotic_sequences if phylum_frequent in seq_taxonomy[s]]
        # Check if the interspersing criteria hold among these prokaryotic sequences
        if 1 - (len(p_sequences) / len(prokaryotic_sequences)) <= interspersing_proportion:
            # If they are true, try to find the lowest level among them   
            clade, rank = get_lowest_level_prokaryotes(p_sequences, seq_taxonomy)
            # Found a clade at phylum level, taking into account interspersing sequences - check if it meets the interspersing proportion criterium, taking into account all sequences in the set (not just prokaryotic)
            type_sequences = select_type_sequences(prokaryotic_sequences, seq_taxonomy, clade, rank)
            if 1 - (len(type_sequences) / sequence_number) <= interspersing_proportion:
                return clade, rank 
            else:
                clade, rank = "", ""
    # If the phylum level, with interspersing, didn't work, try at the domain level
    if clade == "":
        domain_counts = get_counts_taxonomic_level(prokaryotic_sequences, seq_taxonomy, rank='d__')
        # Find the domain that occurs most often
        domain_frequent = [l[0] for l in sorted(domain_counts.items(), key=operator.itemgetter(1), reverse=True)][0]
        d_sequences = [s for s in prokaryotic_sequences if domain_frequent in seq_taxonomy[s]]
        # Check if the interspersing criteria hold among these prokaryotic sequences
        if 1 - (len(d_sequences) / len(prokaryotic_sequences)) <= interspersing_proportion:
            # If they are true, try to find the lowest level among them   
            clade, rank = get_lowest_level_prokaryotes(d_sequences, seq_taxonomy)
            # Found a clade at domain level, taking into account interspersing sequences - check if it meets the interspersing proportion criterium, taking into account all sequences in the set (not just prokaryotic)
            type_sequences = select_type_sequences(prokaryotic_sequences, seq_taxonomy, clade, rank)
            if 1 - (len(type_sequences) / sequence_number) <= interspersing_proportion:
                return clade, rank 
            else:
                clade, rank = "", ""
    # If the clade and rank were still undefined, simply designate the group as 'prokaryotes' (paraphyletic designation)
    if clade == "":
        return "prokaryotes", ""
    else:
        sys.exit(f"no proper taxonomic annotation found for prokaryotic sequences: {prokaryotic_sequences}")


def annotate_eukaryotes(eukaryotic_sequences, seq_taxonomy, sequence_number, interspersing_proportion, eukaryotic_supergroups):
    """Annotate a predominantly eukaryotic group of sequences"""
    clade, rank = get_lowest_level_eukaryotes(eukaryotic_sequences, seq_taxonomy)
    if clade != 2759:
        # Found a clade at <= supergroup level - check if it meets the interspersing proportion criterium, taking into account all sequences in the set (not just eukaryotic)
        type_sequences = select_type_sequences(eukaryotic_sequences, seq_taxonomy, clade, rank, prokaryotes=False)
        if 1 - (len(type_sequences) / sequence_number) <= interspersing_proportion:
            return clade, rank 
        # Reset; try later at supergroup level with interspersing value
        else:
            clade, rank = "", ""
    # Reset; try later at supergroup level with interspersing value
    else:
        clade, rank = "", ""
    if clade == "":
        for taxid in eukaryotic_supergroups.keys():
            s_sequences = [s for s in eukaryotic_sequences if taxid in seq_taxonomy[s]]
            if 1 - (len(s_sequences) / len(eukaryotic_sequences)) <= interspersing_proportion:
                clade, rank = get_lowest_level_eukaryotes(s_sequences, seq_taxonomy)
                # Found a clade at supergroup or lower level, taking into account interspersing sequences - check if it meets the interspersing proportion criterium, taking into account all sequences in the set (not just eukaryotic)
                type_sequences = select_type_sequences(eukaryotic_sequences, seq_taxonomy, clade, rank, prokaryotes=False)
                if 1 - (len(type_sequences) / sequence_number) <= interspersing_proportion:
                    return clade, rank 
                else:
                    clade, rank = "", ""
    # If the clade and rank were still undefined, simply designate the clade as 'eukaryotes'
    if clade == "":
        return 2759, ncbi.get_rank([2759])[2759]
    else:
        sys.exit(f"no proper taxonomic annotation found for eukaryotic sequences: {eukaryotic_sequences}")


def split_lateral_phylum(lateral, seq_species, taxonomy_prokaryotes_extended, eukaryotic_supergroups):
    """Split the dataframe containing sequences with hgt_index > 1 into subdataframes best on their best hit annotation"""
    lateral['phylum'] = ''
    for i, j in lateral.iterrows():
        hit = j['best_alien_hit']
        species = seq_species[hit]
        phylum, _ = get_phylum_domain_lineage(taxonomy_prokaryotes_extended, eukaryotic_supergroups, species, taxid=True)
        lateral.loc[i, 'phylum'] = phylum
    subsets = create_dict_dataframe(lateral, 'phylum')
    return subsets


def annotate_alien(lateral, seq_taxonomy, interspersing_proportion, eukaryotic_supergroups):
    """Find the LCA of the species representing the best alien hits"""
    best_alien_hits = [h for h in lateral.best_alien_hit.tolist() if pd.isna(h) == False]
    if len(best_alien_hits) == 0:
        return np.nan, np.nan
    prokaryotic_hits = [h for h in best_alien_hits if 'd__Archaea' in seq_taxonomy[h]]
    prokaryotic_hits.extend([h for h in best_alien_hits if 'd__Bacteria' in seq_taxonomy[h]])
    eukaryotic_hits = [h for h in best_alien_hits if 2759 in seq_taxonomy[h]]
    hit_no = len(best_alien_hits)
    # The default values are "mix" and ""
    clade, rank = "", ""
    if 1 - (len(prokaryotic_hits) / hit_no) <= interspersing_proportion:
        # Node likely prokaryotic, so use GTDB
        clade, rank = annotate_prokaryotes(prokaryotic_hits, seq_taxonomy, hit_no, interspersing_proportion)
    elif 1 - (len(eukaryotic_hits) / hit_no) <= interspersing_proportion:
        # Node likely eukaryotic, so use NCBI
        clade, rank = annotate_eukaryotes(eukaryotic_hits, seq_taxonomy, hit_no, interspersing_proportion, eukaryotic_supergroups)
    if clade == "":
        clade, rank = "mix", ""
    return clade, rank


def get_overlaps(mylist):
    overlaps=[]
    for i in range(0,len(mylist)-1):
        for j in range(i+1, len(mylist)):
            overlap = len(set(mylist[i]).intersection(set(mylist[j])))
            overlaps.append(overlap)
    return overlaps


def merge_clusters(mylist):
    new_clusters = []
    merged = []
    for i in range(0,len(mylist)-1):
        for j in range(i+1, len(mylist)):
            overlap = len(set(mylist[i]).intersection(set(mylist[j])))
            if overlap > 0 and mylist[i] not in merged and mylist[j] not in merged:
                merge = mylist[i] + mylist[j]
                new_clusters.append(merge)
                merged.append(mylist[i])
                merged.append(mylist[j])
    for n in mylist:
        if n not in merged:
            new_clusters.append(n)
    return new_clusters


def cluster_self(sequence_ids, self_indices, self_index_threshold):
    """Create a list of lists (clusters) with sequences likely forming a monophyletic group based on their self indices"""
    best_hits = {}
    self_indices_index = self_indices.set_index('query', inplace=False)
    for s in sequence_ids:
        if s in self_indices_index.index:
            best_self_hit = self_indices_index.loc[s, 'best_alien_hit']
            self_index = self_indices_index.loc[s, 'hgt_index1']
            if self_index > self_index_threshold and pd.isna(best_self_hit) == False:
                best_hits[s] = best_self_hit
        else:
            sys.exit(f"Sequence not found in self_indices: {s}")
    clusters_initial = []
    # Create initial clusters based on chains of best hits 
    for s in best_hits.keys():
        current_cluster = [s]
        addition = best_hits[s]
        while addition not in current_cluster:
            current_cluster.append(addition)
            if addition in best_hits:
                addition = best_hits[addition]
            else:
                break
        clusters_initial.append(current_cluster)
    # Check for overlap, merge if there is overlap
    print('initial_clusters', clusters_initial)
    clusters = clusters_initial.copy()
    overlap = get_overlaps(clusters)
    while all(x == 0 for x in overlap) == False:
        clusters = merge_clusters(clusters)
        overlap = get_overlaps(clusters)
    # Add the unclustered sequences to the clusters of their best hits; alternatively, make them their own
    print('merged_clusters', clusters)
    sequences_clustered = [seq for cluster in clusters for seq in cluster]
    sequences_unclustered = [seq for seq in sequence_ids if seq not in sequences_clustered] 
    for s in sequences_unclustered:
        if s in best_hits:
            for i, c in enumerate(clusters):
                if best_hits[s] in c:
                    clusters[i].append(s)
                    continue
        else:
            clusters.append([s])
    print(clusters)
    # Clean clusters from duplicates (not likely)
    clusters = [list(set(c)) for c in clusters]
    print(clusters)
    return clusters


def genome_native_check_sequences(rhizarian_node_sequences, vertical_sequences, protein_contig):
    """Check if at least on sequence on the scaffold on which the protein is encoded has a (direct) vertical provenance"""
    sequences_same_contig = [] 
    for seq in rhizarian_node_sequences:
        contig = protein_contig.loc[seq]['scaffold']
        seq_con = protein_contig[protein_contig['scaffold'] == contig].index.tolist()
        sequences_same_contig.extend(seq_con)
    if len(set(sequences_same_contig).intersection(set(vertical_sequences))) > 0:
        return True
    else:
        return False


def classify_sequence_origins(dataset, clf, prediction_results_sequences):
    """Use a supervised machine learning classifier to predict the origins of the sequences and the probabilities"""
    sequences = dataset.index.tolist()
    features = np.array(dataset)
    predictions = clf.predict(features)
    probabilities = clf.predict_proba(features)
    classes = clf.classes_
    predictions_dict = {j:predictions[i] for i, j in enumerate(sequences)}
    probabilities_dict = {j:probabilities[i] for i, j in enumerate(sequences)} 
    for s in sequences:
        prediction_results_sequences.loc[s, 'rhizarian_node_origin_predicted_seq_level'] = predictions_dict[s]
        for i, class_name in enumerate(classes):
            class_probability = probabilities_dict[s][i]
            prediction_results_sequences.loc[s, f'rhizarian_node_origin_probability_{class_name}'] = class_probability
    return prediction_results_sequences, predictions_dict, probabilities_dict, classes


def update_pred_sequences(prediction_results_sequences, rhizarian_node_sequences, og, rhizarian_node, origin):
    """Append values to the prediction_results_sequences dataframe"""
    for s in rhizarian_node_sequences:
        prediction_results_sequences.loc[s, 'orthogroup'] = og
        prediction_results_sequences.loc[s, 'rhizarian_node'] = rhizarian_node
        prediction_results_sequences.loc[s, 'rhizarian_node_origin_predicted_node_level'] = origin 
    return prediction_results_sequences


def get_origin_mean_class_probabilities(classes, probabilities_dict, mean_probability_threshold):
    """For each class, calculate the mean probability across samples (=sequences) and return the class (=origin) that meets the threshold"""
    class_mean_probabilities = {j:np.mean([item[i] for item in list(probabilities_dict.values())]) for i, j in enumerate(classes)}
    for k, v in class_mean_probabilities.items():
        if v > mean_probability_threshold:
            return k
    return "ND"
