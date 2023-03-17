#!/usr/bin/env python3

import argparse
import os
import sys
import pandas as pd

from ete3 import TreeStyle
sys.path.insert(0, 'scrips/')
from gene_tree_operations import *
from data_preparations import *
from configparser import ConfigParser
ts = TreeStyle()
sys.setrecursionlimit(10**6)

# Parse arguments
parser = argparse.ArgumentParser(description = "This script identifies the origins ('lateral', 'vertical' or 'invention') of rhizarian genes from phylogenies")
parser.add_argument("-l", metavar = "list", type=str, help ="list of OGs for which phylogenies are searched for, and, if found, analysed")
parser.add_argument("-i", metavar = "inputdir", type=str,  help = "directory in which for each listed OG a corresponding treefile should be searched")
parser.add_argument('-o', metavar = 'outputdir', type=str, default='output_detect', help='directory for output files (DEFAULT: output_detect)')
parser.add_argument('-c', metavar = 'config', type=str, default='config_files/detect_config.ini', help='path to configuration file (DEFAULT: config_files/detect_config.ini)')
args = parser.parse_args()

# Parse configuration settings from the config file in the current working directory
config_object = ConfigParser()
config_object.read(args.c)
settings = config_object["DEFAULT"]

# Prepare arguments and configurations
oglist = get_oglist(args.l)
inputdir = os.path.abspath(args.i)
outputdir = args.o
mkdir_and_cd(outputdir)

taxonomy_prokaryotes = get_full_taxonomies(settings["taxonomy_prokaryotes"]) 
taxonomy_prokaryotes_extended = get_full_taxonomies_all_levels(taxonomy_prokaryotes) 
eukaryotic_supergroups = get_supergroups_taxids(settings["supergroups_eukaryotes"]) 
sar_species_translation = get_translation_abbr_taxids(settings["taxonomy_identifiers_sar"]) 
species_names, phyla, data_types, contig_collections = get_metadata(settings["data_rhizaria"]) 
newick_rhizaria = settings["newick_rhizaria"]
root = settings["root"]
interspersing_proportion_tree = float(settings["interspersing_proportion_tree"])
interspersing_sequences_node = int(settings["interspersing_sequences_node"])
maximum_polyphyly = int(settings["maximum_polyphyly"])
support_max = int(settings["support_max"])
support_threshold = float(settings["support_threshold"])
monophyly_separation_threshold = int(settings["monophyly_separation_threshold"])
origin_decision_nested = pd.read_csv(settings["origin_decision_scheme_nested"], sep="\t")
origin_decision_root = pd.read_csv(settings["origin_decision_scheme_root"], sep="\t")
do_genome_check = True if settings["genome_check"] == 'True' else False
do_contamination_check = True if settings["contamination_check"] == 'True' else False
contamination_sequences = get_contaminants(settings["contamination_signal_diamond"])
sister_size_threshold = int(settings["sister_size_threshold"])
do_notung = True if settings["notung"] == 'True' else False
notung_path = settings["notung_path"]
branch_length_threshold = float(settings["branch_length_threshold"])

# Create or load the output dataframes
if os.path.isfile(f'results_per_og.csv') is False:
    results_per_og = pd.DataFrame(columns=["orthogroup", "treefile", "rhizarian_nodes_merged", "viruses", "acceptors"])
    results_per_acquisition = pd.DataFrame(columns=["orthogroup", "rhizarian_node", "rhizarian_node_sequences", "tree", "tree_parents", "tree_reconciled", "data_check", "contamination_check", "lca", "origin", "parent1", "parent2", "donor", "donor_phylum", "donor_domain", "inference", "median_distance_rhizarian_node", "median_distance_sister1", "median_distance_sister2", "duplication", "loss", "acceptor_proportion_node", "virus_proportion_tree"])
else:
    results_per_og = pd.read_csv("results_per_og.csv")
    results_per_acquisition = pd.read_csv("results_per_acquisition.csv")
# Load Rhizaria species tree as PhyloTree object with events as attributes
species_tree = open_species_tree(newick_rhizaria)

# Analysis per OrthoGroup (per tree)
for i, og in enumerate(oglist):
    # Print the so far collected results
    results_per_og.to_csv("results_per_og.csv", index=False)
    results_per_acquisition.to_csv("results_per_acquisition.csv", index=False)
    # Check if the current og has already been investigated - if so do not study it   
    if og in results_per_og.orthogroup.values:
        continue 
    print("\n\n"+og)
    treefile = f'{inputdir}/{og}.iqtree.treefile'
    results_per_og = results_per_og.append({"orthogroup" : og}, ignore_index=True)
    if not os.path.exists(treefile):
        print(f"no treefile found")
        continue
    results_per_og.loc[i, ['treefile']] = f'{og}.iqtree.treefile'
    tree = open_tree(treefile, support_max)
    annotate_leaves(tree, taxonomy_prokaryotes, sar_species_translation)

    # Verify that the tree has rhizarian sequences, in that case go to the next OG (no analysis)
    if len(tree.search_nodes(group="rhizarians")) == 0:
        print(f"no rhizarian sequences in tree")
        continue
    
    # Remove viruses and store the original tree
    viruses_tree, tree_original, virus_proportion_tree = prune_ignore(tree, 'viruses')
    results_per_og.loc[i, ['viruses']] = viruses_tree

    # Remove the unknown leaves
    unknown, _, _ = prune_ignore(tree, 'unknown')

    # Determine how the tree should be analysed
    rhizarian_leaves, outgroup_leaves = get_outgroup_leaves(tree, "rhizarians")
    
    # Verify if the tree might actually be rhizarian-specific, despite the presence of non-rhizarian sequences ('post-acquisition' analysis only)
    acceptors = False
    acceptor_proportion_tree = 0
    if len(outgroup_leaves) / (len(rhizarian_leaves) + len(outgroup_leaves)) < interspersing_proportion_tree:
        # Remove the putative acceptors and empty the outgroup leaves list 
        acceptors, _, acceptor_proportion_tree = prune_ignore(tree, outgroup_leaves)
        outgroup_leaves = []
    results_per_og.loc[i, ['acceptors']] = acceptors #should be 'True'

    # Start the analysis for trees with rhizarian and non-rhizarian sequences: determine the origins of the rhizarian sequences
    if len(outgroup_leaves) > 0:

        # Find the initial ancestral rhizarian nodes (note that after this, the tree is rooted on a random outgroup leaf)
        anc_nodes = get_anc_nodes(tree, outgroup_leaves, interspersing_sequences_node)

        # Merge the initial ancestral rhizarian nodes unless separated with good support
        tree, anc_nodes = find_and_attach_dispersed_relatives(tree, anc_nodes, support_threshold, monophyly_separation_threshold)
        results_per_og.loc[i, ['rhizarian_nodes_merged']] = len(anc_nodes)

        # Check if there are not too many 
        if len(anc_nodes) > maximum_polyphyly:
            print(f"too many anc rhiz nodes: {len(anc_nodes)}")
            # No reporting in results_per_acquisition
            continue

        # Analyse each (merged) ancestral rhizarian node to find its origin
        for n in anc_nodes:
            # Per ancestral rhizarian node, manipulate a copy of the tree (the working_tree), not the tree itself
            working_tree = tree.copy(method="cpickle")
            nodename = n.anc
            print(og,nodename)
            n = working_tree.search_nodes(anc=nodename)[0]

            # Reroot and annotate the tree
            working_tree.set_outgroup(n)
            annotate_nodes(working_tree, interspersing_proportion_tree, interspersing_sequences_node, support_threshold, eukaryotic_supergroups)
            reroot_tree(working_tree, n, sister_size_threshold)
            working_tree.write(format=0, outfile=f"{og}_{n.anc}_working_tree.nw")

             # Reannotate the nodes after rerooting
            annotate_nodes(working_tree, interspersing_proportion_tree, interspersing_sequences_node, support_threshold, eukaryotic_supergroups)

            # Specifically annotate the anc_rhiz_node (different rules, different input data) and replace "clade" and "rank"
            annotate_anc_node(n, species_names, phyla)

            # Collect names of the Rhizaria sequences in the node
            rhizarian_node_leaves, outgroup_leaves_cluster = get_outgroup_leaves(n, "rhizarians")
            rhizarian_node_sequences = [l.name for l in rhizarian_node_leaves]

            # Perform data or contamination check: either assess whether sequences from a single transcriptome species are contaminants ('contamination_check') or perform a regular 'data_check'
            data_check, contamination_check = verify_contamination(rhizarian_node_sequences, contamination_sequences, data_types, do_contamination_check)
            if (do_contamination_check == True and contamination_check == False) or (do_contamination_check == False and data_check == False):
                print(f"{og}, {n.anc} might contain contaminant sequences, see {rhizarian_node_sequences}")
                results_per_acquisition = results_per_acquisition.append({"orthogroup" : og, "rhizarian_node": n.anc, "rhizarian_node_sequences": ';'.join(rhizarian_node_sequences), "data_check": data_check, "contamination_check": contamination_check}, ignore_index=True)
                continue
        
            # Find the nodes that will be used to deduce the anc_rhiz_nodes origins (the parents)
            parents = find_parents(working_tree, n, support_threshold)

            # Search for the origin via the parents
            origin = np.nan
            donor = np.nan
            donor_phylum = np.nan
            donor_domain = np.nan
            clade_first_parent = np.nan
            clade_second_parent = np.nan
            decision = np.nan
            lca = np.nan
            median_distance_rhizarian_node = np.nan
            median_distance_sister1 = np.nan
            median_distance_sister2 = np.nan
            if len(parents) == 2:
                ##Identify both parents
                clade_first_parent, group_first_parent = identify_parent(parents[0], n, interspersing_proportion_tree, interspersing_sequences_node, eukaryotic_supergroups)
                clade_second_parent, group_second_parent = identify_parent(parents[1], n, interspersing_proportion_tree, interspersing_sequences_node, eukaryotic_supergroups)
                origin, donor, donor_phylum, donor_domain = detect_origin_nested(origin_decision_nested, group_first_parent, group_second_parent, clade_first_parent, clade_second_parent, taxonomy_prokaryotes_extended, eukaryotic_supergroups)
                decision = "second_sister"
                upper_parent = parents[1]
                median_distance_rhizarian_node, median_distance_sister1, median_distance_sister2 = collect_median_lengths(n, parents[0], parents[1])
            elif len(parents) == 1:
                clade_first_parent, group_first_parent = identify_parent(parents[0], n, interspersing_proportion_tree, interspersing_sequences_node, eukaryotic_supergroups)
                origin, donor, donor_phylum, donor_domain = detect_origin_root(origin_decision_root, group_first_parent, clade_first_parent, taxonomy_prokaryotes_extended, eukaryotic_supergroups)
                decision = "first_sister"
                upper_parent = parents[0]
                median_distance_rhizarian_node, median_distance_sister1, median_distance_sister2 = collect_median_lengths(n, parents[0], None)
            # No parents, no origin inferred, exit
            elif len(parents) == 0:
                decision = 'none'
                results_per_acquisition = results_per_acquisition.append({"inference": decision, "orthogroup" : og, "rhizarian_node": n.anc, "rhizarian_node_sequences": ';'.join(rhizarian_node_sequences), "data_check": data_check, "contamination_check": contamination_check}, ignore_index=True)
                print(f"no parents found for {n.anc}!")
                continue

            # If the median distance to sister1 is larger than the branch length threshold, no origin inferred, exit
            if median_distance_sister1 > branch_length_threshold:
                results_per_acquisition = results_per_acquisition.append({'median_distance_sister1': median_distance_sister1, "orthogroup" : og, "rhizarian_node": n.anc, "rhizarian_node_sequences": ';'.join(rhizarian_node_sequences), "data_check": data_check, "contamination_check": contamination_check}, ignore_index=True)
                print(f"median distances to sisters too high for {n.anc}, with median distance {median_distance_sister1}")
                continue
        
            # Convert donor and parents to taxon names, instead of taxids 
            donor = taxid_to_name(donor)
            clade_first_parent = taxid_to_name(clade_first_parent)
            clade_second_parent = taxid_to_name(clade_second_parent)

            print(f"{n.anc} - {n.clade} - {origin} - {donor} - {decision}")

            # Generate a tree visualisation
            n.add_face(TextFace(n.anc, fgcolor = 'green'), column = 0, position = 'branch-top')
            n.add_features(origin=origin, donor=donor, donor_phylum=donor_phylum, donor_domain=donor_domain)
            # (Re)label the parents; the clade now matches the clade identity as parent, so omitting the anc_rhiz_node
            if len(parents) > 0:
                for idx, p in enumerate(parents):
                    try:
                        parent_id = p.parent
                        p.add_face(TextFace(parent_id, fgcolor='red'), column=0, position='branch-top')
                        if idx == 0:
                            p.clade = clade_first_parent
                        elif idx == 1:
                            p.clade = clade_second_parent
                    except:
                        continue
            # In the tree figure, annotate all nodes with their taxonomic affiliation
            for descendant in working_tree.traverse():
                try:
                    clade = descendant.clade
                    clade = taxid_to_name(clade)
                    descendant.add_face(TextFace(clade, fgcolor = 'black'), column = 0, position = 'branch-bottom')
                except:
                    continue
            # Create layout and image
            ts.layout_fn = lambda x: None
            ts.show_leaf_name = True
            ts.show_branch_support = True
            result_visualisation = f'{og}_{n.anc}_annotated_tree.pdf'
            working_tree.render(result_visualisation, tree_style = ts)

            # Obtain the cluster with only focal sequences
            acceptors_node = False
            acceptor_proportion_node = 0
            if len(outgroup_leaves_cluster) >= 1:
                acceptors_node, _, acceptor_proportion_node = prune_ignore(n, outgroup_leaves_cluster)
            lca = find_lca(rhizarian_node_leaves, species_tree, species_names)

            # Save the pruned anc_rhiz_node (cluster) as NHX; same for the upper parent including eventual non-rhiz sequences in the rhiz-node
            n.add_features(lca=lca)
            node_nhx = f"{og}_{nodename}.NHX"
            n.write(features=[], outfile=node_nhx, format_root_node=True)
            parent_nhx = f"{og}_{nodename}.parent.NHX"
            upper_parent.write(features=[], outfile=parent_nhx, format_root_node=True)

            # Apply Notung 
            duplication = 0
            loss = 0
            speciesresult = None
            tree_reconciled = np.nan
            if (do_notung == True) and (len(rhizarian_node_leaves) >= 2):
                globalresult, speciesresult, tree_reconciled, _ = notung(node_nhx, og, n.anc, newick_rhizaria, support_threshold, notung_path=notung_path)
                duplication = globalresult['nD'][0]
                loss = globalresult['nL'][0]
                lca = globalresult[' root(G)'][0]
                
            # Report the results for this anc_rhiz_node 
            results_per_acquisition = results_per_acquisition.append(
                {
                "orthogroup" : og,
                "rhizarian_node" : n.anc,
                "rhizarian_node_sequences" : ';'.join(rhizarian_node_sequences),
                "tree" : node_nhx,
                "tree_parents": parent_nhx,
                "tree_reconciled" : tree_reconciled,
                "data_check" : data_check,
                "contamination_check" : contamination_check,
                "lca": lca,
                "origin" : origin,
                "parent1": clade_first_parent,
                "parent2": clade_second_parent,
                "donor" : donor,
                "donor_phylum": donor_phylum,
                "donor_domain": donor_domain,
                "inference" : decision,
                "median_distance_rhizarian_node" : median_distance_rhizarian_node,
                "median_distance_sister1": median_distance_sister1,
                "median_distance_sister2" : median_distance_sister2,
                "duplication": duplication,
                "loss" : loss,
                "acceptor_proportion_node" : acceptor_proportion_node, 
                "virus_proportion_tree": virus_proportion_tree
                }
                ,
                ignore_index=True)  
            
            # Append the results for this anc_rhiz_node to the species tree
            project_events(species_tree, lca, origin, speciesresult)
             
    # Analyse the rhizaria-only tree
    else:
        print("no outgroup leaves - rhizaria-only tree")
        annotate_anc_node(tree, species_names, phyla)
        anc_id = f'anc.1'

        # Collect the Rhizaria sequence names
        rhizarian_node_sequences = tree.get_leaf_names()

        # Perform data or contamination check: either assess whether sequences from a single transcriptome species are contaminants ('contamination_check') or perform a regular 'data_check'
        data_check, contamination_check = verify_contamination(rhizarian_node_sequences, contamination_sequences, data_types, do_contamination_check)
        if (do_contamination_check == True and contamination_check == False) or (do_contamination_check == False and data_check == False):
            print(f"{og} might contain contaminant sequences, see {rhizarian_node_sequences}")
            results_per_acquisition = results_per_acquisition.append({"orthogroup" : og, "rhizarian_node": "anc.1", "rhizarian_node_sequences": ';'.join(rhizarian_node_sequences), "data_check": data_check, "contamination_check": contamination_check}, ignore_index=True)
            continue

        origin = "invention"

        # Apply Notung (note that the non-rhizarian sequences have been removed already)
        duplication = 0
        loss = 0
        lca = find_lca(tree.get_leaves(), species_tree, species_names)
        speciesresult = None
        tree_rooted = None
        tree_reconciled = np.nan
        if (do_notung == True) and (len(tree.get_leaves()) >= 2):
            # Save the first version of the tree
            node_nhx_prerec = f"{og}.prerec.NHX"
            tree.write(features=[], outfile=node_nhx_prerec, format_root_node=True)
            globalresult, speciesresult, tree_reconciled, _ = notung(node_nhx_prerec, og, n.anc, newick_rhizaria, support_threshold, notung_path=notung_path)
            duplication = globalresult['nD'][0]
            loss = globalresult['nL'][0]
            lca = globalresult[' root(G)'][0]

        # Reroot the tree based on Notung's output, name the tree with an anc node identifier, and save the tree with all node features as NHX
        if tree_rooted != None:
            reroot_tree_template(tree, tree_rooted)
        node_nhx = f"{og}.NHX"
        tree.add_features(anc=anc_id, anc_status=True, origin=origin, lca=lca)
        tree.write(features=[], outfile=node_nhx, format_root_node=True)

        # Report the results for this anc_rhiz_node 
        results_per_acquisition = results_per_acquisition.append(
            {
            "orthogroup" : og,
            "rhizarian_node" : anc_id,
            "rhizarian_node_sequences" : ';'.join(rhizarian_node_sequences),
            "tree": node_nhx,
            "tree_reconciled" : tree_reconciled,
            "data_check" : data_check,
            "contamination_check": contamination_check,
            "lca": lca,
            "origin" : origin,
            "duplication": duplication, 
            "loss" : loss, 
            "acceptor_proportion_node": acceptor_proportion_tree,
            "virus_proportion_tree": virus_proportion_tree
            }
            ,
            ignore_index=True)

        # Append the results for this rhizaria-only tree to the results of the OG
        results_per_og.loc[i, ['rhizarian_nodes_merged']] = 1

        # Append the results for this rhizaria-only tree to the species tree
        project_events(species_tree, lca, origin, speciesresult) 
    


# After all trees have been inspected:
# Check sequences with lateral origin from a single species with genomes - note that this can only be done after having the complete set
# Are these laterally transferred candidates on a scaffold with vertically inherited ones?
if do_genome_check == True:
    results_per_acquisition["genome_check"] = ""
    for species, data in data_types.items():
        if data == 'genome':
            print(f'Starts checking genome of: {species}')
            protein_contig = pd.read_csv(contig_collections[species], sep="\t", index_col='protein_og')
            results_per_acquisition = genome_native_check(species, results_per_acquisition, protein_contig)

# Print the final results to tables
results_per_og.to_csv("results_per_og.csv", index=False)
results_per_acquisition.to_csv("results_per_acquisition.csv", index=False)

# Print the global results of and to the species tree
species_tree_results = create_event_table(species_tree)
species_tree_results.to_csv("species_tree_results.csv", index=False)
species_tree.write(features=[], outfile="species_tree_results.NHX", format_root_node=True)

visualise_species_tree(species_tree, ts)