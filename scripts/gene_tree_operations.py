#!/usr/bin/env python3

# Functions for performing operations on rhizarian gene trees

import sys
import re
from ete3 import PhyloTree, NCBITaxa, TextFace, PieChartFace, COLOR_SCHEMES, AttrFace, faces
import collections
import operator
import pandas as pd
import numpy as np
from io import StringIO
from functions import *
ncbi = NCBITaxa()


def open_species_tree(tree_file_path):
    """Opens the species tree file and adds relevant attributes to the internal and terminal nodes"""
    species_tree = PhyloTree(tree_file_path, format=1, sp_naming_function = None)
    for node in species_tree.traverse():
        node.add_features(acquisitions=0, vertical=0, lateral=0, invention=0, unknown=0, loss=0, duplication=0)
    return species_tree


def open_tree(tree_file_path, support_max):
    """Opens treefile and assigns the correct type of support values to nodes"""
    tree = PhyloTree(tree_file_path, sp_naming_function = None)
    if support_max != 1:
        for node in tree.traverse():
            try:
                s = node.support
                if s == 1:
                    node.support = 1.000
                else:
                    node.support = round(s / support_max, 3)
            except:
                continue
    return tree


def annotate_leaves(tree, taxonomy_prokaryotes, sar_species_translation):
    """Finds the species' taxonomy and group (prokaryotes, non-sar eukaryotes, viruses, halvaria and rhizaria) and labels the leaves with them"""
    for leaf in tree:
        name = leaf.name
        taxonomy = []
        group = ""
        group_colours = {"rhizarians": "green", "viruses": "grey", "eukaryotes": "purple", "halvarians": "blue", "prokaryotes":"red", "unknown": "black"}
        if "__" in leaf.name:
            sp = name.split("__")[0]
            if sp in taxonomy_prokaryotes.keys():
                taxonomy = taxonomy_prokaryotes[sp]
                group = "prokaryotes"
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
        elif re.search("_OG\d{7}$", name) != None:
            sp = protein_id_to_species_abbr(name)
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
        leaf.add_features(taxonomy = taxonomy)
        leaf.add_features(group = group)
        leaf.add_face(TextFace(group, fgcolor = group_colours[group]), column = 0, position = 'branch-right')


def get_taxonomies_by_leaves(leaves):
    taxids = []
    for leaf in leaves:
        taxids.append(leaf.taxonomy[-1])
    return taxids


def get_lowest_level_eukaryotes(taxids=[], leaves=[]):
    if bool(taxids) == False:
        taxids = get_taxonomies_by_leaves(leaves)
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


def get_lowest_level_prokaryotes(leaves):
    clade = ""
    rank = ""
    all_taxonomies = []
    for leaf in leaves:
        all_taxonomies.extend(leaf.taxonomy)
    all_taxonomies_set = set(all_taxonomies)
    hierarchy = collections.OrderedDict({'s__':'species', 'g__':'genus', 'f__':'family', 'o__':'order', 'c__':'class', 'p__':'phylum', 'd__':'domain'})
    for h, j in hierarchy.items():
        hits = [x for x in all_taxonomies_set if h in x]
        if len(hits) == 1:
            r, clade = hits[0].split('__')
            rank = hierarchy[f'{r}__']
            return clade, rank
    return clade, rank


def annotate_eukaryotes(leaves, total_leaf_number, interspersing_value, interspersing_proportion, child_leaves, eukaryotic_supergroups):
    """Annotate a predominantly eukaryotic clade"""
    clade, rank = get_lowest_level_eukaryotes(leaves=leaves)
    if clade != 2759:
        # Found a clade at <= phylum level - check if it makes up the majority of all leaves (not just prokaryotic) and is found in all children
        type_leaves = select_type_leaves(leaves, clade, rank, prokaryotes=False)
        representative = clade_representation(type_leaves, total_leaf_number, interspersing_value, interspersing_proportion, child_leaves)
        if representative == True:
            return clade, rank 
        # Reset; try later at supergroup level with interspersing value
        else:
            clade, rank = "", ""
    # Reset; try later at supergroup level with interspersing value
    else:
        clade, rank = "", ""
    if clade == "" and interspersing_value > 0:
        for taxid in eukaryotic_supergroups.keys():
            s_leaves = [l for l in leaves if taxid in l.taxonomy]
            if interspersing_threshold(len(s_leaves), len(leaves), interspersing_value, interspersing_proportion) == True:
                clade, rank = get_lowest_level_eukaryotes(leaves=s_leaves)
                type_leaves = select_type_leaves(leaves, clade, rank, prokaryotes=False)
                representative = clade_representation(type_leaves, total_leaf_number, interspersing_value, interspersing_proportion, child_leaves)
                if representative == True:
                    return clade, rank 
                else:
                    clade, rank = "", ""
    # If the clade and rank were still undefined, simply designate the clade as 'prokaryotes' (paraphyletic designation)
    if clade == "":
        if children_groups("non-prokaryotes", child_leaves) == True:
            return 2759, ncbi.get_rank([2759])[2759]
        else:
            return "mix", ""    
    return(clade, rank)


def get_counts_taxonomic_level(leaves, rank='p__'):
    taxon_counts = {}
    for l in leaves:
        # Find the phylum
        taxa = [t for t in l.taxonomy if rank in t]
        # In principle, there should be only one phylum
        if len(taxa) > 0:
            if taxa[0] not in taxon_counts:
                taxon_counts[taxa[0]] = 1
            else:
                taxon_counts[taxa[0]] += 1
    return taxon_counts


def select_type_leaves(leaves, clade, rank, prokaryotes=True):
    if prokaryotes == True:
        select_clade = rank[0]+"__"+clade
    else:
        select_clade = clade
    type_leaves = [l for l in leaves if select_clade in l.taxonomy]
    return type_leaves


def annotate_prokaryotes(leaves, total_leaf_number, interspersing_value, interspersing_proportion, child_leaves):
    """Annotate a predominantly prokaryotic clade"""
    clade, rank = get_lowest_level_prokaryotes(leaves)
    # If a clade, with a rank lower than 'domain' was found
    if clade != "" and rank != "domain":
        # Found a clade at <= phylum level - check if it makes up the majority of all leaves (not just prokaryotic) and is found in all children
        type_leaves = select_type_leaves(leaves, clade, rank)
        representative = clade_representation(type_leaves, total_leaf_number, interspersing_value, interspersing_proportion, child_leaves)
        if representative == True:
            return clade, rank 
        # Reset; try later at phylum level wither interspersing value
        else:
            clade, rank = "", ""
    # Reset; try later at phylum level wither interspersing value
    else: 
        clade, rank = "", ""
    # If the lowest level was either not found, or only found at the domain level, try at phylum level by taking into account interspersing sequences
    if clade == "" and interspersing_value > 0:
        phylum_counts = get_counts_taxonomic_level(leaves)
        # Find the phylum that occurs most often
        phylum_frequent = [l[0] for l in sorted(phylum_counts.items(), key=operator.itemgetter(1), reverse=True)][0]
        p_leaves = [l for l in leaves if phylum_frequent in l.taxonomy]
        # Check if the interspersing criteria hold among these prokarotic leaves
        if interspersing_threshold(len(p_leaves), len(leaves), interspersing_value, interspersing_proportion) == True:
            # If they are true, try to find the lowest level among them   
            clade, rank = get_lowest_level_prokaryotes(leaves=p_leaves)
            # Found a clade at phylum level, taking into account interspersing sequences - check if it makes up the majority of all leaves (not just prokaryotic) and is found in all children
            type_leaves = select_type_leaves(leaves, clade, rank)
            representative = clade_representation(type_leaves, total_leaf_number, interspersing_value, interspersing_proportion, child_leaves)
            if representative == True:
                return clade, rank 
            else:
                clade, rank = "", ""
    # If the phylum level, with interspersing, didn't work, try at the domain level
    if clade == "" and interspersing_value > 0:
        domain_counts = get_counts_taxonomic_level(leaves, rank='d__')
        # Find the domain that occurs most often
        domain_frequent = [l[0] for l in sorted(domain_counts.items(), key=operator.itemgetter(1), reverse=True)][0]
        d_leaves = [l for l in leaves if domain_frequent in l.taxonomy]
        # Check if the interspersing criteria hold among these prokarotic leaves
        if interspersing_threshold(len(d_leaves), len(leaves), interspersing_value, interspersing_proportion) == True:   
            clade, rank = get_lowest_level_prokaryotes(leaves=d_leaves)
            type_leaves = select_type_leaves(leaves, clade, rank)
            # Found a clade at domain level, taking into account interspersing sequences - check if it makes up the majority of all leaves (not just prokaryotic) and is found in all children
            representative = clade_representation(type_leaves, total_leaf_number, interspersing_value, interspersing_proportion, child_leaves)
            if representative == True:
                return clade, rank 
            else:
                clade, rank = "", ""
    # If the clade and rank were still undefined, simply designate the clade as 'prokaryotes' (paraphyletic designation)
    if clade == "":
        if children_groups("prokaryotes", child_leaves) == True:
            return "prokaryotes", ""
        else:
            return "mix", ""


def children_groups(focal_group, child_leaves):
    """Check if all children have members with the group annotation of the focal group; if one of them has not, return false"""
    group_check = True
    if bool(child_leaves):
        for child, current_child_leaves in child_leaves.items():
            groups = []
            for l in current_child_leaves:
                g = l.group
                if g != "prokaryotes":
                    g = "non-prokaryotes"
                groups.append(g) 
            if focal_group not in groups:
                group_check = False
    return group_check


def clade_representation(type_leaves, leaf_no, interspersing_value, interspersing_proportion, child_leaves):
    """Verify that the type leaves are found in all direct children (so on two sides of the split)"""
    representative = True
    # First check: do the type leaves meet the interspersing threshold?
    if interspersing_threshold(len(type_leaves), leaf_no, interspersing_value, interspersing_proportion) == True:
        # Second check: are they found in all daughters?
        if bool(child_leaves):
            for child, current_child_leaves in child_leaves.items():
                # What if the type leaves are not found in this daughter? The clade is not representative
                if len(set(type_leaves).intersection(set(current_child_leaves))) == 0:
                    representative = False 
                    break
    else:
        representative = False
    return representative


def interspersing_threshold(leaf_number_focus, leaf_number, interspersing_value, interspersing_proportion):
    """Checks if the interspersing criterion is met for a given set of leafs and the sum of leaves, given any type of threshold (proportional or absolute)"""
    if interspersing_proportion == True:
        if (leaf_number_focus / leaf_number) >= (1-interspersing_value):
            return True
        else:
            return False
    elif interspersing_proportion == False:
        if leaf_number_focus >= (leaf_number - interspersing_value):
            return True
        else:
            return False
    else:
        print(f"Could not assess whether the node meets the threshold, with {leaf_number_focus}, {leaf_number}, {interspersing_value}, {interspersing_proportion}")
        sys.exit() 


def annotate_node(leaves, interspersing_value, interspersing_proportion, child_leaves, eukaryotic_supergroups):
    """Returns the appropriate node annotation based on an interspersing proportion and a representation in all children"""
    prokaryotic_leaves = [l for l in leaves if 'd__Archaea' in l.taxonomy]
    prokaryotic_leaves.extend([l for l in leaves if 'd__Bacteria' in l.taxonomy])
    eukaryotic_leaves = [l for l in leaves if 2759 in l.taxonomy]
    leaf_no = len(leaves)
    # The default values are "mix" and ""
    clade, rank = "", ""
    if interspersing_threshold(len(prokaryotic_leaves), leaf_no, interspersing_value, interspersing_proportion) == True:
        # Node likely prokaryotic, so use GTDB
        clade, rank = annotate_prokaryotes(prokaryotic_leaves, leaf_no, interspersing_value, interspersing_proportion, child_leaves)
    elif interspersing_threshold(len(eukaryotic_leaves), leaf_no, interspersing_value, interspersing_proportion) == True:
        # Node likely eukaryotic, so use NCBI
        clade, rank = annotate_eukaryotes(eukaryotic_leaves, leaf_no, interspersing_value, interspersing_proportion, child_leaves, eukaryotic_supergroups)
    if clade == "":
        clade, rank = "mix", ""
    return clade, rank


def annotate_nodes(tree, interspersing_proportion_tree, interspersing_sequences_node, support_threshold, eukaryotic_supergroups):
    """Visit each internal and terminal node and define its identity on the lowest possible level and its group name - interspersing allowed between major clades - either supergroups (eukaryotes) or Bacteria/Archaea"""
    for node in tree.traverse():
        clade = ""
        rank = ""
        ## Get group labels of all leaves
        if node.support > support_threshold: 
            leaves = node.get_leaves()
            child_leaves = {child:child.get_leaves() for child in node.get_children()}
            # Only use the interspersing_sequences_node if the node has at least 3 leaves
            if 2 < len(leaves) < (1 / interspersing_proportion_tree):
                clade, rank = annotate_node(leaves=leaves, interspersing_value=interspersing_sequences_node, interspersing_proportion=False, child_leaves=child_leaves, eukaryotic_supergroups=eukaryotic_supergroups)
            else:
                clade, rank = annotate_node(leaves=leaves, interspersing_value=interspersing_proportion_tree, interspersing_proportion=True, child_leaves=child_leaves, eukaryotic_supergroups=eukaryotic_supergroups)
        node.add_features(clade = clade, rank = rank)


def prune_ignore(tree, ignore):
    """Removes uninteresting (e.g., viral) sequences from the tree, reports their presence (if found) and returns a copy of the original tree. Uses as input either a string (mostly the group name) or a list of Node or PhyloNode objects"""
    original = tree.copy(method="cpickle")
    if isinstance(ignore, str):
        ignore_leaf_names = [leaf.name for leaf in tree.iter_search_nodes(group=ignore)]
    elif isinstance(ignore, list):
        ignore_leaf_names = [leaf.name for leaf in ignore]
    remove_status = True if len(ignore_leaf_names) > 0 else False
    proportion_removed = len(ignore_leaf_names) / len(tree.get_leaf_names())
    retained_leaf_names = [x for x in tree.get_leaf_names() if x not in ignore_leaf_names]
    tree.prune(retained_leaf_names, preserve_branch_length = True)
    return remove_status, original, proportion_removed


def identify_clade(leaves):
    """Finds the clade that unites the species given by the leaves based on their taxonomies - no interspersing allowed"""
    groups = set([leaf.group for leaf in leaves])
    ##If only eukaryotes, halvarians, rhizarians and viruses, use NCBI taxonomy and return taxonomy ID
    if "prokaryotes" not in groups:
        clade, rank = get_lowest_level_eukaryotes(leaves=leaves)
    ##If prokaryotes and eukaryotes, halvarians, rhizarians or viruses, clade is undefined
    elif ("prokaryotes" in groups) and (len(groups) > 1):
        clade = "mix"
        rank = ""
    ##If prokaryotes only, use their taxonomies
    else:
        clade, rank = get_lowest_level_prokaryotes(leaves)
    return clade, rank


def get_outgroup_leaves(node, identity):
    """Gets all leaves (Node/PhyloNode objects) belonging to the ingroup and to the outgroup based on a group identity assigned to terminal nodes (leaves)"""
    ingroup_leaves = [n for n in node.search_nodes(group=identity) if n.is_leaf()]
    outgroup_leaves = list(set(node.get_leaves()).difference(set(ingroup_leaves)))
    return ingroup_leaves, outgroup_leaves


def get_subset_leaves(node, node_ignore):
    """Gets all leaves (Node/PhyloNode objects) belonging to the ingroup and to the outgroup based on a group identity assigned to terminal nodes (leaves)"""
    ignore_leaves = [n for n in node_ignore.get_leaves()]
    all_leaves = node.get_leaves()
    subset_leaves = list(set(all_leaves).difference(set(ignore_leaves)))
    return subset_leaves


def get_anc_nodes(tree, outgroup_leaves, interspersing_sequences_node):
    """Detects those nodes in the trees that most likely represent the ancestor of (a group of) rhizarian sequences in the tree, allowing for some interspersing sequences as it's child and returns them in the order of the number of rhizarian species (largest first)"""
    # Root on an outgroup leaf: find the first one that has a single sister, which itself is also an outgroup leaf
    outgroup_leaves_names = [l.name for l in outgroup_leaves]
    outgroup_leaves_names.sort()
    for o_name in outgroup_leaves_names:
        o = tree&o_name
        sisters = o.get_sisters()
        # Reroot
        tree.set_outgroup(o)
        # Make root permanent
        if len(sisters)==1 and sisters[0] in outgroup_leaves:
            break
    # Find all rhizarian monophylies
    node_candidates={}
    for x in tree.get_monophyletic(values=["rhizarians"], target_attr='group'):
        node_candidates[x]=True
    node_ancestral = []
    for x, status in node_candidates.items():
        if status == True:
            focal_node = x
            # Iterate over the parents for as long as the interspersing_sequences_node is not exceeded
            for parent in x.iter_ancestors():
                ingroup, outgroup = get_outgroup_leaves(parent, "rhizarians")
                if len(outgroup) > interspersing_sequences_node:
                    break
                else:
                    focal_node = parent
            if focal_node == x:
                node_ancestral.append(x)
            # Focal node is a candidate for redefining the ancestral rhizarian node, but only if it merges to a still availalble other ancestral rhizarian node
            else:
                merged_nodes = False
                focal_descendants = focal_node.get_descendants()
                for y in node_candidates.keys():
                    if y in focal_descendants and y != x and node_candidates[y] == True:
                        merged_nodes = True
                        node_candidates[y] = False
                if merged_nodes == True:
                    node_ancestral.append(focal_node)
                else:
                    node_ancestral.append(x)
    # Order the detected ancestral rhizarian nodes according to their number of rhizarian leafs
    node_ancestral_dict = {}
    for node in node_ancestral:
        ingroup = [n for n in node.search_nodes(group="rhizarians") if n.is_leaf()]
        node_ancestral_dict[node] = len(ingroup)
    node_ancestral = [l[0] for l in sorted(node_ancestral_dict.items(), key=operator.itemgetter(1), reverse=True)]
    # Name the ancestral rhizarian nodes using node attribute 'anc'
    if len(node_ancestral) > 0:
        for i, n in enumerate(node_ancestral):
            anc_id = f'anc.{i+1}'
            n.add_features(anc=anc_id, anc_status=True)
    return(node_ancestral)


def get_supports_merge(tree, merge_sorted_reversed, target):
    merge_supports={}
    for i, m in enumerate(merge_sorted_reversed):
        nodelist = merge_sorted_reversed[i:]
        nodelist.append(target)
        unite = tree.get_common_ancestor(nodelist)
        support = unite.support
        merge_supports[m] = support
    return merge_supports


def reorder_nodes_by_distance(tree, nodes):
    """Finds the topological tip - to -node distances for all anc nodes and orders them accordingly - largest first"""
    nodes_distance_dict = {n:tree.get_distance(n, topology_only=True) for n in nodes}
    nodes_ordered_dict = collections.OrderedDict(sorted(nodes_distance_dict.items(), key=lambda kv: kv[1], reverse=True))
    for n, d in nodes_ordered_dict.items():
        print(f"node: {n.anc}, distance: {d}")
    nodes_ordered = list(nodes_ordered_dict.keys())
    return nodes_ordered


def find_and_attach_dispersed_relatives(tree, anc_nodes, support_threshold, monophyly_separation_threshold):
    """Examines each ancestral rhizarian node to determine if it has sufficient support to be kept as a separate monophyletic ancestral rhizarian node; if not it attaches it to a nearby (larger) ancestral rhizarian node"""
    anc_nodes_check = anc_nodes.copy() ## List ordered by size, used to keep track of the nodes that are not yet examined, makes sure a node is not both a node to move and a target node
    anc_nodes = reorder_nodes_by_distance(tree, anc_nodes)
    print("anc_nodes ordered:", anc_nodes)
    nodes_to_target={}
    # Start with the furthest ancestral rhizarian node
    for n in anc_nodes:
        if n in anc_nodes_check:
            anc_nodes_check.remove(n)
            for parent in n.iter_ancestors():
                # Check all the ancestral rhizarian nodes that themselves aren't (inspected for) relocated yet
                for other_anc_node in anc_nodes_check:
                    if other_anc_node in parent.search_nodes(anc_status=True):
                        well_supported_up = 0
                        well_supported_down = 0
                        # How many well-supported nodes are in between? Check the support values for the parents of subject and target
                        for up in n.iter_ancestors():
                            if up == parent:
                                break
                            elif up.support > support_threshold:
                                well_supported_up += 1
                        for down in other_anc_node.iter_ancestors():
                            if down == parent:
                                break
                            elif down.support > support_threshold:
                                well_supported_down += 1
                        if well_supported_up + well_supported_down < monophyly_separation_threshold:
                            nodes_to_target[n]=other_anc_node
                            break
                if n in nodes_to_target.keys():
                    break
    if len(nodes_to_target) > 0:
        print("nodes merged")
        print(nodes_to_target)
    # Check nestedness: any node that is targeted to another node that is itself targeted? If so, directly target the first to the last
        nodes_to_merge = collections.defaultdict(list)
        for n, target in nodes_to_target.items():
            while target in nodes_to_target.keys():
                target = nodes_to_target[target]
            nodes_to_merge[target].append(n)
        print(nodes_to_merge)
        # Detach and attach all nodes to merge for each target node
        for target, merge in nodes_to_merge.items():
            print(f'target: {target.anc}', target)
            for m in merge:
                print(f'merge: {m.anc}', m)
            # Sort the nodes to merge according to their 'closeness' (in terms of the number of internal nodes separating it from the target) to the target
            merge_dist = {m:int(target.get_distance(m, topology_only=True)) for m in merge}
            merge_sorted = [m[0] for m in sorted(merge_dist.items(), key=operator.itemgetter(1))]
            U = target.up
            dist_fraction = target.dist / (len(merge) + 1)
            anc_id_target = target.anc
            # Check if any or more of the nodes to merge are direct child of the same parent as the target itself. If so, designate the parent of the target as the new target
            for parent in target.iter_ancestors():
                target_shift = False
                for move in merge_sorted:
                    p_relocate = move.up 
                    if p_relocate == parent:
                        merge.remove(move)
                        target_shift = True 
                if target_shift == True:
                    print('target shifted for this original anc node target:', anc_id_target)
                    target = parent 
                    U = target.up
                    dist_fraction = target.dist / (len(merge) + 1)
                else:
                    break
            # Start by adding the nodes to merge with the node (the 'move') that's furthest away from the target (in terms of nodes separating it from the target)
            merge_sorted_reversed = merge_sorted[::-1]
            # Use the support values of the new target
            merge_supports = get_supports_merge(tree, merge_sorted_reversed, target)
            for i, move in enumerate(merge_sorted_reversed):
                if move in merge:
                    support = merge_supports[move]
                    p_relocate = move.up
                    relocate = move.detach()
                    p_relocate.delete()
                    dist = (i+1) * dist_fraction
                    new_internal = f"{anc_id_target}.add.{i+1}"
                    U.add_child(name=new_internal, dist=dist, support=support)
                    N = tree&new_internal
                    N.add_child(child=relocate)
                    U = N
            target_detached = target.detach()
            U.add_child(child=target_detached, dist=dist_fraction)
            # Name the new parent according to the target node and remove all anc identifiers of the former nodes
            new_parent = tree.get_common_ancestor(merge_sorted + [target])
            #print('new parent: ', new_parent)
            merged_nodes = new_parent.search_nodes(anc_status=True)
            for f in merged_nodes:
                print(f'going to remove: {f.anc}')
                f.add_features(anc_status=False, relocated=True)
                f.del_feature("anc")
            new_parent.add_features(anc=anc_id_target, anc_status=True, merge=True)
    anc_nodes = tree.search_nodes(anc_status=True)
    return tree, anc_nodes


def check_identity_clade(node, search_identifiers):
    identity = node.clade
    check = False
    if isinstance(identity, int):
        try:
            taxonomy = ncbi.get_lineage(identity)
            if len(set(search_identifiers).intersection(set(taxonomy))) >= 1:
                check = True
        except ValueError:
            print("taxonomy not found", identity)
    return check


def farthest_leaf_sister(tree, sister_dict):
    """Find the leaf with the longest branch, based on a pruned version of the tree containing only the relevant sisters"""
    leaf_names_sisters = []
    for sis in sister_dict.keys():
        for leaf_name in sis.get_leaf_names():
            leaf_names_sisters.append(leaf_name)
    tree_copy = tree.copy(method="cpickle")
    tree_copy.prune(leaf_names_sisters, preserve_branch_length=True)
    farthest_leaf_copy = tree_copy.get_farthest_leaf()[0]
    farthest_leaf_name = farthest_leaf_copy.name
    farthest_leaf = tree&farthest_leaf_name
    return farthest_leaf


def get_root_leaf_non_relatives(tree, anc_node, sister_relatives_proportion, sister_sizes, sister_size_threshold):
    """Among leaves of all the non-SAR sister nodes, pick the best leaf to root the tree on"""
    sister_sizes_check = True
    # First check if all sisters meet the minimal number of sequences to compare their SAR proportions
    for size in sister_sizes.values():
        if size < sister_size_threshold:
            sister_sizes_check = False 
    if sister_sizes_check == True:
        # Take the sister with the fewest SAR percentage in its leaves
        root_sister = [s[0] for s in sorted(sister_relatives_proportion.items(), key=operator.itemgetter(1))][0]
        root_leaf = root_sister.get_farthest_leaf()[0]
    # If one of the sisters contains very few sequences, do not compare their SAR proportions but directly infer 
    else:
        root_leaf = farthest_leaf_sister(tree, sister_relatives_proportion)
    print(f"name of the root_leaf: {root_leaf.name}")
    return root_leaf


def reroot_tree(tree, anc_node, sister_size_threshold):
    """Rerootes the tree in order to facilitate the determination of the ancestral node's origin. The inputs are the annotated tree rooted on the anc_rhiz_node and the anc_rhiz_node itself"""
    daughters = tree.get_children()
    sister = ""
    for d in daughters:
        if d != anc_node:
            sister = d
            break
    sisters = sister.get_children()
    ##Check if one of the sisters is Stramenopile, Alveolata or SAR
    search_identifiers = [2698737, 33634, 33630]
    found_relatives = []
    sister_relatives_proportion = {}
    sister_sizes = {}
    for s in sisters:
        # First check if the sister itself is a SAR
        if check_identity_clade(s, search_identifiers) == True:
            found_relatives.append(s)
        # Then count the relative number of SAR in its leaves
        s_relatives_count = 0
        for l in s.get_leaves():
            if check_identity_clade(l, search_identifiers) == True:
                s_relatives_count += 1
        sister_relatives_proportion[s] = s_relatives_count / len(s.get_leaves())
        sister_sizes[s] = len(s.get_leaves())
    # If at least one of the sisters is SAR, remove these from the sisters on which one potentially can root
    if len(found_relatives) >= 1:
        for f in found_relatives:
            sister_relatives_proportion.pop(f)
            sister_sizes.pop(f)
    # If all sisters are SAR, root on the anc node
    if len(found_relatives) == len(sisters):
        root_leaf = anc_node
    else:
        root_leaf = get_root_leaf_non_relatives(tree, anc_node, sister_relatives_proportion, sister_sizes, sister_size_threshold)
    tree.set_outgroup(root_leaf)


def annotate_anc_node(n, species_names, subgroups):
    """Specifically annotates the Rhizaria node, not considering any interspersing sequences and using a different (non-NCBI) classification"""
    leaves = n.get_leaves()
    leaves_rhizaria = [l for l in leaves if 543769 in l.taxonomy]
    rhizaria_sequences = [l.name for l in leaves_rhizaria]
    rhizaria_abbr = []
    for name in rhizaria_sequences:
        sp = protein_id_to_species_abbr(name)
        rhizaria_abbr.append(sp)
    species = set([species_names[a] for a in rhizaria_abbr])
    genera = set([a.split(" ")[0] for a in species])
    phyla = set([subgroups[a] for a in rhizaria_abbr])
    if len(species) == 1:
        clade = list(species)[0]
        rank = "species"
    elif len(genera) == 1:
        clade = list(genera)[0]
        rank = "genus"
    elif len(phyla) == 1:
        clade = list(phyla)[0]
        rank = "phylum"
    else:
        clade = "Rhizaria"
        rank = "clade"
    n.add_features(clade=clade, rank=rank)


def verify_contamination(leaf_names, contamination_sequences, data_types, do_contamination_check):
    """Checks whether the ancestral rhizarian node contains sequences from a single rhizarian species likely to be contaminant; it returns first the data_check and then the contamination_check (only if applicable)"""
    # First assess whether the node contains sequences from a single transcriptome (with contamination check) or from at least two species (without contamination check)
    data_check = False
    if do_contamination_check == True:
        data_check = verify_data(leaf_names, data_types)
    elif do_contamination_check == False:
        data_check = verify_multispecies(leaf_names, data_types)
    if data_check == False:
        # If so, check if at least one of those is hits an outgroup sequence with unusually high identity; if this is true, it might be contamination
        if len(set(leaf_names).intersection(set(contamination_sequences.keys()))) >= 1:
            return False, False
        else:
            return False, True
    else:
        return True, "NA"


def verify_data(leaf_names, data_types):
    """Checks whether the ancestral rhizarian node meets the data requirements for further investigation into its origins: in at least two transcriptomes or in at least one genome"""
    genomes = set()
    transcriptomes = set()
    for l in leaf_names:
        if re.search("_OG\d{7}$", l) != None:
            sp = protein_id_to_species_abbr(l)
            if sp in data_types.keys():
                data_type = data_types[sp]
                if data_type == 'transcriptome':
                    transcriptomes.add(sp)
                if data_type == 'genome':
                    genomes.add(sp)
    if len(genomes) > 0:
        return True
    elif len(transcriptomes) > 1:
        return True
    else:
        return False


def verify_multispecies(leaf_names, data_types):
    """Checks whether the ancestral rhizarian node meets the data requirements for further investigation into its origins: in at least two species"""
    species_set = set()
    for l in leaf_names:
        if re.search("_OG\d{7}$", l) != None:
            sp = protein_id_to_species_abbr(l)
            if sp in data_types.keys():
                species_set.add(sp)
    if len(species_set) > 1:
        return True
    else:
        return False


def find_parents(tree, node, support_threshold):
    """Finds the two parents in the ancestors that have a good support and that are taxonomically annotated (have a specified 'clade') - not rhizarian"""
    parents = []
    count = 1
    root = tree.get_tree_root()
    while len(parents) < 2 and node != root:
        current_node = node.up
        check_rhizaria = check_identity_clade(current_node, [543769])
        if (current_node.support > support_threshold) and (current_node.clade != "") and check_rhizaria == False:
            parents.append(current_node)
            identifier = f'parent.{count}'
            current_node.add_features(parent=identifier)
            count += 1
        node = current_node
    return(parents)


def determine_sister_group(clade):
    """Finds out whether the sister node belongs to either SAR, non-SAR eukaryotes, prokaryotes or mix (undefined is mix)"""
    if isinstance(clade, int):
        taxonomy = ncbi.get_lineage(clade)
        if 2698737 in taxonomy:
            group = "SAR"
        elif 2759 in taxonomy:
            group = "eukaryotes"
        else:
            group = "mix"
    elif (clade == "prokaryotes") or (clade == "mix"):
        group = clade
    else:
        group = "prokaryotes"
    return(group)


def numeric_ranks(rank, clade):
    """Turns ranks into either one of seven fixed levels by increasing until one is found, subsequently convert to numbers (1 to 7 for species to domain). No rank found? Set it to 8"""
    rank_numbers = {'species':1, 'genus':2, 'family':3, 'order':4, 'class':5, 'phylum':6, 'domain':7}
    if rank in rank_numbers.keys():
        return rank_numbers[rank]
    elif isinstance(clade, int):
        taxonomy = ncbi.get_lineage(clade)
        for t in reversed(taxonomy):
            new_rank = ncbi.get_rank([t])[t]
            if new_rank in rank_numbers.keys():
                return rank_numbers[new_rank]
        return 8
    else:
        return 8


def validate_eukaryotes(sisters):
    """Checks if the specific clade (combined) identity of both of the sisters; if eukaryotic"""
    high_rank_taxids = [2759]
    valid_eukaryotes = True
    eukaryotic_sisters = [s for s in sisters if isinstance(s,int)]
    if len(eukaryotic_sisters) > 0:
        eukaryotic_sisters_lowest_level, _ = get_lowest_level_eukaryotes(taxids=eukaryotic_sisters)
        if eukaryotic_sisters_lowest_level in high_rank_taxids:
            valid_eukaryotes = False
    else:
        valid_eukaryotes = False
    return(valid_eukaryotes)


def detect_origin_nested(scheme, group_first_parent, group_second_parent, clade_first_parent, clade_second_parent, taxonomy_prokaryotes_extended, eukaryotic_supergroups):
    """Based on the combination of the sister clade identities and their mutual rank, this function indicates whether the rhizarian sequences are likely of vertical or lateral origin, or if it is undetermined"""
    filtered_scheme = scheme[(scheme['group_first_parent'] == group_first_parent) & (scheme['group_second_parent'] == group_second_parent)]
    if len(filtered_scheme) > 1:
        if len(filtered_scheme[filtered_scheme['clade_second_parent'] == clade_second_parent]) == 1:
            filtered_scheme = filtered_scheme[filtered_scheme['clade_second_parent'] == clade_second_parent]
        elif len(filtered_scheme[filtered_scheme['clade_first_parent'] == clade_first_parent]) == 1:
            filtered_scheme = filtered_scheme[filtered_scheme['clade_first_parent'] == clade_first_parent]
        elif len(filtered_scheme[filtered_scheme['clade_second_parent'] == 0]) == 1:
            filtered_scheme = filtered_scheme[filtered_scheme['clade_second_parent'] == 0]
        elif len(filtered_scheme[filtered_scheme['clade_first_parent'] == 0]) == 1:
            filtered_scheme = filtered_scheme[filtered_scheme['clade_first_parent'] == 0]
    if len(filtered_scheme) == 1:
        origin = filtered_scheme.iloc[0]['origin']
        donor_inference = filtered_scheme.iloc[0]['donor']
        if donor_inference == 'clade_first_parent':
            donor = clade_first_parent
        elif donor_inference == 'clade_second_parent':
            donor = clade_second_parent
        else:
            donor = donor_inference
        donor_phylum, donor_domain = get_phylum_domain_lineage(taxonomy_prokaryotes_extended, eukaryotic_supergroups, donor, taxid=True)
    else:
        sys.exit(f"current combination not found! {group_first_parent}, {group_second_parent}, {clade_second_parent}")
    return(origin, donor, donor_phylum, donor_domain)


def detect_origin_root(scheme, group_first_parent, clade_first_parent, taxonomy_prokaryotes_extended, eukaryotic_supergroups):
    """Based on the sister clade identity and its rank, this function indicates whether the rhizarian sequences are likely of vertical or lateral origin, or if it is undetermined"""
    filtered_scheme = scheme[scheme['group_first_parent'] == group_first_parent]
    if len(filtered_scheme) > 1:
        if len(filtered_scheme[filtered_scheme['clade_first_parent'] == clade_first_parent]) == 1:
            filtered_scheme = filtered_scheme[filtered_scheme['clade_first_parent'] == clade_first_parent]
        else:
            filtered_scheme = filtered_scheme[filtered_scheme['clade_first_parent'] == 0]
    if len(filtered_scheme) == 1:
        origin = filtered_scheme.iloc[0]['origin']
        donor_inference = filtered_scheme.iloc[0]['donor']
        if donor_inference == 'clade_first_parent':
            donor = clade_first_parent
        else:
            donor = donor_inference
        donor_phylum, donor_domain = get_phylum_domain_lineage(taxonomy_prokaryotes_extended, eukaryotic_supergroups, donor, taxid=True)
    else:
        sys.exit(f"current combination not found! {group_first_parent}, {clade_first_parent}")
    return(origin, donor, donor_phylum, donor_domain)


def notung(clusterfile, og, nodename, newick_rhizaria, support_threshold, rooting=False):
    """First generates Notung input, then applies Notung rearrange to the rhizarian (sub)tree"""
    notung_output = f"{og}_{nodename}_notung"
    globalresult=None
    speciesresult=None
    tree_rooted=""
    if rooting == True:
        notung_root = f"{og}_{nodename}_notung_root"
        cmd = f"java -jar /home/jolien/tools/Notung-2.9.1.5/Notung-2.9.1.5.jar -g {clusterfile} -s {newick_rhizaria} --absfilenames --root --infertransfers false --speciestag prefix --treeoutput newick --outputdir {notung_root}"
        status = bash_command(cmd)
        if status:
            clusterfile=f"{notung_root}/{clusterfile}.rooting.0"
            tree_rooted = clusterfile
        else:
            print(f"{notung_root}: notung rooting failed")
            return globalresult, speciesresult
    cmd = f"java -jar /home/jolien/tools/Notung-2.9.1.5/Notung-2.9.1.5.jar -g {clusterfile} -s {newick_rhizaria} --absfilenames --rearrange --infertransfers false --speciestag prefix --threshold {support_threshold} --edgeweights name --treeoutput nhx --outputdir {notung_output} --events --parsable --homologtablecsv"
    status = bash_command(cmd)
    if status:
        resultfile = f'{notung_output}/{clusterfile}.rearrange.0.parsable.txt'
        tree_reconciled = f'{notung_output}/{clusterfile}.rearrange.0'
        with open(resultfile, 'r') as file_object:
            lines = file_object.readlines()
            globallines = [lines[1], lines[0]]
            globalresult = pd.read_csv(StringIO("".join(globallines)), sep="\t")
            specieslines = [l for l in lines if "#S\t" in l]
            speciesresult = pd.read_csv(StringIO("".join(specieslines)), sep="\t", index_col=1)
    else:
        print(f"{notung_output}: notung failed")
    return globalresult, speciesresult, tree_reconciled, tree_rooted


def find_lca(rhizarian_leaves_cluster, species_tree, species_names):
    """Finds the 'last common ancestor' of the rhizarian sequences in a gene tree based on the species tree"""
    leaf_names = [l.name for l in rhizarian_leaves_cluster]
    sp_names = []
    for name in leaf_names:
        sp = protein_id_to_species_abbr(name)
        sp_names.append(sp)
    sp_nodes = set([species_tree&sp for sp in sp_names])
    if len(sp_nodes) > 1:
        lca = species_tree.get_common_ancestor(list(sp_nodes)).name
    else:
        lca = sp_names[0]
    return lca
    

def identify_parent(parent, node, interspersing_proportion_tree, interspersing_sequences_node, eukaryotic_supergroups):
    """Finds the specific identity ('clade') and group of a parent after removing one of its descendant nodes; after that applies node-based annotation"""
    subset_leaves = get_subset_leaves(parent, node)
    subset_leaves_names = [l.name for l in subset_leaves]
    parent_copy = parent.copy(method="cpickle")
    parent_copy.prune([parent_copy&n for n in subset_leaves_names])
    annotate_nodes(parent_copy, interspersing_proportion_tree, interspersing_sequences_node, 0, eukaryotic_supergroups)
    parent_identity = parent_copy.clade
    parent_group = determine_sister_group(parent_identity)
    return parent_identity, parent_group


def project_events(species_tree, lca, origin, speciesresult):
    """Projects the events collected for an ancestral rhizarian node or a rhizaria-only tree onto the species tree by summing the events"""
    lca_node = species_tree&lca
    lca_node.acquisitions += 1
    if origin == "vertical":
        lca_node.vertical += 1
    elif origin == "lateral":
        lca_node.lateral += 1
    elif origin == "invention":    
        lca_node.invention += 1
    else:
        lca_node.unknown += 1
    if isinstance(speciesresult, pd.DataFrame):
        for i, j in speciesresult.iterrows():
            taxon_node = species_tree&i
            num_dup = j['Dups']
            num_loss = j['Losses']
            taxon_node.duplication += num_dup
            taxon_node.loss += num_loss


def create_event_table(species_tree):
    """Summarises the results across the species tree in a dataframe"""
    species_tree_results = pd.DataFrame(columns=["node", "acquisitions", "vertical", "lateral", "invention", "unknown", "loss", "duplication"])
    for node in species_tree.traverse(strategy='preorder'):
        species_tree_results = species_tree_results.append({
            "node" : node.name, 
            "acquisitions" : node.acquisitions, 
            "vertical" : node.vertical, 
            "lateral" : node.lateral, 
            "invention" : node.invention, 
            "unknown" : node.unknown, 
            "loss" : node.loss, 
            "duplication" : node.duplication
            }
            , 
            ignore_index=True)
    return(species_tree_results)


def my_layout(node):
    if node.is_leaf():
         # If terminal node, draws its name
         name_face = AttrFace("name")
    else:
         # If internal node, draws label with smaller font size
         name_face = AttrFace("name", fsize=10)
    # Adds the name face to the image at the preferred position
    faces.add_face_to_node(name_face, node, column=0, position="branch-right")

def visualise_species_tree(species_tree, ts, outfile="species_tree_results.pdf"):
    """Modifies TreeStyle and Face objects and a pdf visualisation of the events projected onto the species node"""
    ts.show_leaf_name = False
    ts.show_branch_support = False
    ts.layout_fn = my_layout
    ts.margin_left=1
    ts.margin_right=1
    ts.margin_bottom=1
    ts.margin_top=1
    for node in species_tree.traverse(strategy="preorder"):
        node.add_face(TextFace(str(node.acquisitions), fgcolor='blue'), column = 0, position = "branch-top")
        total_acq = node.vertical + node.lateral + node.invention + node.unknown
        if total_acq > 0:
            percents = [(c/total_acq) * 100 for c in [node.vertical, node.lateral, node.invention, node.unknown]]
            node.add_face(PieChartFace(percents, 15, 15, colors=COLOR_SCHEMES['accent'][0:4]), column = 1, position = 'branch-top')
        node.add_face(TextFace(str(node.duplication), fgcolor='green'), column = 0, position = "branch-bottom")
        node.add_face(TextFace(str(node.loss), fgcolor='red'), column = 1, position = "branch-bottom")
    species_tree.render(outfile, tree_style = ts)

def reroot_tree_template(tree, rooted_tree):
    template = PhyloTree(rooted_tree, format=0)
    daughter1 = template.get_children()[0]
    daughter1_leaf_names = daughter1.get_leaf_names()
    root_leaf_nodes=[]
    for leaf_name in daughter1_leaf_names:
        try:
            leaf_node = tree&leaf_name
            root_leaf_nodes.append(leaf_node)
        except:
            print(f"no corresponding node found for {leaf_name}")
    try:
        root_node = tree.get_common_ancestor(root_leaf_nodes)
        tree.set_outgroup(root_node)
    except:
        print(f"no common ancestor found for {root_leaf_nodes}")

def genome_native_check(species, results_per_acquisition, protein_contig):
    results_to_check = results_per_acquisition[(results_per_acquisition['lca'] == species) & (results_per_acquisition['origin'] == 'lateral')]
    results_for_verification = pd.concat([results_per_acquisition, results_to_check]).drop_duplicates(keep=False)
    sequences_verification = []
    # Use the sequences of ancestral rhizarian nodes with an inferred origin
    for i, j in results_for_verification.iterrows():
        rhizarian_node_sequences = j['rhizarian_node_sequences']
        origin = j['origin']
        if (pd.isna(rhizarian_node_sequences) == False) and (pd.isna(origin) == False):
            sequences_verification.extend(rhizarian_node_sequences.split(';'))
    sequences_verification_set = set(sequences_verification)
    # Visit each result row that contains an origins data point from the current species
    for i, j in results_to_check.iterrows():
        og = j['orthogroup']
        rhiz_node = j['rhizarian_node']
        sequences = j['rhizarian_node_sequences'].split(';')
        print(f"\nChecking genomic locations for: {og}_{rhiz_node}, with sequences {sequences}")
        # Check if any of the sequences is on a scaffold that is vertically inherited
        sequences_same_contig = [] 
        for seq in sequences:
            contig = protein_contig.loc[seq]['scaffold']
            print(f"The sequence {seq} is on {contig}")
            seq_con = protein_contig[protein_contig['scaffold'] == contig].index.tolist()
            sequences_same_contig.extend(seq_con)
        # Remove the focal sequences themselves for the sequence set ont the same scaffold
        sequences_same_contig_set = set(sequences_same_contig).difference(set(sequences))
        if len(sequences_verification_set & sequences_same_contig_set) == 0:
            values_dict = results_per_acquisition.iloc[i].to_dict()
            values_dict['lca'] = np.nan
            values_dict['origin'] = np.nan
            values_dict['parent1'] = np.nan
            values_dict['parent2'] = np.nan
            values_dict['donor'] = np.nan
            values_dict['donor_phylum'] = np.nan
            values_dict['donor_domain'] = np.nan
            values_dict['inference'] = np.nan
            values_dict['median_distance_rhizarian_node'] = np.nan
            values_dict['median_distance_sister1'] = np.nan
            values_dict['median_distance_sister2'] = np.nan
            values_dict['duplication'] = np.nan
            values_dict['loss'] = np.nan
            values_dict['acceptors_node'] = np.nan
            values_dict['genome_check'] = False
            results_per_acquisition.iloc[i] = pd.Series(values_dict)
            print(f"removed origin after scaffold verification ('genome_check' is False): {og}, {rhiz_node}, with sequences {sequences}")
        else:
            print(f"Sequences on same scaffold with verification: ", sequences_verification_set & sequences_same_contig_set)
            values_dict = results_per_acquisition.iloc[i].to_dict()
            values_dict['genome_check'] = True
            results_per_acquisition.iloc[i] = pd.Series(values_dict)
    return results_per_acquisition


def get_median_branch_length(node1members, node2members):
    """"Calculate all distances from a given set of nodes to another set of nodes; then calculate their median"""
    branch_lengths = []
    for a in node1members:
        for b in node2members:
            branch_lengths.append(a.get_distance(b))
    median_branch_length = np.median(branch_lengths)
    return median_branch_length


def collect_median_lengths(node, parent1, parent2):
    """Collects medians of branch lenghts of the rhizarian node leaves to the node itself, and to parent leaves"""
    node_leaves = node.get_leaves()
    # Get median_distance_rhizarian_node
    median_distance_rhizarian_node = get_median_branch_length(node_leaves, [node])
    # Get median_distance_sister1
    sister1_leaves = get_subset_leaves(parent1, node)
    median_distance_sister1 = get_median_branch_length(node_leaves, sister1_leaves)
    if parent2 == None:
        return median_distance_rhizarian_node, median_distance_sister1, np.nan
    else:
        # Get median_distance_sister2; only take those leaves under parent2 that are not under parent1 (as opposed to parent usage in origin identification)
        sister2_leaves = get_subset_leaves(parent2, parent1)
        median_distance_sister2 = get_median_branch_length(node_leaves, sister2_leaves)
        return median_distance_rhizarian_node, median_distance_sister1, median_distance_sister2


def taxid_to_name(lineage):
    """Translate lineage/clade from taxid to name - only if it appears a taxid"""
    if isinstance(lineage, int):
        try:
            lineage = ncbi.get_taxid_translator([lineage])[lineage]
        except KeyError:
            pass
    return lineage


def get_phylum_domain_lineage(taxonomy_prokaryotes_extended, eukaryotic_supergroups, lineage, taxid=False):
    """Uses taxonomies of a lineage to find its phylum (prokaryotes), domain (eukaryotes; variable is 'phylum') and its domain (Bacteria, Archaea or Eukaryota)"""
    phylum = ''
    domain = ''
    # if the lineage itself is NaN, also return NaN for phylu and domain 
    if pd.isna(lineage):
        return np.nan, np.nan
    # Prokaryotic lineage
    if lineage in taxonomy_prokaryotes_extended.keys():
        taxonomy = taxonomy_prokaryotes_extended[lineage]
        for p in taxonomy:
            if 'p__' in p:
                phylum = p.split("__")[1]
            if 'd__'in p:
                domain = p.split("__")[1]
    # Eukaryotic / viral lineage
    else:
        # Check if the (eukaryotic/viral) lineage is given with a taxid; if not, use the name translator
        if taxid == False:
            try:
                taxonomy = ncbi.get_lineage(ncbi.get_name_translator([lineage])[lineage][0])
                for p in taxonomy:
                    if p in eukaryotic_supergroups.keys():
                        phylum = eukaryotic_supergroups[p]
                        domain = 'Eukaryota'
                        break
            except KeyError:
                phylum = lineage
                print(f'taxid not found: {lineage}')
            except ValueError:
                phylum = lineage
                print(f'lineage not found: {lineage}')
        elif taxid == True:
            try:
                taxonomy = ncbi.get_lineage(lineage)
                for p in taxonomy:
                    if p in eukaryotic_supergroups.keys():
                        phylum = eukaryotic_supergroups[p]
                        domain = 'Eukaryota'
                        break
            except KeyError:
                phylum = lineage
                print(f'taxid not found: {lineage}')
            except ValueError:
                phylum = lineage
                print(f'taxid not found: {lineage}')
    phylum = phylum if phylum != '' else lineage
    domain = domain if domain != '' else lineage
    return phylum, domain


def get_phylum_donor_lineage_dataframe(taxonomy_prokaryotes_extended, eukaryotic_supergroups, result_dataframe, taxid=False):
    """Adds a column to a result dataframe converting the donor to a phylum (prokaryotes) or supergroup (eukaryotes), if possible; also indicate if for eukaryotic donors the name (taxid=False) or the taxids should be used"""
    result_dataframe['donor_domain'] = np.nan
    result_dataframe['donor_phylum'] = np.nan
    for i, j in result_dataframe.iterrows():
        donor = j['donor']
        if pd.isna(donor) == False:
            phylum, domain = get_phylum_domain_lineage(taxonomy_prokaryotes_extended, eukaryotic_supergroups, donor, taxid=taxid)
            result_dataframe.loc[i, 'donor_phylum'] = phylum
            result_dataframe.loc[i, 'donor_domain'] = domain
    return result_dataframe


def confer_all_annotations(node1, node2):
    """Collect all annotations from a given node and place them onto another node"""
    annotations = node1.features
    for a in annotations:
        ldic = locals()
        exec(f"value = node1.{a}", globals(), ldic)
        value = ldic["value"]
        exec("node2.add_features({} = '{}')".format(a, value))