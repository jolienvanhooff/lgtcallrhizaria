#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import re
import sys
import collections
import subprocess
from scipy.stats import gaussian_kde, mannwhitneyu, linregress, chi2_contingency
import seaborn as sns
import matplotlib.pyplot as plt
from statsmodels.graphics.mosaicplot import mosaic
from Bio import SeqIO
import gc
from lgt_rhizaria_variables import non_terminal_and_genomic_taxa
import argparse

def select_non_terminal_branches(dataframe, select=False):
    """If only to analyse ancestral and genomic taxa (regardless of origin type)"""
    if select == True:
        print(f"The number of ancestral rhizarian node with origins: {len(dataframe)}")
        dataframe = dataframe[dataframe['lca'].isin(non_terminal_and_genomic_taxa)]
        print(f"The number of sequences with origins after selecting: {len(dataframe)}")
    return dataframe


def create_histogram_dict(minimum=0, maximum=20, interval=0.2, values=[]):
    """Get value counts for a list of values, based on a range and an interval size, returns an ordered dict with as keys the lower - inclusive - border of the bin and as values the counts within the bin"""
    splits = np.arange(minimum, maximum, interval).tolist()
    hist = collections.OrderedDict()
    hist = {s:0 for s in splits}
    for val in values:
        if isinstance(val, float) or isinstance(val, int):
            for h in hist.keys():
                if h <= val < h+interval:
                    hist[h] += 1
                    break  
            # Check if the value is larger than or equal to the final bin + the interval; if so, add it to the final bin
            bin_final = list(hist.keys())[-1]
            if val >= bin_final+interval:
                hist[bin_final] += 1
    return hist


def create_dict_dataframe(df_complete, column, value_list=[]):
    """Slices a dataframe based on a value in a given column and returns them in a dictionary, with the different values in that column as keys and the dataframe slices as values"""
    df_dict = {}
    # If no value_list is passed, use all non-NA value types
    if len(value_list) == 0:
        value_list = df_complete[column].dropna().unique().tolist()
    for x in value_list:
        df_sub = df_complete[df_complete[column] == x]
        df_dict[x] = df_sub
    return df_dict


def bash_command(cmd):
    """Function to run bash commands in a shell"""
    command_run = subprocess.call(cmd, shell=True, executable='/bin/bash')#,
                                #  stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
    if command_run == 0:
        return True
    else:
        return False


def mkdir_and_cd(directory_name):
    try:
        os.mkdir(directory_name)
    except FileExistsError:
        pass
    os.chdir(directory_name)


def checktrailingslash(directory_name):
    if directory_name.endswith("/") == False:
        directory_name += "/"
    return directory_name


def get_og_members(infile):
    with open(infile, 'r') as f:
        oglist = [x.rstrip("\n") for x in f.readlines() if x != ""]
        ogdict = {l.split(": ")[0]:l.split(": ")[1].split(" ") for l in oglist}
    return ogdict


def protein_id_to_species_abbr(identifier):
    sp = ""
    if re.search("^.+_OG\d{7}", identifier) != None:
        name, og = re.findall("^(.+)_(OG\d{7})$", identifier)[0]
        sp = name[0:-6]
    elif re.search("^[\w-]+\d{6}", identifier) != None:
        sp = identifier[0:-6]
    # LamDig is a merger of LamDig_1 and LamDig_2
    if "LamDig_" in sp:
        return "LamDig"
    elif sp != "":
        return sp 
    else:
        return "unknown"
        #print(f"species abbreviation not found from sequence {identifier}")


def protein_id_without_og(identifier):
    """Returns the short ("OG-less") version of the inhouse sequence identifier"""
    if "_OG" in identifier:
        name, og = re.findall("^(.+)_(OG\d{7})$", identifier)[0]
        return name
    elif len(identifier) > 6:
        return identifier
    else:
        print(f"improper identifier: {identifier}")


def local_densities(data=[], minimum=0, maximum=100, sample_size=10, interval=0.1):
    """Creates a gaussian and searches for local densities in a sample size, returning an ordered dictionary (ordered according to the density)"""
    data = np.array(data)
    data = data[~np.isnan(data)]
    kde = gaussian_kde(data)
    samples = np.arange(minimum, maximum, sample_size)
    # for each sample, obtain the value with the highest density
    samples_max_dict = {}
    for idx, x in np.ndenumerate(samples):
        subsample = np.arange(x, x+sample_size, interval)
        probs = kde.evaluate(subsample)
        maximum_index = probs.argmax()
        maximum = subsample[maximum_index]
        samples_max_dict[round(maximum, 2)] = probs.max()
    samples_dict_ordered = collections.OrderedDict(sorted(samples_max_dict.items(), key=lambda x: x[1], reverse=True)) 
    return samples_dict_ordered

def convert_hmmscan_tblout_to_dataframe(pfam_table_path):
    """Create a pandas df from an HMMscan tabular output file"""
    hmmscan_col_names = ['target_name', 'target_accession', 'query_name', 'query_accession', 'e-value-full', 'score-full', 'bias-full', 'e-value-best-domain', 'score-best-domain', 'bias-best-domain', 'exp', 'reg', 'clu', 'ov', 'env', 'dom', 'rep', 'inc', 'description_of_target']
    # Artifically add columns to deal with the target descriptions containing spaces
    for i in range(0,21,1):
        addcol = f'col_add_{i}'
        hmmscan_col_names.append(addcol)
    pfam_hmm_df = pd.read_csv(pfam_table_path, skiprows=3, skipfooter=10, engine='python', names=hmmscan_col_names, delim_whitespace=True)
    return pfam_hmm_df

def index_in_list(a_list, index):
    check = index < len(a_list)
    return check


def plot_distributions(taxon, pattern, hue, dataframe, datatype, maximum=None):
    """For a given taxonomic level, plot the distributions for a given evolutionary pattern, with separate distributions for different kinds of the 'hue', such as origin"""
    # Drop NA values
    dataframe_non_na = dataframe.dropna(subset=[pattern], axis=0)
    # Check if the values have a distribution
    if len(set(dataframe_non_na[pattern])) > 1:
        minval = dataframe_non_na[pattern].min()
        maxval = dataframe_non_na[pattern].max()
        # If a maximum value is specified, generate a bin list
        if maximum != None and maxval > maximum:
            bins = list(range(0, maximum, int(maximum/50)))
        else:
            maximum = maxval
            bins='auto'
        if datatype == 'discrete':
            dataframe_non_na = dataframe_non_na.astype({pattern: 'int'})
            sns.histplot(x=pattern, hue=hue, data=dataframe_non_na, stat='probability', common_norm=False, multiple='dodge', discrete=True, shrink=0.8)
            plt.xticks(np.arange(minval, maxval+1, 1), fontsize=10)
        elif datatype == 'continuous':
            sns.histplot(x=pattern, hue=hue, data=dataframe_non_na, stat='probability', common_norm=False, multiple='layer', bins=bins, kde=True)
            plt.xlim([0, maximum])
        plt.ylabel('relative frequency')
        plt.title(f"{taxon}_{pattern}")
        figure_name = f"{taxon}_{pattern}.pdf"
        figure_name = figure_name.replace('/', '')
        plt.savefig(figure_name, bbox_inches='tight')
        plt.clf()
        plt.close('all')
        gc.collect()


def obtain_statistics_column(pattern, dataframe, column, val1, val2, mean=False):
    """Perform Mann-Whitney-U test on distributions for a given evolutionary pattern; compare for a given column in the dataframe two values (val1, val2); also calculate and return the median (or mean) values"""
    dataframe_1 = dataframe[dataframe[column] == val1]
    dataframe_2 = dataframe[dataframe[column] == val2]
    if mean==True:
        statistic1 = dataframe_1[pattern].mean()
        statistic2 = dataframe_2[pattern].mean()
    else:
        statistic1 = dataframe_1[pattern].median()
        statistic2 = dataframe_2[pattern].median()
    # Check if the values have a distribution
    if len(set(dataframe_1[pattern])) > 1 and len(set(dataframe_2[pattern])) > 1:
        p_value = mannwhitneyu(dataframe_1[pattern], dataframe_2[pattern]).pvalue
    else:
        p_value = "NA"
    return(statistic1, statistic2, p_value)


def scatter_plot_columns(dataframe, col1, col2, name=None, labels=True):
    """Create scatter plots of 2 types of events with (ancestral/current-day) species as data and the counts of these events as their coordinates; also perform linear regression"""
    figname = f"{name}_{col1}-vs-{col2}.pdf" if name!=None else f"{col1}-vs-{col2}.pdf"
    slope, intercept, rvalue, pvalue, stderr = linregress(dataframe[col1], dataframe[col2])
    ax = sns.regplot(x=col1, y=col2, data=dataframe, color='b', scatter_kws={'s':3}, line_kws={'label':"y={0:.1f}x+{1:.1f}\nr={2:.4f}\np={3:.4f}".format(slope,intercept,rvalue,pvalue)})
    if labels==True:
        # Place labels given by indices to the left side of the dot
        for idx in dataframe.index:
            ax.text(dataframe[col1][idx]+3, dataframe[col2][idx], idx, horizontalalignment='left', size='x-small', color='black', weight='semibold')
    ax.legend()
    plt.title(f"{col1}-vs-{col2} instances")
    plt.savefig(figname, bbox_inches='tight')
    plt.close()
    return rvalue, pvalue


def select_best_hit_per_query(df, select='evalue'):
    """Selects the best hits per query in a dataframe based on a diamond blast output"""
    new_df = pd.DataFrame(columns=df.columns)
    for q in df.qseqid.unique():
        df_query = df[df.qseqid == q]
        if select == 'evalue':
            df_query.sort_values(by=[select], ascending=True, inplace=True, ignore_index=True)
        else:
            df_query.sort_values(by=[select], ascending=False, inplace=True, ignore_index=True)
        select_row = df_query.iloc[0]
        new_df = new_df.append(select_row, ignore_index=True)
    return new_df

def plot_mosaic(taxon, pattern, hue, data, data_in_df=False):
    """For a given taxonomic level, a.k.a. name, plot a mosaic for a contingency table"""
    # Check if the values have a distribution
    if data_in_df == True:
        datadict = {}
        for x in data[pattern].unique():
            for y in data[hue].unique():
                datadict[(x, y)] = len(data[(data[pattern] == x) & (data[hue] == y)])
    else:
        datadict = data
    title = f"{taxon}_{pattern}"
    figure_name =  f"{taxon}_{pattern}.pdf"
    figure_name = figure_name.replace('/', '')
    mosaic(datadict, title=title, labelizer=lambda k: datadict[k])
    plt.savefig(figure_name, bbox_inches='tight')
    plt.close()

def get_protein_lengths_from_fasta(fasta_location):
    """From a fasta file, collects the lengths of all proteins and includes them in a list (returns this list)"""
    records = SeqIO.index(fasta_location, 'fasta')
    lengths = [len(records[rec]) for rec in records]
    del records
    return lengths

def preprocess_evolutionary_histories(dataframe, duplication_origin=False, lateral_multilevel=False, implied=True, ancestral=False, levelsep=True):
    """Selects rows from the evolutionary_histories_subtree dataframe in order to study patterns (e.g. functional, evolutionary)"""
    # Only maintain the rows with an inferred origin
    dataframe.dropna(subset=['rhizarian_node_origin'], inplace=True)
    # Remove subtrees whose parent is a duplication node: these are already included at this taxonomic level 
    if duplication_origin == False:
        dataframe = dataframe[dataframe.subtree_origin != 'duplication']
    # For rhizarian_nodes with lateral or invention origins, we may limit the analysis of this data point to the timepoint of that acquisition/invention
    if lateral_multilevel == False:
        dataframe = dataframe[(dataframe.rhizarian_node_origin == 'vertical') | (dataframe.rhizarian_node_status == True)]
    # Ignore subtrees that are inferred only (vertical origin with 'pushback')
    if implied == False:
        dataframe = dataframe[~dataframe.subtree.str.contains('O\.\d|I\.\d')]
    # Only keep the subtrees that correspond to an ancestral (non-current-day) species, or genomic species
    if ancestral == True:
        dataframe = dataframe[dataframe['subtree_lca'].isin(non_terminal_and_genomic_taxa)]
    # If a given feature or pattern is not studied separately on each taxonomic level, only keep the subtrees that reflect the ancestral rhizarian nodes
    if levelsep == False:
        dataframe = dataframe[dataframe.rhizarian_node_status == True]
    return dataframe

def add_row(dataframe, row_list):
    """Add a list of values as a row to a dataframe, provided that the length of the row equals the number of columns in the dataframe"""
    if len(dataframe.columns) == len(row_list):
        a_series = pd.Series(row_list, index = dataframe.columns)
        dataframe = dataframe.append(a_series, ignore_index=True)
    else:
        sys.exit(f"number of columns ({len(dataframe.columns)}) does not correspond to the number of values in the new row ({len(row_list)})")
    return dataframe


def get_mean_distribution_subset(dataframe, column_select, column_value, column_distribution):
    """"Create a subset of a dataframe based on a certain value in a given column, and use that to calculate a mean value of another column"""
    dataframe_sub = dataframe[dataframe[column_select] == column_value]
    mean_subset = dataframe_sub[column_distribution].mean(skipna=True)
    return mean_subset


def perform_mannwhitneyu_subsets(dataframe, column_select, column_value1, column_value2, column_distribution):
    """"Generate dataframe subsets and compare the distributions of these subsets for another column"""
    dataframe_sub1 = dataframe[dataframe[column_select] == column_value1]
    dataframe_sub2 = dataframe[dataframe[column_select] == column_value2]
    # Remove NaN values
    values1 = dataframe_sub1[column_distribution].dropna()
    values2 = dataframe_sub2[column_distribution].dropna()
    if len(set(values1)) > 1 and len(set(values2)) > 1:
        p_value = mannwhitneyu(values1, values2).pvalue
    else:
        p_value = "NA"
    return p_value

def plot_confusion_matrix(cnf_matrix, class_names, filename):
    fig, ax = plt.subplots()
    # create heatmap
    sns.heatmap(pd.DataFrame(cnf_matrix), annot=True, cmap="YlGnBu" ,fmt='g')
    ax.xaxis.set_label_position("top")
    tick_marks = [x + 0.5 for x in list(np.arange(len(class_names)))]
    plt.xticks(ticks=tick_marks, labels=class_names)
    plt.yticks(ticks=tick_marks, labels=class_names)
    plt.tight_layout()
    plt.title('Confusion matrix', y=1.1)
    plt.ylabel('Actual label')
    plt.xlabel('Predicted label')
    plt.savefig(f"{filename}.cnf_matrix.pdf", bbox_inches='tight')
    plt.clf()
    plt.close('all')

def flatten(mylist, n):
    """Flatten a list of lists (of lists), give the number of layers"""
    if n == 0:
        return mylist
    return flatten([j for i in mylist for j in i], n - 1)

def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def get_median_value_numeric_data(dataframe, protein_set, val):
    """Used to calculate median values for specific features from (sub)sets of proteins - used for intron and genomic context analysis"""
    proteins = protein_set.index.tolist()
    # Calculate the median, ignoring NA values of the mean distance
    dist_median = dataframe[val][proteins].median(skipna=True)
    return dist_median


def calculate_stats_pairwise(dataframe, protein_set1, protein_set2, val):
    """Used to obtain p-values for specific features from subsets of proteins - used for intron and genomic context analysis"""
    proteins1 = protein_set1.index.tolist()
    proteins2 = protein_set2.index.tolist()
    values1 = dataframe[val][proteins1]
    values2 = dataframe[val][proteins2]
    if all(i > 1 for i in [len(set(values1)), len(set(values2))]):
        p_value = mannwhitneyu(values1, values2, use_continuity=True).pvalue
        return p_value
    else:
        return np.nan


def get_counts_proportions_binary(dataframe, protein_set, val):
    """Used to get presences/absences of specific features for (sub)sets of proteins - used for intron analysis"""
    proteins = protein_set.index.tolist()
    predictions_subset = dataframe[dataframe.index.isin(proteins)]
    predictions_subset_non_zero = predictions_subset[predictions_subset[val] > 0]
    positive = len(predictions_subset_non_zero)
    negative = len(predictions_subset) - len(predictions_subset_non_zero)
    proportion = np.nan
    if len(predictions_subset) > 0:
        proportion = positive / len(predictions_subset)
    return positive, negative, proportion


def calculate_stats_pairwise_binary(dataframe, protein_set1, protein_set2, val):
    """Used to obtain p-values for specific binary features from subsets of proteins - used for intron analysis"""
    positive1, negative1, _ = get_counts_proportions_binary(dataframe, protein_set1, val)
    positive2, negative2, _ = get_counts_proportions_binary(dataframe, protein_set2, val)
    if all(i > 0 for i in [positive1, negative1, positive2, negative2]):
        chi2, p_value, dof, ex = chi2_contingency([[positive1, negative1], [positive2, negative2]], correction=True)
        return p_value
    else:
        return np.nan