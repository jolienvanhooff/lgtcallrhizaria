#!/usr/bin/env python3

import numpy as np
import pandas as pd
import os
import re
import subprocess
from lgt_rhizaria_variables import non_terminal_and_genomic_taxa

def select_non_terminal_branches(dataframe, select=False):
    """If only to analyse ancestral and genomic taxa (regardless of origin type)"""
    if select == True:
        print(f"The number of ancestral rhizarian node with origins: {len(dataframe)}")
        dataframe = dataframe[dataframe['lca'].isin(non_terminal_and_genomic_taxa)]
        print(f"The number of sequences with origins after selecting: {len(dataframe)}")
    return dataframe

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
