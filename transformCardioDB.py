#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 12:15:48 2023
@author: Charlotte Bodenmüller

DESCRIPTION: Reads the genomic variation database "cardioDBwithREF.csv" and 
    transforms into a tsv file that is readable by NIRVANA. Only substitution 
    variants are considered. A corrected column for the reference allel was
    added to the database beforehand.
    
    This code was a result of a project in the UPV class 2023 "Design and 
    management of genomic information systems"
"""

import pandas as pd
import re
from natsort import natsort_keygen

# =============================================================================
# Data cleaning
# =============================================================================


# Import cardioDB
cardioDB = pd.read_csv("cardioDBwithREF.csv", delimiter = ",")


# Create the database that will be transformed into the tsv file in the end
cardioDB_prep = cardioDB.copy()


# Keeping only needed columns and renaming them according to Nirvana standard
cardioDB_prep = cardioDB[['Nucleotide.Change', 'OMGL.class', 'LMM.class', 'Phenotype', 'Location.GRCh37.', 'correct_ref']]
cardioDB_prep = cardioDB_prep.rename(columns={"OMGL.class": "OmglClass", "LMM.class" : "LmmClass","correct_ref": "REF"}, errors="raise")


# Dropping rows that do not contain a substitution
error_rows = []
for i in range(len(cardioDB_prep["Nucleotide.Change"])):
    if not re.search("^c\.[*-]*[0-9]+[+-_]*[0-9]*[ACGT]>[ACGT]+$", cardioDB_prep["Nucleotide.Change"][i]):
        error_rows.append(i)

cardioDB_prep = cardioDB_prep.drop(error_rows)


# Collecting info for CHROM and POS
cardioDB_prep[["CHROM","POS"]] = cardioDB_prep["Location.GRCh37."].str.split(":", expand=True)
cardioDB_prep = cardioDB_prep.drop(columns = "Location.GRCh37.")


# Splitting the values of column Nucleotide.Change
cardioDB_prep[["wrong_reference", "ALT"]] = cardioDB_prep["Nucleotide.Change"].str.split(">", expand=True)
cardioDB_prep = cardioDB_prep.drop(columns = ["wrong_reference", "Nucleotide.Change"])


# Dropping rows where important values are nan and replace nan values of other columns
cardioDB_prep = cardioDB_prep.dropna(subset=["CHROM", "POS", "REF", "ALT"])


# Dropping double variations
list_columns = ["OmglClass","LmmClass","Phenotype"]
cardioDB_prep = cardioDB_prep.groupby(by=["CHROM", "POS", "REF", "ALT"])[list_columns].agg(list)
cardioDB_prep = cardioDB_prep.reset_index()


# Cleaning new list columns
for col in list_columns:
    for i in range(len(cardioDB_prep[col])):
        cardioDB_prep[col][i] = [x for x in cardioDB_prep[col][i] if not pd.isna(x)] # Deleting nans
        cardioDB_prep[col][i] = list(dict.fromkeys(cardioDB_prep[col][i])) # Deleting doubles
        cardioDB_prep[col][i] = ", ".join(str(list_el) for list_el in cardioDB_prep[col][i])
        

# Sorting
cardioDB_prep = cardioDB_prep.sort_values(by = ['CHROM', 'POS'], key=natsort_keygen()) 
cardioDB_prep = cardioDB_prep.reset_index(drop=True)
        

# =============================================================================
# Write tsv file
# =============================================================================


with open("cardioDB.tsv", "w") as f:

    columns = ["CHROM","POS","REF","ALT", "OmglClass","LmmClass","Phenotype"]
    
    
    # Writing the header    
    f.write("#title=CardioDB\n")
    f.write("#assembly=GRCh37\n")
    f.write("#matchVariantsBy=allele\n")
    f.write("#" + "\t".join(str(el) for el in columns)  + "\n")
    f.write("#categories"   + "\t."*6                   + "\n")
    f.write("#descriptions" + "\t."*6                   + "\n")
    f.write("#type"         + "\tstring"*6              + "\n")
    
    
    # Writing the data
    for i in cardioDB_prep.index:
        line = []
        for col in columns:
            el = cardioDB_prep[col][i]
            if el:
                line.append(el)
            else:
                line.append(".")            
            
        f.write("\t".join(str(el) for el in line) + "\n")
    

