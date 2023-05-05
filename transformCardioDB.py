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

# drop double variants


import pandas as pd
import re

# =============================================================================
# Data cleaning
# =============================================================================

cardioDB = pd.read_csv("cardioDBwithREF.csv", delimiter = ",")

# Dropping columns that will not be used
cardioDB = cardioDB.drop(columns=['Unnamed: 0','Gene', 'Protein.Change', 'Consequence','Type'])

# Dropping rows that do not contain a substitution
error_rows = []
for i in range(len(cardioDB["Nucleotide.Change"])):
    if not re.search("^c\.[*-]*[0-9]+[+-_]*[0-9]*[ACGT]>[ACGT]+$", cardioDB["Nucleotide.Change"][i]):
        error_rows.append(i)

cardioDB = cardioDB.drop(error_rows)

# Splitting the values of column Location.GRCh37.
cardioDB[['Location.GRCh37.chromosome','Location.GRCh37.position']] = cardioDB['Location.GRCh37.'].str.split(':', expand=True)

# Splitting the values of column Nucleotide.Change
cardioDB[['Nucleotide.Change.reference', 'Nucleotide.Change.variation']] = cardioDB['Nucleotide.Change'].str.split('>', expand=True)

# Dropping double variations
cardioDB = cardioDB.drop_duplicates(subset=['Location.GRCh37.','correct_ref', 'Nucleotide.Change.variation'])


# =============================================================================
# Write tsv file
# =============================================================================

with open("cardioDB.tsv", "w") as f:

    
    # Writing the header
    f.write("#title=CardioDB\n")
    f.write("#assembly=GRCh37\n")
    f.write("#matchVariantsBy=allele\n")
    f.write("#CHROM	POS	REF	ALT OmglClass LmmClass Phenotype\n")
    f.write("#descriptions" + "\t."*6       + "\n")
    f.write("#categories"   + "\t."*6       + "\n")
    f.write("#type\t"       + "string\t"*6  + "\n")
    
    # Writing the data
    for row_idx in cardioDB.index:
        line = [cardioDB['Location.GRCh37.chromosome'][row_idx]]
        line.append(cardioDB['Location.GRCh37.position'][row_idx]) 
        line.append(cardioDB['correct_ref'][row_idx])
        line.append(cardioDB['Nucleotide.Change.variation'][row_idx])
        line.append(cardioDB['OMGL.class'][row_idx])
        line.append(cardioDB['LMM.class'][row_idx])
        line.append(cardioDB['Phenotype'][row_idx])
        f.write("\t".join(str(el) for el in line) + "\n")
    

