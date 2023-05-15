#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  6 19:09:47 2023

@author: charlotte
"""
import json
import pandas as pd

f = open('VCF_3_annotated.json')
  
# returns JSON object as a dictionary
vcf3_annotated = json.load(f)


#%% Calculating number of variants

variants_num = 0

for i in range(len(vcf3_annotated["positions"])):
    variants_num  += len(vcf3_annotated["positions"][i]["variants"])

# Print result
print("The patient got " + str(variants_num) + " variants\n")

#%% Analyis of CardioDB variants 


cardioDBvariations = []
cardiodb_only_variants = 0
clinvar_cardiodb_shared_variants = 0


for i in range(len(vcf3_annotated["positions"])):
    for j in range(len(vcf3_annotated["positions"][i]["variants"])):

        keys = vcf3_annotated["positions"][i]["variants"][j].keys()
        if "CardioDB" in keys:
            cardioDBvariations.append(vcf3_annotated["positions"][i]["variants"][j])
            if "clinvar" in keys:
                clinvar_cardiodb_shared_variants += 1
            else:
                cardiodb_only_variants += 1
                
# Print results
print("ClinVar and CardioDB appear " + str(clinvar_cardiodb_shared_variants) + " time(s) in a variant simultaniously\n")
print("There are " + str(cardiodb_only_variants) + " variants that only appear in CardioDB\n")
print("The following variants from CardioDB are found in the vcf file: \n")
for i in range(len(cardioDBvariations)):
    print(cardioDBvariations[i]["CardioDB"])

#%% Looking for pathogenic variants

patho_var = []
patho_var_phenotypes = []
patho_var_genes = []

# Looping through positions and through their variatons
for i in range(len(vcf3_annotated["positions"])):    
    for j in range(len(vcf3_annotated["positions"][i]["variants"])):
        not_yet_patho = 1
        keys = vcf3_annotated["positions"][i]["variants"][j].keys()
        
        # Analysis of significance with clinvar
        if "clinvar" in vcf3_annotated["positions"][i]["variants"][j].keys():
            for k in range(len(vcf3_annotated["positions"][i]["variants"][j]["clinvar"])):
                if "pathogenic" in vcf3_annotated["positions"][i]["variants"][j]["clinvar"][k]["significance"]:
                    
                    # Collect pathogenic variant
                    if not_yet_patho:
                        patho_var.append(vcf3_annotated["positions"][i]["variants"][j])
                        not_yet_patho = 0
                        
                    # Collect all phenotypes related to pathogenic variations
                    if "phenotypes" in vcf3_annotated["positions"][i]["variants"][j]["clinvar"][k].keys():
                        patho_var_phenotypes.append(vcf3_annotated["positions"][i]["variants"][j]["clinvar"][k]["phenotypes"])
                        
                    # Collect all genes related to pathogenic variations
                    patho_var_genes.append(vcf3_annotated["positions"][i]["variants"][j]["transcripts"][0]["hgnc"])

# Cleaning data, dropping duplicates
patho_var_genes = list(dict.fromkeys(patho_var_genes))
patho_var_genes = sorted(patho_var_genes)
patho_var_phenotypes = sum(patho_var_phenotypes, [])
patho_var_phenotypes = list(dict.fromkeys(patho_var_phenotypes))
if "not provided" in patho_var_phenotypes:
    patho_var_phenotypes.remove("not provided")
    
# Create dataframe with potentially pathogenic variants, basic informations
patho_var_df = pd.DataFrame(columns = ["CHROM", "BEGIN","REF","ALT"])

for i in range(len(patho_var)):
    var = patho_var[i]
    row = pd.DataFrame([[var["chromosome"], var["begin"], var["refAllele"], var["altAllele"]]], columns = ["CHROM", "BEGIN","REF","ALT"])
    patho_var_df = pd.concat([patho_var_df, row], ignore_index=True)


# Print results
print("The patient has " + str(len(patho_var)) + " variant(s) could be related to pathogenic variations listed by clinvar\n")
print("These variations are:\n")
print(patho_var_df)
#print("Latex table: \n")
#print(patho_var_df.to_latex(index=False))
print("\nThe following phenotypes are related to these variations: \n")
print("\n".join(patho_var_phenotypes))
print("\n")
print("And the following genes: \n")
print(", ".join(x for x in patho_var_genes))

#%% Analysis of the genes which belong to the pathogenic variations

patho_var_genes_diseases = []
patho_var_genes_function = []


for el in patho_var_genes:
    for gene in vcf3_annotated["genes"]:
        if gene["name"] == el:
            
            # Checking related diseases to the genes
            if "clingenGeneValidity" in gene.keys():
                patho_var_genes_diseases.append(gene["clingenGeneValidity"][0]["disease"])
            
            # Checking gene description with omim
            if "omim" in gene.keys() and "description" in gene["omim"][0].keys():
                patho_var_genes_function.append(el + ": " + gene["omim"][0]["description"])
            else:
                patho_var_genes_function.append(el + ": no omim description")

# Print results
print("The function of these genes are described by omim as\n")
print("\n\n".join(patho_var_genes_function))
print("\n")
print("The following diseases are related to these genes: \n")
print("\n".join(patho_var_genes_diseases))
