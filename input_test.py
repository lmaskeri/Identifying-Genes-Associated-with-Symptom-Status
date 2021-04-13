# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 22:49:03 2021

@author: lcmas
"""

import pandas as pd


full = pd.read_csv("full_presence_absence.csv")
# print(full.head())
full_new = full.transpose() #transposing the matrix so that the strains becomes the column headers
# print(full_new.head)
homologous_genes = full.columns.values[2:] #putting all homologous into a list
# print(full_new.head())

#I need to replace the header row after transposing it because the header consists of indexing 0-67 instead of the strain numbers
new_header = full_new.iloc[0] #taking the first row that contains the strain numbers 
full_new = full_new[1:] #remake the dataframe without the top row (the indexs 0-67)
full_new.columns = new_header #set the header row as the strain numbers in order to access their column information using the strain number as the key

#pulling out the strain numbers from the transposed matrix
strain_numbers = full_new.columns.values[0:]

def strain_info(matrix, strain, response):
    genes_present = []
    genes_absent = []
    pres_abs_vals = matrix[strain] #here, matrix[strain] pulls out the entire column of pres/abs values for that strain
    #Now I can go through each gene in the list of homologous genes and pull out the value from the column that has the gene name as the row 
    for gene in homologous_genes:
        if pres_abs_vals[gene] == 1:
            genes_present.append(gene)
        if pres_abs_vals[gene] == 0:
            genes_absent.append(gene)
    if response == "present":
        return genes_present   
    elif response == "absent":
        return genes_absent
    
    
#Summary statistics for the strain that the user wishes to see 
#includes: 
#what patient the strain was extracted from 
#the number of total symptom samples for that strain symptom 
#
def summary(matrix, strain):
    print("THIS WILL BE THE SUMMARY INFORMATION")
    

def compare_two_strains(matrix, strain1, strain2):
    print("THIS WILL COMPARE TWO STRAINS")
    #first grab the columns from both strains and then compare gene by gene which gene values are present in one and absent in another, or absent in both or present in both
    
    
    
def compare_two_symptoms(matrix, sympt1, sympt2):
    print("THIS WILL COMPARE TWO SYMPTOM GROUPS")
    
    
    
##User##   

inital_response = input("Are you interested in looking at a specific strain ('strain')', comparing two strains (compare two strains), or comparing two symptom groups ('compare symptom groups')?")

if inital_response == 'strain':
    print("These are the strain numbers from the dataset that can be used: \n")
    print(strain_numbers)
    strain = input("Which strain number are you interested in? Please enter the strain number into the console. ")
    summary_response = input("Would you like the summary information about this strain? Please enter 'yes' or 'no'")
    
    if summary_response == 'yes':
        summary(full_new, strain)
    
    pres_or_abs = input("Would you like the list of genes present or absent within this strain? If yes, please enter 'present' or 'absent'. Otherwise, enter 'no' ")

    if pres_or_abs != "no":
        strain_info = strain_info(full_new, int(strain), pres_or_abs)
        print("There are ", str(len(strain_info)), "homologous genes that are", pres_or_abs, "within the genome of this strain.")
        print("The genes ", pres_or_abs, "in the strain that you selected are: ")
        print(strain_info[1:]) #excluding the "Strain Symptom" value

if inital_response == 'compare two strains':
    print("These are the strain numbers from the dataset that can be used: \n")
    print(strain_numbers)
    strain1 = input("Please enter the number for the first strain. ")
    strain2 = input("Please enter the numer for the second strain that you wish to compare the first to. ")
    compare_two_strains(full_new, strain1, strain2)
    
if inital_response == 'compare symptom groups':
    print("These are the possible symptom groups: '0' (no LUTS), '2' (OAB), '3' (UTI), '4' (UUI)")
    symptom_group1 = input("Please enter the integer for the first symptom group. ")
    symptom_group2 = input("Please enter the integer for the second symptom group. ")
    compare_two_symptoms(full_new, symptom_group1, symptom_group2)


#add on the ability to compare two different strains to see which genes are present in one but absent in the other 
#allow the user to be able to do this with entire symptom groups as well -- like compare which genes are present in one symptom and not in the other accross ALL strains for that particular symptom


