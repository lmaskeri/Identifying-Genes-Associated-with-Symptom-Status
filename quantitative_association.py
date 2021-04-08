
#Python script to do Quantitative Analysis of each pres/abs symptom matrix


# import numpy as np
import pandas as pd # import pandas
import csv

report = open("Quantitative_Association_Summary_Results.txt","w")
strain_genes_present = open("Genes_Present_Within_Strains.txt","w")

full = pd.read_csv("full_presence_absence.csv")
no_luts = pd.read_csv("no_LUTS_pres_abs.csv") 
oab = pd.read_csv("OAB_pres_abs.csv") 
uti = pd.read_csv("UTI_pres_abs.csv") 
uui = pd.read_csv("UUI_pres_abs.csv") #without index as the first column

homologous_genes = uui.columns.values[1:] #putting all homologous into a list


def pres_abs_genes_symptom(matrix, symptom):
    presence_percent = {} #making dictionary to hold each gene and it's presence percentage
    genes_not_in_symptom = [] #list to hold gene names that are NOT present in samples with symptom
    genes_in_symptom = [] #list to hold gene names that are present in samples with symptom
    for gene in homologous_genes:
        percent_present = sum(matrix[gene].values) / len(matrix[gene]) * 100 #adding the 1's in each gene column and then dividing by the total # of strains 
        presence_percent[gene] = percent_present #dictionary of each gene for this matrix
        if percent_present == 0.0:
            genes_not_in_symptom.append(gene)
        if percent_present != 0.0:
            genes_in_symptom.append(gene)
            
    report.write("There are "+ str(len(genes_in_symptom)) + " genes associated with patient samples with " + symptom + "\n")
    report.write("There are "+ str(len(genes_not_in_symptom)) + " genes NOT associated with patient samples with " + symptom + "\n")

    return presence_percent, genes_in_symptom, genes_not_in_symptom

full_no_luts, in_no_luts, not_in_no_luts = pres_abs_genes_symptom(no_luts, "no LUTS")
full_oab, in_oab, not_in_oab = pres_abs_genes_symptom(oab, "OAB")
full_uti, in_uti, not_in_uti = pres_abs_genes_symptom(uti, "UTI")
full_uui, in_uui, not_in_uui = pres_abs_genes_symptom(uui, "UUI")


#for each strain (row), where the column value is == 1, write the gene name into a list
def strain_genes(matrix):
    for i in range(len(matrix.values)): #for the number of rows (which would indicate the number of strains)
        strain_genes_present.write("\nThe genes present in strain " + str(matrix.values[i][0]) + " are: \n")
        strain_genes_present.write(str(matrix.columns[(matrix == 1).iloc[i]])) #for each row (using iloc), pull out the column headers where the value == 1 
        strain_genes_present.write("\n")
        
strain_genes(full)

report.close()
strain_genes_present.close()