
#Python script to generate a presence absence matrix for each strain and the homologous genes that they either have or do not have 

from Bio import SeqIO
import csv

acc_numbers = open("./ecoliAccessions.txt").read().split("\n") #pulling out accession numbers for each strain

strain_symptoms = []
with open("./strain_symptom_samples.csv", "r") as in_file: #grabbing the strain # and symptom from the csv file
    reader = csv.reader(in_file, delimiter = ',')
    i = False
    for row in in_file:
        if i != False: #avoiding the first "Strain # and Participant Symtom" heading
            strain_symptoms.append(row.strip("\n").split(",")) #gives ["strain_number", "Symtom status"] 
        i = True
        

# Getting a list of only the homologous gene names
take_out = [ "[Proteobacteria]", "[Escherichia]", "[Escherichia coli]", "[Enterobacterales]","[Enterobacteriaceae]", "[Bacteria]", "[Gammaproteobacteria]"] #Names for E.coli classes that we want to take out of the gene name
homologous_genes = set() #set for homologous gene names (needs to be set because there are multiple "hypothetical proteins")
with open("./nr.fasta") as centroids: #open centroid file
    for record in SeqIO.parse(centroids, "fasta"):
        gene_name = record.description[15:] #pull out description (not using id because we don't want the accession number)
        if "MULTISPECIES:" in gene_name: #taking out multispecies
            gene_name = gene_name[14:]
        for out in take_out: #taking out the E.coli class
            if out in gene_name:
                end = len(out) + 1
                gene_name = gene_name[:-end]
        homologous_genes.add(gene_name) #appending the final homologous gene name to a list

#making a presence absence dictionary for each strain according to if they have the homologous gene present (1) or absent (0)
pres_abs = [] #list to hold each strains individual list for presence
for acc in acc_numbers:
    filename = acc + "_protein.faa" #creating fasta file name in folder
    strain_pres_abs = {} #initalizing a dictionary to hold the presence absence values for every homologus gene for each individual strain
    for hom_gene in homologous_genes:
        strain_pres_abs[hom_gene] = 0 #placing the homologous gene name as the key and the value as 0 for now
    with open("./ProteinSeqs/" + filename) as strain_file: #opening the protein file for the strain
        for record in SeqIO.parse(strain_file, "fasta"): #for each record
            strain_gene = record.description[15:] #pull out description (not using id because we don't want the accession number)
            if "MULTISPECIES:" in gene_name: #taking out multispecies
                strain_gene = strain_gene[14:]
            for out in take_out: #taking out the E.coli class
                if out in strain_gene:
                    end = len(out) + 1
                    strain_gene = strain_gene[:-end] #this gives us the final gene name to check homologus genes with
            if strain_gene in strain_pres_abs.keys(): #checking if that particular gene is in the homologous genes list
                strain_pres_abs[strain_gene] = 1 #if it is in the homolgous genes keys make the value 1 for present
                #if it is not in the homolgous genes keys the value will stay 0 for absent
        pres_abs.append(strain_pres_abs) #appending the dictionary of presence absence values for that strain to the main presence absence list
    

'''
Symptom Key:
0 = No Lower Urinary Tract Symptoms (no LUTS)

1 = Overactive bladder (OAB)

2 = Urinary Tract Infection (UTI)

3 = Urgency Urinary Incontinence (UUI)
'''
#opening the 5 csv files to write out to
no_luts =  open("no_LUTS_pres_abs.csv", "w", newline = "")
oab = open("OAB_pres_abs.csv", "w", newline = "") 
uti = open("UTI_pres_abs.csv", "w", newline = "") 
uui = open("UUI_pres_abs.csv", "w", newline = "")
out_file =   open("full_presence_absence.csv", "w", newline = "")  

#making csv writers for each to write to the file
writer1 = csv.writer(no_luts)
writer2 = csv.writer(oab)
writer3 = csv.writer(uti)
writer4 = csv.writer(uui)
writer5 = csv.writer(out_file)

#making the header with all of the homologous gene names
header = ["Strain Number"]
header.append("Strain Symptom")
for gene in homologous_genes:
    header.append(gene)

#writing the header row to each csv file
writer1.writerow(header)
writer2.writerow(header)
writer3.writerow(header)
writer4.writerow(header)
writer5.writerow(header)

for i in range(len(strain_symptoms) - 1): 
    strain_row = [strain_symptoms[i][0]] #placing the strain number as the first item in the row
    strain_row.append(strain_symptoms[i][1]) #placing the strain symptom associated with that strain from the data
    strain_pres_abs = pres_abs[i] #pulling out the presence absence dictionary for this particular strain
    for gene in homologous_genes: #for each gene in the homologous genes list
        strain_row.append(strain_pres_abs[gene]) #placing the value for each homologous gene in the presence absence dictionary for this particular strain
    writer5.writerow(strain_row) #writing the presence absence information from this strain to the csv file as a row
    #Here I can now separate out the strains based on their symptom status and write out to separate matricies
    if strain_symptoms[i][1] == "0":
            writer1.writerow(strain_row)
    if strain_symptoms[i][1] == "1":
            writer2.writerow(strain_row)
    if strain_symptoms[i][1] == "2":
            writer3.writerow(strain_row)
    if strain_symptoms[i][1] == "3":
            writer4.writerow(strain_row)

no_luts.close()
oab.close()
uti.close()
uui.close()
out_file.close()






