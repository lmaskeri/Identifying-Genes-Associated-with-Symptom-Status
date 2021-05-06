from Bio import SeqIO
import csv
import os
import pandas as pd
from scipy.stats import kendalltau

def transpose(file): #method for transposing the pres/abs matrix
    m = pd.read_csv(file)
    rez = m.transpose()
    new_header = rez.iloc[0]
    rez =rez[1:]
    rez.columns =new_header #removing header made during the transposition and replacing with the first line of the matrix
    
    return rez


def kendall(names, matrix,code):
    # matrix = pd.read_csv("transposed.csv") #open transposed pres_abs matrix
    
    correlations = {} #initialize dictionary for keys as cluster names and values as correlation coefficients 
    correlations_p = {} #empty dictionary to hold p-values from correlation significance tests

    results1 = open(code+"_correlation_results.txt","w") #results text file, will state the highest positive corr. cluster and highest neg. corr. cluster
    results2 = open(code+"_correlation_matrix.csv","w") # csv file that will hold all clusters correlation coefficients
    y = matrix["Strain Symptom"].to_list() #make the symptom status into a list 
    for i in names: #for i in cluster name list
        x = matrix[i].to_list() #make the clusters column into a list
        corr,p= kendalltau(x,y) #scipy function for the kendalls tau test
        if str(corr) == "nan": #replace "nan" correlation coefficients with 0
            correlations[i]= 0
            correlations_p[i]=p #append p value to dictionary as value
        else:
            correlations[i]= corr #append key and value to correlations dictionary
            correlations_p[i]=p #append p value to dictionary as value
    v=list(correlations.values()) #make a list of correlations values
    k=list(correlations.keys()) #make a list of the correlations jeys 
    positive = k[v.index(max(v))] #find key with max value
    negative = k[v.index(min(v))] #find key with min value
    pos_corr = correlations[positive] #getting positive correlation value
    neg_corr = correlations[negative] #getting negative correlation value 


    results1.write("The clusters whose correlation coefficients are statistically significant are: \n")
    count = 0
    for key,value in correlations_p.items():
        if value <=0.05: #significance level is 0.05
            results1.write(key + "\n") #if p-value above 0.05, reject null hypothesis, the correlation value is statistically significant and there is likely an association present
            count += 1 #counting how many clusters are significantly correlated
    results1.write("The cluster with the highest positive correlation to symptom status is " + positive + "\n") #key with max value is the highest positively correlated cluster
    results1.write("The cluster with the highest negative correlation to symptom status is " + negative) #key with min value is the highest negatively correlated cluster
    writer = csv.writer(results2)
    for key,value in correlations.items():
        writer.writerow([key,value]) #write keys and values to a matrix in csv format
    results2.close()
    return correlations,[positive,pos_corr,negative,neg_corr],count 

def findGeneName(name):
    genes = []
    take_out = [ "[Proteobacteria]", "[Escherichia]", "[Escherichia coli]", "[Enterobacterales]","[Enterobacteriaceae]", "[Bacteria]", "[Gammaproteobacteria]"] #Names for E.coli classes that we want to take out of the gene name
    
    with open("./cluster_dir/"+name) as cluster: #open centroid file
        for record in SeqIO.parse(cluster, "fasta"):
            gene_name = record.description[15:] #pull out description (not using id because we don't want the accession number)
            if "MULTISPECIES:" in gene_name: #taking out multispecies
                gene_name = gene_name[14:]
            for out in take_out: #taking out the E.coli class
                if out in gene_name:
                    end = len(out) + 1
                    gene_name = gene_name[:-end]
        genes.append(gene_name) #appending the getting a list of all the gene names found in the clutser to a list
    counter = 0
    num = genes[0]

    for i in genes:#going through the list of gene names to find the most common
        curr_freq = genes.count(i)
        if(curr_freq>counter):
            counter=curr_freq
            num = i
    return num #returning most common gene name
def makeResultsList(code,pos_neg,pos_gene_name,neg_gene_name,count): #code is symptom name, count is number of correlated clusters
    pos_name = pos_neg[0] #grabbing necessary items from list returned in kendall() method
    pos_corr = pos_neg[1] 

    neg_name = pos_neg[2]
    neg_corr = pos_neg[3]

    results_list = [code, pos_name, pos_corr,pos_gene_name, neg_name, neg_corr,neg_gene_name,count] #making a list of needed components for each symptom for use in final results table
    return results_list

def makeResultsTable(list1,list2,list3):
    data =[list1,list2,list3] #making a list of lists of all symptom lists made in makeResultsList method
    df = pd.DataFrame(data,columns = ["Symptom","Positive Name","Positive Correlation Value", "Positive Gene Name","Negative Name","Negative Correlation Value", "Negative Gene Name","Number of Significantly Correlated Clusters"])
    df.to_csv("ResultsTable.csv",index=False) #turning into a pandas dataframe and then into a csv file

#pulling out the names of all clusters using file directory
cluster_values = [] #empty list to append cluster filenTeddames to
for cluster_filename in os.listdir("./cluster_dir"):
    cluster_values.append(cluster_filename) #getting a list of cluster file names
    
r = transpose("cluster_pres_abs_matrix.csv") #transpose full matrix

#Splitting the matrix up into each symptom with the added control no LUTS

control = r.loc[r['Strain Symptom'] == "0"] #pulling out all of the control patient sample rows
symp1 = r.loc[r['Strain Symptom'] == "1"] #pulling out all of the symptom 1 patient sample rows
symp2 = r.loc[r['Strain Symptom'] == "2"] #pulling out all of the symptom 2 patient sample rows
symp3 = r.loc[r['Strain Symptom'] == "3"] #pulling out all of the symptom 3 patient sample rows

# #combining a new dataframe with each symptom and the control group for correlation testing

final_symp1 = pd.concat([symp1, control], ignore_index=True, sort=False)
final_symp2 = pd.concat([symp2, control], ignore_index=True, sort=False)
final_symp3 = pd.concat([symp3, control], ignore_index=True, sort=False)

#running kendalls tau for each matrix
symp1_corr,symp1_pos_neg,symp1_count= kendall(cluster_values, final_symp1,"symp1")
symp2_corr,symp2_pos_neg, symp2_count = kendall(cluster_values,final_symp2,"symp2")
symp3_corr,symp3_pos_neg,symp3_count = kendall(cluster_values,final_symp3,"symp3")

#finding gene names for each positive and negative cluster for each symptom
symp1_pos_gene_name = findGeneName(symp1_pos_neg[0])
symp1_neg_gene_name = findGeneName(symp1_pos_neg[2])

symp2_pos_gene_name = findGeneName(symp2_pos_neg[0])
symp2_neg_gene_name = findGeneName(symp2_pos_neg[2])

symp3_pos_gene_name = findGeneName(symp3_pos_neg[0])
symp3_neg_gene_name = findGeneName(symp3_pos_neg[2])

#making ResultsList for each symptom
symp1_list = makeResultsList("Symptom 1",symp1_pos_neg,symp1_pos_gene_name,symp1_neg_gene_name,symp1_count)
symp2_list = makeResultsList("Symptom 2",symp2_pos_neg,symp2_pos_gene_name,symp2_neg_gene_name,symp2_count)
symp3_list = makeResultsList("Symptom 3",symp3_pos_neg,symp3_pos_gene_name,symp3_neg_gene_name,symp3_count)

makeResultsTable(symp1_list,symp2_list,_list)



