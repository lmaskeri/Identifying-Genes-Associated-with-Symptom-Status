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
    results1.write("The clusters whose correlation coefficients are statistically significant are: \n")
    for key,value in correlations_p.items():
        if value <=0.05: #significance level is 0.05
            results1.write(key + "\n") #if l-value above 0.05, reject null hypothesis, the correlation value is statistically significant and there is likely an association present
    results1.write("The cluster with the highest positive correlation to symptom status is " + positive + "\n") #key with max value is the highest positively correlated cluster
    results1.write("The cluster with the highest negative correlation to symptom status is " + negative) #key with min value is the highest negatively correlated cluster
    writer = csv.writer(results2)
    for key,value in correlations.items():
        writer.writerow([key,value]) #write keys and values to a matrix in csv format
    results2.close()
    return correlations

#pulling out the names of all clusters using file directory
cluster_values = [] #empty list to append cluster filenames to
for cluster_filename in os.listdir("./cluster_dir"):
    cluster_values.append(cluster_filename) #getting a list of cluster file names
    
r = transpose("cluster_pres_abs_matrix.csv") #transpose full matrix

#Splitting the matrix up into each symptom with the added control no LUTS

no_luts = r.loc[r['Strain Symptom'] == "NoLUTS"] #pulling out all of the control patient sample rows
oab = r.loc[r['Strain Symptom'] == "OAB"] #pulling out all of the OAB patient sample rows
uti = r.loc[r['Strain Symptom'] == "UTI"] #pulling out all of the UTI patient sample rows
uui = r.loc[r['Strain Symptom'] == "UUI"] #pulling out all of the UUI patient sample rows

# #combining a new dataframe with each symptom and the control group for correlation testing

final_oab = pd.concat([oab, no_luts], ignore_index=True, sort=False)
final_uti = pd.concat([uti, no_luts], ignore_index=True, sort=False)
final_uui = pd.concat([uui, no_luts], ignore_index=True, sort=False)

#running kendalls tau for each matrix
kendall(cluster_values, final_oab,"oab")
kendall(cluster_values,final_uti,"uti")
kendall(cluster_values,final_uui,"uui")


