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
    transposed = open("transposed.csv","w")

    writer = csv.writer(transposed)
    rez.to_csv(transposed,index=False,header=True)



def kendall(names):
    matrix = pd.read_csv("transposed.csv") #open transposed pres_abs matrix
    correlations = {} #initialize dictionary for keys as cluster names and values as correlation coefficients 
    results1 = open("correlation_results.txt","w") #results text file, will state the highest positive corr. cluster and highest neg. corr. cluster
    results2 = open("correlation_matrix.csv","w") # csv file that will hold all clusters correlation coefficients
    y = matrix["Strain Symptom"].to_list() #make the symptom status into a list 
    for i in names: #for i in cluster name list
        x = matrix[i].to_list() #make the clusters column into a list
        corr,_= kendalltau(x,y) #scipy function for the kendalls tau test
        if str(corr) == "nan": #replace "nan" correlation coefficients with 0
            correlations[i]= 0
        else:
            correlations[i]= corr #append key and value to correlations dictionary
    v=list(correlations.values()) #make a list of correlations values
    k=list(correlations.keys()) #make a list of the correlations jeys 
    positive = k[v.index(max(v))] #find key with max value
    negative = k[v.index(min(v))] #find key with min value 
    results1.write("The cluster with the highest positive correlation to symptom status is " + positive + "\n") #key with max value is the highest positively correlated cluster
    results1.write("The cluster with the highest negative correlation to symptom status is " + negative) #key with min value is the highest negatively correlated cluster
    writer = csv.writer(results2)
    for key,value in correlations.items():
        writer.writerow([key,value]) #write keys and values to a matrix in csv format
    results2.close()
    return correlations
cluster_values = [] #empty list to append cluster filenames to
for cluster_filename in os.listdir("./cluster_dir"):
    cluster_values.append(cluster_filename) #getting a list of cluster file names
    
transpose("cluster_pres_abs_matrix.csv") #transpose matrix
kendall(cluster_values) #run correlation tests
