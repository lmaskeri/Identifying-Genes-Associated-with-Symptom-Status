
'''
Due 05/06/2021
COMP BIO 383 Final Project
Main Wrapper File
'''
import os
import argparse
    

#Time 21.78 minutes for test file
def prokka(flag):
    main_dir = os.getcwd()
    if flag == "test":
        os.system("mkdir Prokka_Annotations")# Making a folder to store all the annotations from prokka
        os.chdir("./Test_Data/Assemblies")# Changing directory to Assemblies for test data
        os.system('ls *.fna | parallel --verbose "prokka {} --prefix {.}_out"')# Run prokka and start annotate all the .fna input 
        os.chdir(main_dir)# Back to the main directory to organize prokka output files
        command1 = "mv ./Test_Data/Assemblies/*_out Prokka_Annotations"# move prokka annotation files into new folder
        os.system(command1)
        os.system("mkdir ProteinSeqs")
        command2 = "find ./Prokka_Annotations -name \*.faa -exec cp {} ProteinSeqs \;" #find all .faa protein files and store them into ProteinSeqs file
        os.system(command2)
    if flag == "user":
        os.system("mkdir Prokka_Annotations")# Making a folder to store all the annotations from prokka
        os.chdir("./Assemblies")# Changing directory to Assemblies for USER data
        os.system('ls *.fna | parallel --verbose "prokka {} --prefix {.}_out"')# Run prokka and start annotate all the .fna input 
        os.chdir(main_dir)# Back to the main directory to organize prokka output files
        command1 = "mv ./Assemblies/*_out Prokka_Annotations"# move prokka annotation files into new folder
        os.system(command1)
        os.system("mkdir ProteinSeqs")
        command2 = "find ./Prokka_Annotations -name \*.faa -exec cp {} ProteinSeqs \;" #find all .faa protein files and store them into ProteinSeqs file
        os.system(command2)
        

def concat_protein():
    command6 = "cat ./ProteinSeqs/*.faa > ProteinSeqs.faa" #merges protein sequences in ProteinSeqs folder into one multi fasta file
    os.system(command6) 
    
#Time 30.0 seconds for test file
def usearch():
    os.system("gunzip usearch_linux") # Before using Usearch, the binary Usearch downloaded from https://drive5.com/usearch/download.html, must be unzipped using this command
    os.system("chmod +x ./usearch_linux") #Then this command must be run to gain access to the binary file
    os.system("mkdir cluster_dir") #a directory must be made for usearch to place the cluster files into
    #Obtaining clusters using Usearch 
    #Takes all .faa files from the ProteinSeqs folder and generates cluster files
    os.system("./usearch_linux -cluster_fast ProteinSeqs.faa -id 0.90 -centroids nr.fasta -clusters cluster_dir/c_")

#Time for test file ~30 minutes
def cluster_pres_abs(flag):
    if flag == "test":
        os.system("python3 generate_cluster_pres_abs.py test")
    if flag == "user":
        os.system("python3 generate_cluster_pres_abs.py user")

#Time for test file ~10 minutes
def calculate_stats():
        os.system("python3 statistics.py")
    



####TO RUN####

#getting working directory of the repo after cloning and moving into repo folder
main_dir = os.getcwd()


#First, setting up argparser in order to take in the argument from the command line
#to either run the test dataset or to run the user dataset
parser = argparse.ArgumentParser(description = "Determine what data to run through script.")
#adding an argument for the specific flag needed to determine which dataset the user wants to run
parser.add_argument('dataset_to_run', type = str, help = "Specify which dataset you'd like to run: test or user")
#parsing the arguments
argument = parser.parse_args()
#argument.dataset_to_run now holds the string from the commandline that will determine if the whole dataset is run or the test dataset is run

#naming the parsed argument, flag
flag = argument.dataset_to_run



##Running Scripts##

#After starting prokka-env ... 

#1. Run Prokka To obtain Annotation Files
print("Starting Prokka to Obtain Annotation Files...")
print("Note: For 66 assemblies, this should take ~25-30 minutes.")
prokka(flag)
print("Prokka Complete.")

#2.
print("Concatenating Protein Files...")
concat_protein()

#3. Generating Clusters using Usearch
print("Generating Clusters using UClust...")
usearch()
print("Clusters Generated Successfully.")

#4. Constructing Presence Absence Cluster Matrix
print("Constructing Presence Absence Cluster Matrix")
print("Note: For ~14000 clusters, this should take ~30 minutes.")
cluster_pres_abs(flag)
print("Presence Absence Cluster Matrix Generated Successfully.")

#5. Running Kendalls Tau Correlation Tests on Clusters
print("Running Kendalls Tau Correlation Tests on Clusters...")
print("Note: For a pres/abs Matrix of ~14000 clusters, this should take ~10 minutes.")
calculate_stats()
print("Correlation Testing Complete. Please Refer to Output Files for Analysis!")
print("Thank you for using our pipeline :) Please deactivate the prokka environment by using the command 'conda deactivate'")



