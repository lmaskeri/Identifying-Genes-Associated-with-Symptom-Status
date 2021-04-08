# Identifying Genes Associated with Symptom Status

## Installation

In order to run this code from your working directory, use this git command to clone this repository to your workspace:
```
git clone https://github.com/lmaskeri/Identifying-Genes-Associated-with-Symptom-Status
```
Then, change working directories in order to access all files from the cloned repo:
```
cd Identifying-Genes-Associated-with-Symptom-Status
```

## Directions

### 1. Download the Dataset Using NCBI's "datasets" Tool

First, run the curl command to obtain "datasets" from NCBI:
```
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets'
```

Next, run the chmod command to get permission and enable execution of "datasets" for Linux:
```
chmod +x datasets
```
To obtain the test dataset, run the getAnnotations python script which will use the ecoliAccessions.txt file to pull all assembly and annotation files for each test sample:
```
python3 getAnnotations.py
```
### Output:
Individual annotations for each strain in Annotations folder
ProteinSeqs folder with a file of protein sequences of every CDS for each strain
ProteinSeqs.faa file that contains all protein sequences for every strain (for use in uClust)
### 2. Obtain Clusters and Centroid Sequences Using Usearch

Run the usearch command to obtain clusters: 
```
./usearch_lin -cluster_fast proteinSeqs.faa -id 0.90 -centroids nr.fasta -clusters cluster_dir/c_
```
Usearch uses the cluster_fast method in order to produce similar gene clustering accross all sample strains. The centroids used to define these clusters are the homologous genes that we will be looking at.

#### Output:
* Cluster files containing genes that cluster together using the cluster_fast method
* Centroid sequence file (nr.fasta) which contains the centroid sequences for each cluster

### 3. Generate Presence/Absence Matrix:

To generate a presence/absence matrix for the entire dataset, as well as separate matricies for each symptom sample, run the generate_pres_abs_matrix python script:
```
python3 generate_pres_abs_matrix.py
```
This script creates a list of homologous genes from the centroid fasta file, then parses through each protein file per sample to identify which homologous gene names are present or absent for each sample.

#### Output:
* Samples with no_luts presence absence matrix
* Samples with OAB presence absence matrix
* Samples with UTI presence absence matrix
* Samples with UUI presence absence matrix
* All Samples presence absence matrix
 
