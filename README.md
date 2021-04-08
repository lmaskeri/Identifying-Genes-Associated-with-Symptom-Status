# Identifying Genes Associated with Symptom Status

## Installation

In order to run this code from your working directory, use this git command to clone this repository to your workspace:
```
git clone https://github.com/lmaskeri/Identifying-Genes-Associated-with-Symptom-Status
.git
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
### 2. Obtain Clusters and Centroid Sequences Using Usearch
Run the usearch command to obtain clusters 
```
./usearch_lin -cluster_fast proteinSeqs.faa -id 0.90 -centroids nr.fasta -clusters cluster_dir/c_
```
### 3. Generate Presence/Absence Matrix:
