# Identifying Genes Associated with Symptom Status

## Introduction

Pathogenic bacteria often rely on virulence factors encoded in their genome to infect a host and cause disease. Our knowledge about virulence factors and how they affect the pathogenicity of bacteria remains limited. For example, we do not yet understand how some bacterial strains with genes encoding for known virulence factors do not cause disease, while others without these genes do. Studies in comparative genomics can address this conundrum by comparing symptomatic and asymptomatic genomes  in order to identify genes associated with the presence or absence of symptoms. Therefore, our goal is to  develop a cohesive pipeline that will annotate the genomes from several strains of a single species of bacteria, cluster the strains’ genes based on sequence similarity, and analyze how the homologous genes shared between them can be indicative to the presence or absence of a certain disease. This tool will give the user the ability to quickly, efficiently and accurately perform analysis on genomes of their choice.

## Project Objectives

![Project 5 Overview](https://github.com/lmaskeri/Identifying-Genes-Associated-with-Symptom-Status/blob/main/wiki-images/project5_broad_overview.jpg)


## Test Dataset

This pipeline was tested using the genome assemblies of 66 E. coli strains retrieved from the urinary tract of 66 female patients presenting with four different symptom statuses: overactive bladder (OAB, n = 5), urinary incontinence (UUI, n = 13), urinary tract infection (UTI, n = 42), and no lower urinary tract symptoms (noLUTS, n = 6). This data was taken from the study, “Genomic Survey of E. coli From the Bladders of Women With and Without Lower Urinary Tract Symptoms” (Garretto et al.).
https://www.frontiersin.org/articles/10.3389/fmicb.2020.02094/full


## Programming Language
#### **Python**

## Operating System
#### **Linux**

## Software Required

### **Prokka** 
#### For Genome Annotation
- https://github.com/tseemann/prokka
### 2. **Usearch** 
#### For Sequence Clustering Based on Homology
- https://www.drive5.com/usearch/

## Python Packages

#### os
- https://docs.python.org/3/library/os.html
#### Argparse
- https://docs.python.org/3/library/argparse.html
#### Scipy.stats
- **From Scipy.stats**
  - kendalltau
  - https://www.statisticshowto.com/kendalls-tau/
 #### Pandas
* https://pandas.pydata.org/pandas-docs/stable/getting_started/install.html
 #### Biopython
 - https://biopython.org/wiki/Documentation
   - **From Biopython:**
     - SeqIO
     - https://biopython.org/docs/1.75/api/Bio.SeqIO.html


## Scripts Included

### 1. generate_cluster_pres_abs.py
Python script that generates a presence/absence matrix for genes found in each strain for each cluster generated by UClust.
### 2. statistics.py 
Python script that calculates the correlation coefficents between each cluster's presence/absence values and the control group and generates result files.
### 3. project_wrapper.py
Python script to run all os commands and scripts for the pipeline.

## Output Files
### 1. Three _correlation_results.txt files corresponding to each symptom group
A list of significantly correlated clusters for each symptom and the most positively and most negatively correlated cluster names.
### 2. Three _correlation_matrix.csv files corresponing to each symptom group + control
A matrix holding all the cluster names and their associated correlation value for each symptom and control group.
### 3. ResultsTable.csv
Summary information about the most negatively and most positively correlated clusters, as well as the genes associated with those clusters.

## Installation

In order to run this code from your working directory, use this git command to clone this repository to your workspace:
```
git clone https://github.com/lmaskeri/Identifying-Genes-Associated-with-Symptom-Status
```
Then, change working directories in order to access all files from the cloned repo:
```
cd Identifying-Genes-Associated-with-Symptom-Status
```
**1.	Download Conda**

https://docs.conda.io/en/latest/miniconda.html#linux-installers

Follow this link and download the Miniconda installer for Linux --  match the version of python that you have installed. In this case, we used the Miniconda3-latest-Linux-x86_64.sh installer.

To install, download the installer and place into your working directory. Then run this command and proceed through installation of Miniconda, accepting all defaults:

```
bash Miniconda3-latest-Linux-x86_64.sh
```
Make sure to initalize your shell for conda using this command if you accidentally chose "no" during installation:
```
conda init
```
Run the ``` conda --version``` command to ensure that conda installed properly. If there are any issues (i.e. conda: command not found), please close out of your terminal, restart and run this command using your own home directory in place of $HOME:
```
source $HOME/miniconda3/bin/activate
```
Run the ``` conda --version``` command to ensure that conda was activated if you were getting any errors before.

**2.	Set up Prokka Environment**

Once conda is installed, set up your Prokka virtual environment by running this command:

```
conda create -n prokka-env -c conda-forge -c bioconda prokka
```

Then activate the environment by running this command:

```
conda activate prokka-env
```

**3.	Install Appropriate Packages for Python Scripts**

  - Run these series of commands to install the packages used for this pipeline:
  
 ```conda install biopython```
 
 ```conda install pandas```
 
 ```conda install scipy```
 
**4.	Finally, run the project_wrapper.py script to proceed through the entire pipeline.**

There are two different flags you can use for this pipeline. 

- The **test** flag
   - Use this flag in order to run the test dataset provided with this pipeline.
   
```
python3 project_wrapper.py test
```

OR 

- The **user** flag
   - Use this flag in order to use your own dataset. To do this, you will need to provide 2 things:
     - A folder named “Assemblies” added to the cloned repo directory that includes all assembly files. 
     - A CSV file containing 3 columns: RefSeq Assembly Accession Numbers, Strain Identifiers, Symptoms
     
```
python3 project_wrapper.py user
```

**Important Note for the user flag: Please name your assembly files within the “Assemblies” folder the same as the first column in the “Data_Information.csv” file. Please see an example of what these input files should look like in the “User_Data_Example” folder. Make sure to place both the "Assemblies" folder and the "Data_Information.csv" file that you create within the main working directory of the cloned repo!**


