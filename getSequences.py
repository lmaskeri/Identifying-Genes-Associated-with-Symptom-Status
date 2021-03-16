import os
from Bio import SeqIO
from Bio import Entrez

os.getcwd()


 

def getList(filename): #input is txt file of SRR numbers
    with open(filename) as f:
        ECA = [line.rstrip() for line in f] #putting SRR numbers into a list
    return ECA

def getSequences(EC):
    for i in EC:
        command = "esearch -db assembly -query "+ i + " | elink -target nucleotide -name \
        assembly_nuccore_refseq | efetch -format fasta > " +i + ".fa"
        os.system(command)

getSequences(getList("EcoliAccessions.txt"))
        
