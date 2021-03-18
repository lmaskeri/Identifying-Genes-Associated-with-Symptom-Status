import os
from Bio import SeqIO
from Bio import Entrez

os.getcwd()


 

def getList(filename): #input is txt file of SRR numbers
    with open(filename) as f:
        ECA = [line.rstrip() for line in f] #putting SRR numbers into a list
    return ECA


def getAnnotations(GCF):
    for i in GCF:
        command ="./datasets download genome accession " + i + ".1 --exclude-rna --exclude-protein --exclude-seq --filename " + i +".zip"
        os.system(command)

getAnnotations(getList("EcoliAccessions.txt"))
