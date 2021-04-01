import os
from Bio import SeqIO
from Bio import Entrez

os.getcwd()

def getList(filename): #input is txt file of SRR accession numbers
    with open(filename) as f:
        ECA = [line.rstrip() for line in f] #putting SRR numbers into a list
    return ECA

def getAnnotations(GCF):
    for i in GCF: #for each accession numbers
        command ="./datasets download genome accession " + i + ".1 --exclude-rna  --filename " + i +".zip" #downloading the annotations f>
        os.system(command)
        command2 = "unzip -o "+i+".zip -d Annotations"
        os.system(command2)


def moveAnnotations(GCF):
    for i in GCF:
        command3 = "mv ./Annotations/ncbi_dataset/data/"+i+".1/"+i+"_ASM"+i[6:12]+"v1_genomic.fna Assemblies" #putting all assemblies  fo>
        os.system(command3)
        command4 = "mv ./Annotations/ncbi_dataset/data/"+ i+ ".1/protein.faa " + i + "_protein.faa" #renamed each protein sequence from t>
        os.system(command4)
        command5 = "mv "+ i+ ".1/" + i + "_protein.faa ProteinSeqs" #moving the protein sequence file into a file called ProteinSeqs for >
        os.system(command5)

os.system("mkdir ProteinSeqs") #making a folder for the protein annotations for each accession number
os.system("mkdir Assemblies") #making a folder for Assemblies

getAnnotations(getList("EcoliAccessions.txt"))
moveAnnotations(getList("EcoliAccessions.txt"))
