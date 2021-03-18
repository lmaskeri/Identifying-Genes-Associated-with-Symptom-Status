



import os
from Bio import SeqIO
from Bio import Entrez

os.getcwd()




def getList(filename): #input is txt file of SRR numbers
    with open(filename) as f:
        ECA = [line.rstrip() for line in f] #putting SRR numbers into a list
    return ECA
def moveAnnotations(GCF):

    for num in GCF:
        command2 = "mv "+ num + ".zip ./Annotations"
        os.system(command2)

moveAnnotations(getList("EcoliAccessions.txt"))
