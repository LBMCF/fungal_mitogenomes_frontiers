#This script is part of supplementary documents of "XXX"
#submitted to Frontiers in Microbiology, section XXX
#Manuscript ID: XXX
#Authors:  Paula L. C. Fonseca, Ruth B. De-Paula, Daniel S. Araújo, Luiz Marcelo Ribeiro Tomé, Thairine Mendes-Pereira, Wenderson Felipe Costa Rodrigues, 
#Luiz-Eduardo Del-Bem, Eric R. G. R. Aguiar, and Aristóteles Góes-Neto

# This script correctly computes GC content and length of a sequence in FASTA format, or multiple sequences in a multi-FASTA file.  

#******************************************************************************#
#                            Run the code in Python3                           #
#******************************************************************************#

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC

# Defines a function that takes an input file containing sequences in the .fasta format,
# calculates the length and GC content of every sequence and saves the result in an output file
def lenGC(input_f,output_f):
    sequences = list(SeqIO.parse(input_f,"fasta")) #Creates a list containing sequences in .fasta format
    for item in sequences: #For every sequence in the list
        print(f"{item.description}\t{len(item.seq)}\t{GC(item.seq)}") #Prints its ID, length, and GC content
        output_f.write(f"{item.description}\t{len(item.seq)}\t{GC(item.seq)}\n") #Writes its ID, length, and GC content in the output file

input_file2= #Takes the name of the input file
input_file=open(f"{input_file2}","r")
output_file=open(f"{input_file2}"+"_lenGC.txt","a") #Opens output file in append mode
lenGC(input_file,output_file) #Runs lenGC function
output_file.close() #Closes output file