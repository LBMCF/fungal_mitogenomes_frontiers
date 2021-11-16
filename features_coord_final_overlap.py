#This script is part of supplementary documents of "Global Characterization of Fungal Mitogenomes: New Insights on Genomic Diversity and Dynamism of Coding Genes and Accessory Elements"
#published in Frontiers in Microbiology, research topic "Mitochondrial Genomes and Mitochondrion Related Gene Insights to Fungal Evolution"
#DOI: 10.3389/fmicb.2021.787283
#Authors:  Paula L. C. Fonseca, Ruth B. De-Paula, Daniel S. Araújo, Luiz Marcelo R. Tomé, Thairine Mendes-Pereira, Wenderson F. Rodrigues, 
#Luiz-Eduardo Del-Bem, Eric R. Aguiar, and Aristóteles Góes-Neto
#External help: Luciano Branco

#This script correctly calculates the total genic/coding region length of each mitogenome studied, already discounting the lengths of two features that spatially overlap.
#The files used are as follow: 
# Create a new folder named "input_folder" and add to this folder all the files containing genomic coordinates of features (formatted like .bed file).

#******************************************************************************#
#                            Run the code in Python3                           #
#******************************************************************************#

# coding=utf-8

import glob
import pandas as pd
import numpy as np
import numpy.matlib
import os

def get_filenames(folder):
    fnames_dir = glob.glob(folder + "/*")
    # fnames = [(fname_dir.split(folder+)[1]) for fname_dir in fnames_dir]
    fnames = [(os.path.split(fname_dir))[1] for fname_dir in fnames_dir]
    
    return fnames, fnames_dir
    

input_folder = 'input'
output_folder = 'output'

fnames, fnames_dir = get_filenames(input_folder)

# fname = fnames[0]
for fname, fname_dir in zip(fnames, fnames_dir):

    input_file_coordinates = pd.read_csv(fname_dir, delimiter = '\t', header=None) 
    input_file_coordinates = input_file_coordinates.to_numpy() # Vectorize for performance
    
    all_coords = []
    
    for coord in input_file_coordinates:
        coord_id = coord[0]
        coord_type = coord[2]
        
        # Ignore anything before ":" in coord_id
        coords = coord_id.split(":")[1].split("_")
        
        try:
            # COnvert elements of list to int
            coords = list(map(int, coords))
        except:
            print("Error processing line: " + str(coord))
            continue            
        
        # Sort coords
        coords.sort()
        
        all_coords.append(coords)
        
    max_coord_length = max(max(all_coords)) +1
    # We need to ignore coord 0 and have the vector
    # from 1 to N. Plus, coords 4:5 have size 5-4+1=2
    
    coords_one_hot = np.zeros(max_coord_length , dtype=np.uint8)
    
    for coord in all_coords:
        coords_one_hot[coord[0] : coord[1]+1] = 1 # Matches = 1
    
    total_sum_coords = sum(coords_one_hot)
    
    # Save output
    output_fname = os.path.join(output_folder, fname)
    f = open(output_fname, "w")
    f.write("%s\t%d" % (fname, total_sum_coords))
    f.close()

