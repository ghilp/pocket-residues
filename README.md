# pocket-residues
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 13 17:20:44 2023

@author: ghilp
"""

import re
import os
import csv
from pymol import cmd
import numpy as np

# Specify notbook name
notebook_path = os.path.abspath("pocket_residues.py")

# Get the path to the directory containing the script
script_dir = os.path.dirname(notebook_path)

# Look for a .pdb file in the script directory
pdb_file = None
for filename in os.listdir(script_dir):
    if filename.endswith('.pdb'):
        pdb_file = os.path.join(script_dir, filename)
        break

cmd.load(pdb_file, object="protein")

# Display name of loaded pdb file
print(f'{filename} loaded')

# Search for .poc files and get the path of the first one found
pocket_poc = None
for folder in os.listdir(script_dir):
    folder_path = os.path.join(script_dir, folder)
    if os.path.isdir(folder_path):
        for file in os.listdir(folder_path):
            if file.endswith(".poc"):
                pocket_poc = os.path.join(folder_path, file)
                break
        if pocket_poc is not None:
            break

# Converts .poc file to .csv
def convert_to_csv(file_path):
    with open(file_path, 'r') as file:
        content = file.read()
        content = re.sub(' +', ',', content)
        content = content.replace('\n', ',\n')

    csv_file_path = file_path.split('.')[0] + '.csv'
    with open(csv_file_path, 'w') as file:
        file.write(content)

# Create a .csv file
convert_to_csv(pocket_poc)

# Search for .csv file with same name as .poc file in the same directory
pocket_csv = None
if pocket_poc is not None:
    poc_filename = os.path.basename(pocket_poc)
    csv_filename = os.path.splitext(poc_filename)[0] + ".csv"
    csv_path = os.path.join(os.path.dirname(pocket_poc), csv_filename)
    if os.path.isfile(csv_path):
        pocket_csv = csv_path
        
def extract_residue_info(csv_file_path):
    with open(csv_file_path, 'r') as file:
        lines = file.readlines()

    residue_dict = {}
    for line in lines:
        split_line = line.strip().split(',')
        residue_dict[split_line[5]] = split_line[3]

    residue_numbers = list(residue_dict.keys())
    residue_numbers.sort(key=float)

    return residue_dict, residue_numbers

residue_dict, residue_numbers = extract_residue_info(pocket_csv)

residue_codes = []

# Display total pocket residues
print(f"{len(residue_numbers)} total pocket residues identified.")

# Find the center of mass for the protein and call it COM
COM= cmd.centerofmass("protein", state=1)

# Define a threshold
threshold = 0

# Define a list of inward facing residues
residue_list = []

# Iterate over the residue numbers
for residue in residue_numbers:
    try:
        # Calculate the center of mass for the C-alpha carbon of the residue
        cmd.select("C_alpha_residue", "protein and resi " + str(residue) + " and name CA")
        residue_COM= cmd.centerofmass("C_alpha_residue", state=1 )
        # Calculate the centroid of the residue
        #cmd.select('res',"protein and resi " + str(residue))
        #sidechaincenters ("res","all", name="model_centroid",  method= 'centroid')
        cmd.select('sidechainA', 'protein and resi ' + str(residue) + ' and not name N+C+O+H')

        # Calculate the centroid coordinates of the selected side chain atoms
        residue_centroid = cmd.centerofmass('sidechainA')
    
        # Measure the distance between COM and residue_COM
        d1=np.linalg.norm(np.array(COM)-np.array(residue_COM))
   
        # Measure the distance between COM and residue_centroid
        d2= np.linalg.norm(np.array(COM)-np.array(residue_centroid))
        
        # If the difference between d2 and d1 is less than the threshold, print the residue number
        if d2-d1 < threshold:
            residue_list.append(residue)
    except Exception as e:
        print(f"Error processing residue {residue}: {str(e)}")

# Match residue numbers with names
for residue in residue_list:
    residue_codes.append(residue_dict[residue])

# Create CSV file with final residue list
# Define the header row
header = ["Residue", "Number"]

# Combine the header row and data
rows = [header]
for code, number in zip(residue_codes, residue_list):
    rows.append([code, number])

# Write the data to a CSV file
with open("residues.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerows(rows)
    
print(f"Data successfully saved to residues.csv for {len(residue_list)} residues.")

# Define a threshold for minimum number of residues
min_residue_threshold = 10

# Check if the residue list is below the threshold and print an error message if it is
if len(residue_list) < min_residue_threshold:
    print(f"Error: The number of inward facing residues ({len(residue_list)}) is below the threshold of {min_residue_threshold}.")   
