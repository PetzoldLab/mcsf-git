import os
import sqlite3

import forgi
from forgi.graph import bulge_graph as cgb

# ---------------------------------------------------------------------#
# Script for extracting chemical shift data from the RNA CS database   #
# CSF_DB_vX.X.db. The script requires input lists described below. It  #
# will generate an output file "nuclei.txt" with the following struct: #
#                                                                      #
# Columns:                                                             #
#   PDB ID, BMRB ID, residue ID, base, base before, base after,        #
#   partner base, feature, nuclei type, chemical shift                 #
#                                                                      #
# Appearance:                                                          #
#   2FDT 10018 1 G N G C s H8 8.145                                    #
#                                                                      #
# Dependencies:                                                        #
#   forgi                                                              #
#                                                                      #
# Install forgi (any version, script initially written for version 1.1)#
#   https://viennarna.github.io/forgi/                                 #
#   https://www.tbi.univie.ac.at/RNA/documentation.html                #
#                                                                      #
# Authors: Hampus May-18, Noah Nov-20                                  #
# ---------------------------------------------------------------------#

# Input file 1: structures_imino_20.txt contains PDB id, Put in 'lists/'
# Sequence and lowest energy secondary structure from FRABASE:
# http://rnafrabase.cs.put.poznan.pl
inpf1 = "structures_imino_20.txt"

# Input file 2: bmrb_pdb_15N_20.txt contains a list of BMRB ID:s and the
# corresponding PDB ID:s, extracted from NMR star files (with
# bmrb2pdb.py in the "main" directory). Put in 'lists/'.
# Download BMRB identities from:
# https://bmrb.io/search/query_grid/query_1_20.html
inpf2 = "bmrb_pdb_15N_20.txt"

# Input Database (put in 'db/'):
mcsf_db_path = "MCSF_DB_v3.0.db"

# Nuclei of interest
# interesting_nuclei = ["H6","H8","H2","C6","C8","C2"]
interesting_nuclei = ["H1", "H3", "N1", "N3"]

# Output file:
output_path = "nuclei_imino_20.txt"

# -------- RUN -------- #

# Import data
split_path = os.path.abspath(os.getcwd()).split(os.path.sep)
src_path = os.path.sep
for a in range(len(split_path[1:-1])):
    src_path = src_path + split_path[a + 1] + os.path.sep
    
inpf1 = open(os.path.normpath(os.path.join(src_path, "lists", inpf1)))
inpf2 = open(os.path.normpath(os.path.join(src_path, "lists", inpf2)))
mcsf_db_path = os.path.normpath(os.path.join("db", mcsf_db_path))
output_path = os.path.normpath(os.path.join("lists", output_path))

# create lists to hold pdb id, sequence and dot_bracket structures
pdb_list, seq_list, dot_brack_list = [], [], []
for i in inpf1.readlines():
    line = i.split()
    pdb_list.append(line[0])
    seq_list.append(line[1])
    dot_brack_list.append(line[2])

# dictionary to hold bmrb ids and corresponding pdb ids
bmrb_pdb = list(inpf2.readlines())
d_bmrb_pdb = {}
for j in bmrb_pdb:
    pairs = j.split()
    d_bmrb_pdb[pairs[0]] = pairs[1]

# capitalize all the nucleotides MR
for number in range(len(seq_list)):
    seq_list[number] = seq_list[number].upper()

# connect to sql database
conn = sqlite3.connect(src_path + mcsf_db_path)
c = conn.cursor()

# select all from database table 'entry'
c.execute("SELECT * FROM entry")
entry_table = c.fetchall()

# put bmrb id of wanted pdb structure in bmrb_id list
bmrb_ids = []
for row in entry_table:
    if str(row[0]) in d_bmrb_pdb.keys():
        if d_bmrb_pdb[str(row[0])] in pdb_list:
            bmrb_ids.append(row[0])

# print out bmrb and pdb ids
print(27*"-" + "\nFound BMRB IDs in DB for:\n" + "-"*27 +
      "\nBMRB    PDB")

for bmrb_id in bmrb_ids:
    print(str(bmrb_id) + " = " + d_bmrb_pdb[str(bmrb_id)])

print("\nSummary:\nFound BMRB ID for " + str(len(bmrb_ids)) +
      " out of " + str(len(pdb_list)) + " PDB structures\n")

# collect all nuclei associated with these bmrb/pdb structures
extr_seqs = []
for i in range(len(bmrb_ids)):
    
    # to hold sequence and base letters
    sequence = []
    resid = 0
    dummy_str = ""
    
    # for each individual structure
    c.execute("SELECT * FROM cs WHERE BMRB_entry_ID=(?)",
              (bmrb_ids[i],))
    cs_table = c.fetchall()
    
    # extract sequence from returned bmrb entities
    for row in cs_table:
        if row[3] != resid:
            sequence.append(row[4])
            resid = int(row[3])
    
    # put all extracted sequences in "extr_seqs"
    for letter in sequence:
        dummy_str = dummy_str + letter
    extr_seqs.append(dummy_str)

conn.commit()
conn.close()

# check bmrb vs. fra base
print(37*"-" + "\nCheck match between BMRB and FRAbase:\n" + 37*"-" +
      "\n")

# Control that extracted seqs from bmrb match with seqs from FRAbase
filtered_bmrb_id = []  # bmrb ids of structures where sequences match
for i in range(len(bmrb_ids)):
    
    ind1 = pdb_list.index(d_bmrb_pdb[str(bmrb_ids[i])])
    print(str(bmrb_ids[i]) + " " + d_bmrb_pdb[str(bmrb_ids[i])])
    print("From FRA:   " + seq_list[ind1])
    
    if seq_list[ind1] == extr_seqs[i]:
        print("From BMRB:  " + extr_seqs[i] + "         OK!!!")
        filtered_bmrb_id.append(bmrb_ids[i])
    
    else:
        print("From BMRB:  " + extr_seqs[i] + "       OH NO!!!")
    
    print("            " + dot_brack_list[ind1] + "\n")

print("Summary:\n" + str(len(filtered_bmrb_id)) +
      " sequences OK out of " + str(len(bmrb_ids)) + "\n")

# connect to database again
conn = sqlite3.connect(src_path + mcsf_db_path)
c = conn.cursor()

# Collect all nuclei associated with the filtered bmrb/pdb structures:

col_names = ("\nPDB\t\t\t"
             "BMRB\t\t"
             "Res\t\t"
             "Base\t"
             "BF\t\t"
             "AF\t\t"
             "Pair\t"
             "feat\t"
             "nuc\t\t"
             "cs")

print(27*"-" + "\nChosen nuclei of interest:\n" + "-"*27 + col_names)

# save output to file
text_file = open(src_path + output_path, "w")

for i in range(len(filtered_bmrb_id)):
    
    # for each individual structure
    c.execute("SELECT * FROM cs WHERE BMRB_entry_ID=(?)",
              (filtered_bmrb_id[i],))
    cs_table = c.fetchall()
    
    for nuclei_data in cs_table:
    
        residue_id  = nuclei_data[3]
        nt_type     = nuclei_data[4]
        nuc_type    = nuclei_data[5]
        cs_val      = nuclei_data[6]
        
        fra_seq = seq_list[
            pdb_list.index(d_bmrb_pdb[str(filtered_bmrb_id[i])])
        ]
        dot_brack = dot_brack_list[
            pdb_list.index(d_bmrb_pdb[str(filtered_bmrb_id[i])])
        ]
        
        # if forgi 1.x
        if forgi.__version__.startswith("1"):
            bg = cgb.BulgeGraph()
            bg.from_dotbracket(dot_brack)
            es = bg.to_element_string()
        # if forgi 2.x
        else:
            bg = forgi.graph.bulge_graph.BulgeGraph
            bg = bg.from_dotbracket(dot_brack)
            es = bg.to_element_string()
        
        # define base pairing partner
        partner = ""
        if bg.pairing_partner(residue_id) is None:
            partner = "N"
        else:
            partner = fra_seq[bg.pairing_partner(residue_id) - 1]

        # If first nucleotide in a sequence
        if nuc_type in interesting_nuclei and residue_id == 1:
            
            # PDB ID, BMRB ID, Residue ID, Base type, Base before
            # Base after, partner base, feature, Nuc type, CS value
            print(
                d_bmrb_pdb[str(filtered_bmrb_id[i])] + "\t\t" +
                str(filtered_bmrb_id[i]) + "\t\t" +
                str(residue_id) + "\t\t" +
                nt_type + "\t\t" + "N" + "\t\t" +
                fra_seq[1] + "\t\t" +
                partner + "\t\t" +
                es[0] + "\t\t" +
                nuc_type + "\t\t" +
                str(cs_val)
            )
            
            text_file.write(
                d_bmrb_pdb[str(filtered_bmrb_id[i])] + " " +
                str(filtered_bmrb_id[i]) + " " +
                str(residue_id) + " " +
                nt_type + " " +
                "N" + " " +
                fra_seq[1] + " " +
                partner + " " +
                es[0] + " " +
                nuc_type + " " +
                str(cs_val) + "\n"
            )

        # if not the first or the last nucleotide in sequence
        elif (nuc_type in interesting_nuclei
              and 1 < residue_id < len(fra_seq)):
    
            # PDB ID, BMRB ID, Residue ID, Base type, Base before
            # Base after, partner base, feature, Nuc type, CS value
            print(
                d_bmrb_pdb[str(filtered_bmrb_id[i])] + "\t\t" +
                str(filtered_bmrb_id[i]) + "\t\t" +
                str(residue_id) + "\t\t" +
                nt_type + "\t\t" +
                fra_seq[residue_id - 2] + "\t\t" +
                fra_seq[residue_id] + "\t\t" +
                partner + "\t\t" +
                es[residue_id - 1] + "\t\t" +
                nuc_type + "\t\t" +
                str(cs_val)
            )
            
            text_file.write(
                d_bmrb_pdb[str(filtered_bmrb_id[i])] + " " +
                str(filtered_bmrb_id[i]) + " " +
                str(residue_id) + " " +
                nt_type + " " +
                fra_seq[residue_id - 2] + " " +
                fra_seq[residue_id] + " " +
                partner + " " +
                es[residue_id - 1] + " " +
                nuc_type + " " +
                str(cs_val) + "\n"
            )

        # if last nucleotide in sequence
        elif (nuc_type in interesting_nuclei
              and residue_id == len(fra_seq)):
    
            # PDB ID, BMRB ID, Residue ID, Base type, Base before
            # Base after, partner base, feature, Nuc type, CS value
            print(
                d_bmrb_pdb[str(filtered_bmrb_id[i])] + "\t\t" +
                str(filtered_bmrb_id[i]) + "\t\t" +
                str(residue_id) + "\t\t" +
                nt_type + "\t\t" +
                fra_seq[residue_id - 2] + "\t\t" +
                "N" + "\t\t" +
                partner + "\t\t" +
                es[residue_id - 1] + "\t\t" +
                nuc_type + "\t\t" +
                str(cs_val)
            )
            
            text_file.write(
                d_bmrb_pdb[str(filtered_bmrb_id[i])] + " " +
                str(filtered_bmrb_id[i]) + " " +
                str(residue_id) + " " +
                nt_type + " " +
                fra_seq[residue_id - 2] + " " +
                "N" + " " +
                partner + " " +
                es[residue_id - 1] + " " +
                nuc_type + " " +
                str(cs_val) + "\n"
            )

text_file.close()
conn.commit()
conn.close()
