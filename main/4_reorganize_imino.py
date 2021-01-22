
# Script for reorganizing extracted chemical shift data from "extract.py"
# The script double checks that both proton and nitrogen data are there and
# converts the two chemical shifts into one numerical value
# It reads "nuclei.txt" from extract.py and outputs "nuc_sparse.txt"
# Outliers have to be handled manually
# Author: Hampus May-18

import os
import numpy as np

# open file part
split_path = os.path.abspath(os.getcwd()).split('/')
path_above = '/'

for a in range(len(split_path[1:-1])):
    path_above = path_above + split_path[a+1]+'/'

# open nuclei file
inpf1 = open(path_above+"lists/nuclei_imino_20.txt")

# each line in a list
lines_inpf1 = inpf1.readlines()

# the nuclei we are interested in both H and C
protons, carbons = ["H1", "H3"], ["N1", "N3"]

cv1 = 0

print(45*"-"+"\nThe following residues has both H and N data\n"+45*"-")

H_CS = []
N_CS = []
re_organized = []

# 2O33 15080 4 U U U N i C2 149.808

# loop through file and see if both H and C are there
for line in range(len(lines_inpf1)):

    if line < (len(lines_inpf1)-1):

        columns = lines_inpf1[line].split()
        line_re_org = ""

        if columns[8] == "H1" and line < len(lines_inpf1)-3:
            for j in range(3):
                if lines_inpf1[line+(j+1)].split()[8] == "N1" and columns[0:7] == lines_inpf1[line+(j+1)].split()[0:7]:
                    line_re_org = str(cv1+1)+" "+columns[0]+" "+columns[1]+" "+columns[2] + \
                                  " "+columns[3]+" "+columns[4]+" "+columns[5]+" "+columns[6] + \
                                  " "+columns[7]+" "+columns[8]+" "+lines_inpf1[line+(j+1)].split()[8] + \
                                  " "+columns[9]+" "+lines_inpf1[line+(j+1)].split()[9]

                    H_CS.append(float(columns[9]))
                    N_CS.append(float(lines_inpf1[line+(j+1)].split()[9]))
                    re_organized.append(line_re_org)
                    cv1 = cv1+1

        elif columns[8] == "N1" and line < len(lines_inpf1)-3:
            for j in range(3):
                if lines_inpf1[line+(j+1)].split()[8] == "H1" and columns[0:7] == lines_inpf1[line+(j+1)].split()[0:7]:
                    line_re_org = str(cv1+1)+" "+columns[0]+" "+columns[1]+" "+columns[2] + \
                                  " "+columns[3]+" "+columns[4]+" "+columns[5]+" "+columns[6] + \
                                  " "+columns[7]+" "+columns[8]+" "+lines_inpf1[line+(j+1)].split()[8] + \
                                  " "+columns[9]+" "+lines_inpf1[line+(j+1)].split()[9]

                    N_CS.append(float(columns[9]))
                    H_CS.append(float(lines_inpf1[line+(j+1)].split()[9]))
                    re_organized.append(line_re_org)
                    cv1 = cv1+1

        elif columns[8] == "N3" and line < len(lines_inpf1)-3:
            for j in range(3):
                if lines_inpf1[line+(j+1)].split()[8] == "H3" and columns[0:7] == lines_inpf1[line+(j+1)].split()[0:7]:
                    line_re_org = str(cv1+1)+" "+columns[0]+" "+columns[1]+" "+columns[2] + \
                                  " "+columns[3]+" "+columns[4]+" "+columns[5]+" "+columns[6] + \
                                  " "+columns[7]+" "+columns[8]+" "+lines_inpf1[line+(j+1)].split()[8] + \
                                  " "+columns[9]+" "+lines_inpf1[line+(j+1)].split()[9]

                    N_CS.append(float(columns[9]))
                    H_CS.append(float(lines_inpf1[line+(j+1)].split()[9]))
                    re_organized.append(line_re_org)
                    cv1 = cv1+1

        elif columns[8] == "H3" and line < len(lines_inpf1)-3:
            for j in range(3):
                if lines_inpf1[line+(j+1)].split()[8] == "N3" and columns[0:7] == lines_inpf1[line+(j+1)].split()[0:7]:
                    line_re_org = str(cv1+1)+" "+columns[0]+" "+columns[1]+" "+columns[2] + \
                                  " "+columns[3]+" "+columns[4]+" "+columns[5]+" "+columns[6] + \
                                  " "+columns[7]+" "+columns[8]+" "+lines_inpf1[line+(j+1)].split()[8] + \
                                  " "+columns[9]+" "+lines_inpf1[line+(j+1)].split()[9]

                    H_CS.append(float(columns[9]))
                    N_CS.append(float(lines_inpf1[line+(j+1)].split()[9]))
                    re_organized.append(line_re_org)
                    cv1 = cv1+1

        if line_re_org != "":
            print(line_re_org)

h_cs = np.asarray(H_CS)
n_cs = np.asarray(N_CS)


# 3. filter off outliers

final_nuc_before_filter_outliers = []

for k in range(len(re_organized)):

    final_nuc_before_filter_outliers.append(re_organized[k])

bmrb_id_of_sets_with_outliers = []

for k in final_nuc_before_filter_outliers:

    s_line = k.split()

    H_chem_shift = float(s_line[11])
    N_chem_shift = float(s_line[12])

    if (20.0 < H_chem_shift or H_chem_shift < 6.0) and (s_line[2] not in bmrb_id_of_sets_with_outliers):
        bmrb_id_of_sets_with_outliers.append(s_line[2])

    if (200.0 < N_chem_shift or N_chem_shift < 80.0) and (s_line[2] not in bmrb_id_of_sets_with_outliers):
        bmrb_id_of_sets_with_outliers.append(s_line[2])

absoulutely_final_clean_sensical_list = []

for i in final_nuc_before_filter_outliers:

    split_line = i.split()

    if split_line[2] not in bmrb_id_of_sets_with_outliers:
        absoulutely_final_clean_sensical_list.append(i)

print(bmrb_id_of_sets_with_outliers)

# 4. Now write all (reliable = unchanged + calibrated) structures to output file.

text_file = open(path_above+"lists/imino_nuc_sparse_20.txt", "w")

# first write all reliable
for k in absoulutely_final_clean_sensical_list:

    print(k)
    text_file.write(k+"\n")

text_file.close()
