# This script outputs the BMRB IDs matched with PDB IDs
# extracted from NMR star files downloaded from
# https://bmrb.io/search/query_grid/query_1_20.html
__author__ = 'hamkar'

import os

f_list = os.listdir("./bmrb_entries_15N_RNA_2020/")

for a in f_list:

    f = open("./bmrb_entries_15N_RNA_2020/"+a)
    lf=f.readlines()
    pdb = ""

    for b in range(len(lf)):

        line = str(lf[b])

        if line.find("PDB") != -1 and line.find("BMRB Entry Tracking System") != -1:
            print a[3:-4]+" "+line.split()[1]
