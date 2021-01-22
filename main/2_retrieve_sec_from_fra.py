# This script retrieves the PDB ID, sequence and dot-bracket structures from FRABASE
# for the PDB and BMRB matched IDs
# Authors: Hampus May-18

import urllib2
import time

inpf1 = open('./bmrb_pdb_15N_20.txt')

pdb_list = []

for i in inpf1.readlines():
    line = i.split()
    pdb_list.append(line[1])


for j in range(len(pdb_list)):

    url = 'http://rnafrabase.cs.put.poznan.pl/?act=pdbdetails&id='+pdb_list[j]

    response = urllib2.urlopen(url)
    source = response.read()
    source_lines = source.splitlines()
    statements = source_lines[143].split()

    hold_list = []

    for k in range(len(statements)):

        statement = statements[k]

        if statement.find('((') != -1 and statement.find('))') != -1:
            hold_list.append(statement)

    if len(hold_list) > 0:
        low_e_sec_struc_string = hold_list[0]

        br = [i for i in range(len(low_e_sec_struc_string)) if low_e_sec_struc_string.startswith('<br>', i)]

        sequence = low_e_sec_struc_string[(br[0]+4):(br[1])]
        sec_struc = low_e_sec_struc_string[(br[1]+4):low_e_sec_struc_string.find('<', (br[1]+4))]

        response.close()

        print(pdb_list[j]+" "+sequence+" "+sec_struc)

    time.sleep(1)
