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
#   2JTP 15417 1 G N G C s H1 12.0                                     #
#                                                                      #
# Dependencies:                                                        #
#   forgi                                                              #
#                                                                      #
# Install forgi (any version, script initially written for version 1.1)#
#   https://viennarna.github.io/forgi/                                 #
#   https://www.tbi.univie.ac.at/RNA/documentation.html                #
#                                                                      #
# Authors: Hampus May-18, Noah May-21                                  #
# ---------------------------------------------------------------------#
import operator
import os
import decimal
import sqlite3
import forgi
from forgi.graph import bulge_graph as cgb
import base36

# Input file 1: 210419_source_structures_wo_dash.txt contains PDB id,
# Sequence and lowest energy secondary structure from FRABASE retrieved with script:
# 1A_retrieve_dbn_from_fra.py
# http://rnafrabase.cs.put.poznan.pl
# Put in 'lists/'
inpf1 = "structures_wo_dash_210419.txt"

# Input file 2: BMRB_match_PDB_210419.txt contains the official list of BMRB-PDB matched ID:s
# Download the BMRB-PDB matched list identities from:
# https://bmrb.io/ftp/pub/bmrb/nmr_pdb_integrated_data/adit_nmr_matched_pdb_bmrb_entry_ids.csv
# Put in 'lists/'.
inpf2 = "BMRB_match_PDB_210419.txt"

# Input Database (put in 'db/'):
mcsf_db_path = "MCSF_DB_v5.0.db"

# Nuclei of interest
# interesting_nuclei = ["H6","H8","H2","C6","C8","C2"]
interesting_nuclei = ["H1", "H3", "N1", "N3"]

# Output file:
output_path = "nuclei_imino_21.txt"

# Set to True to print outputs to console
verbose = False


# -------- Functions -------- #


def average_cs_table(_cs_table):
    """ Takes in cs_table from "SELECT * FROM cs WHERE BMRB_entry_ID=(?)"
    Outputs cs_table with averaged values for identical res IDs
    """

    _cs_table = sorted(_cs_table, key=operator.itemgetter(1, 2, 3, 4, 5))

    last_nuc = 0
    cs_table_averaged = []
    cs_res_data = []
    for cs_row in _cs_table:
        bmrb_id = cs_row[0]
        _assembly_id = cs_row[1]
        _entity_id = cs_row[2]
        _res_id = cs_row[3]
        _nt_type = cs_row[4]
        _nuc_type = cs_row[5]

        current_nuc = [bmrb_id, _assembly_id, _entity_id, _res_id, _nt_type, _nuc_type]

        if current_nuc == last_nuc:
            cs_res_data.append(cs_row)
        else:
            if cs_res_data:
                cs_table_averaged = avg_cs(cs_table_averaged, cs_res_data, cs_row)
                cs_res_data = []
            cs_res_data.append(cs_row)

        last_nuc = current_nuc
    cs_table_averaged = avg_cs(cs_table_averaged, cs_res_data, cs_row)
    return cs_table_averaged


def avg_cs(_cs_table_averaged, _cs_res_data, _cs_row):
    shifts = []
    decimals = []
    for cs_res_row in _cs_res_data:
        cs_val = cs_res_row[6]
        cs_dec = abs(decimal.Decimal(str(cs_val)).as_tuple().exponent)
        shifts.append(cs_val)
        decimals.append(cs_dec)
    dec = min(decimals)
    avg_cs = round(sum(shifts)/len(shifts), dec)
    averaged_res = list(_cs_res_data[0][:6]) + [avg_cs] + list(_cs_res_data[0][7:])
    _cs_table_averaged.append(averaged_res)

    _cs_res_data.append(_cs_row)
    return _cs_table_averaged


def max_pdb_id(pdb_ids):
    """Choose maximum PDB ID from iterable collection of PDB ID strings.

    In recent years, all PDB codes are assigned by the PDB from the pool
    of available codes, in sequential ascending order (base-36), without
    reference to the name of the molecule. It may therefore be of
    interest to be able to choose the newest PDB ID out of a list of
    related PDB IDs (newest entry is generally the one with the highest
    PDB ID (in base-36).

    WARNING: This is not guaranteed to work on older PDB entries since
    authors were allowed to pick PDB ID themselves in the past. In many
    cases, however, authors still picked a higher (in base-36) BDB ID
    for their updated depositions.

    Keep this in mind. Use with care.

    Parameters
    ----------
    pdb_ids : Collection[str]
        An iterable collection of strings that follow the standard PDB
        code format: The first character is a numeral in the range 1-9,
        while the last three characters can be either numerals (in the
        range 0-9) or letters (in the range A-Z in the Latin alphabet).

    Returns
    -------
    max_pdb_id_str : str
        The maximum PDB ID of the input collection. In most cases, this
        ID corresponds to the newest deposition.
    """
    pdb_ids_base36 = [base36.loads(_pdb_id) for _pdb_id in pdb_ids]
    max_pdb = max(pdb_ids_base36)
    max_pdb_id_str = base36.dumps(max_pdb).upper()
    return max_pdb_id_str


def reduce_match_list(bmrb_pdb_match_list):
    """Reduce match list into dict with only latest PDB and BMRB entries

    If more than one PDB ID per BMRB ID, choose most recent PDB ID.
    If more than one BMRB ID per PDB ID, choose most recent BMRB ID.

    WARNING: In some edge cases the entry that is picked may not be the
    most recent. In this script, we nevertheless deemed it as the best
    programmatic solution. For more info, read below.

    BMRB IDs
    --------
    Generally, BMRB entries are sequentially chosen. The entry that has
    the highest number as its BMRB ID is therefore deemed the "most
    recent".

    There might be some minor inconsistencies, however (e.g. if entry
    is modified or if a deposition is slightly delayed etc).

    PDB IDs
    -------
    In recent years, all PDB codes are assigned by the PDB from the pool
    of available codes, in sequential ascending order (base-36).

    This is not guaranteed to work on older PDB entries, however, since
    authors back then were allowed to pick PDB ID themselves. In many
    cases, however, authors still picked a higher (in base-36) BDB ID
    for their updated depositions.

    Parameters
    ----------
    BMRB-PDB match table on the format:
        list[list[bmrb_id, pdb_id]]
    Example:
        [['3322', '1E8P'],
         ['3323', '1N72'],
         ...
         ['6342', '5B82']]

        An updated match table can be found at:
        bmrb.io/ftp/pub/bmrb/nmr_pdb_integrated_data/
            adit_nmr_matched_pdb_bmrb_entry_ids.csv

    Returns
    -------
    reduced_match_dict : dict[str: str]
        A dictionary with BMRB IDs as keys and PDB IDs as values. All
        keys and values are unique.
    """
    bmrb_pdb_match_list = sorted(bmrb_pdb_match_list)  # sort by BMRB ID

    pdb_reduced_bmrb_pdb_match_dict = reduce_pdb_ids(bmrb_pdb_match_list)
    bmrb_reduced_pdb_bmrb_match_dict = reduce_bmrb_ids(bmrb_pdb_match_list)

    reduced_match_dict = {}
    for pdb, bmrb in bmrb_reduced_pdb_bmrb_match_dict.items():
        if pdb in pdb_reduced_bmrb_pdb_match_dict.values():
            reduced_match_dict[bmrb] = pdb

    return reduced_match_dict


def reduce_pdb_ids(match_list):
    """If more than one PDB ID per BMRB ID, choose most recent PDB ID.
    """
    bmrb_pdb_match_dict = dict()
    for _bmrb_id, _pdb_id in match_list:
        if _bmrb_id in bmrb_pdb_match_dict:
            bmrb_pdb_match_dict[_bmrb_id].append(_pdb_id)
        else:
            bmrb_pdb_match_dict[_bmrb_id] = [_pdb_id]

    bmrb_pdb_match_dict = {
        bmrb_id: (max_pdb_id(pdb_ids)
                  if (len(pdb_ids) > 1)
                  else pdb_ids[0])
        for bmrb_id, pdb_ids in bmrb_pdb_match_dict.items()
    }

    return bmrb_pdb_match_dict


def reduce_bmrb_ids(match_list):
    """If more than one BMRB ID per PDB ID, choose most recent BMRB ID.
    """
    pdb_bmrb_match_dict = dict()
    for _bmrb_id, _pdb_id in match_list:
        if _pdb_id in pdb_bmrb_match_dict:
            pdb_bmrb_match_dict[_pdb_id].append(_bmrb_id)
        else:
            pdb_bmrb_match_dict[_pdb_id] = [_bmrb_id]

    pdb_bmrb_match_dict = {
        bmrb_id: (max(pdb_ids)
                  if (len(pdb_ids) > 1)
                  else pdb_ids[0])
        for bmrb_id, pdb_ids in pdb_bmrb_match_dict.items()
    }

    return pdb_bmrb_match_dict


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

# dictionary to hold bmrb ids and corresponding pdb ids, generated with
# reduce_match_list function
bmrb_pdb = list(inpf2.readlines())
bmrb_pdb = [line.split() for line in bmrb_pdb]
d_bmrb_pdb = reduce_match_list(bmrb_pdb)

# connect to sql database
conn = sqlite3.connect(src_path + mcsf_db_path)
c = conn.cursor()

# select all from database table 'entry'
c.execute("SELECT * FROM entry")
entry_table = c.fetchall()

# put bmrb id of wanted pdb structure in bmrb_id list,
# check pdb from FRA for dictionary with bmrb-pdb from inpf2 with database
bmrb_ids = []
for row in entry_table:
    if str(row[0]) in d_bmrb_pdb.keys():
        if d_bmrb_pdb[str(row[0])] in pdb_list:
            bmrb_ids.append(row[0])
if verbose:
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
    c.execute("SELECT * FROM rna_seq WHERE BMRB_entry_ID=(?)",
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

if verbose:

    # check bmrb vs. fra base
    print(37*"-" + "\nCheck match between BMRB and FRAbase:\n" + 37*"-" +
          "\n")

# Control that extracted seqs from bmrb match with seqs from FRAbase
filtered_bmrb_id = []  # bmrb ids of structures where sequences match
for i in range(len(bmrb_ids)):

    ind1 = pdb_list.index(d_bmrb_pdb[str(bmrb_ids[i])])
    if verbose:
        print(str(bmrb_ids[i]) + " " + d_bmrb_pdb[str(bmrb_ids[i])])
        print("From FRA:   " + seq_list[ind1])

    if seq_list[ind1].upper() == extr_seqs[i].upper():
        if verbose:
            print("From BMRB:  " + extr_seqs[i] + "         OK!!!")
        filtered_bmrb_id.append(bmrb_ids[i])

    else:
        if verbose:
            print("From BMRB:  " + extr_seqs[i] + "       OH NO!!!")
    if verbose:
        print("            " + dot_brack_list[ind1] + "\n")
if verbose:
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
if verbose:
    print(27*"-" + "\nChosen nuclei of interest:\n" + "-"*27 + col_names)

# save output to file
path = src_path + output_path
text_file = open(path, "w")

for i in range(len(filtered_bmrb_id)):
    
    # for each individual structure
    c.execute("SELECT * FROM cs WHERE BMRB_entry_ID=(?)",
              (filtered_bmrb_id[i],))
    _cs_table = c.fetchall()

    cs_table = average_cs_table(_cs_table)

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
            if verbose:

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
            if verbose:

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
            if verbose:

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
