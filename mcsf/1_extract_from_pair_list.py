#!/usr/local/bin/python
# coding: utf-8

# **************************************************************************** #
# ||                  ~~~  1_extract_from_pair_list.py  ~~~                 || #
# ||                                                                        || #
# ||          Extract CS data from local SQLite3 BMRB RNA DB using          || #
# ||               FRABASE pair list and BMRB-PDB match list.               || #
# ||                                                                        || #
# ||               Written by Noah Hopkins and Magdalena Riad               || #
# ||                         for to the publication                         || #
# ||      "Mutate-and-Chemical-Shift-Fingerprint (MCSF) to characterize     || #
# ||              excited states in RNA using NMR spectroscopy"             || #
# ||                                                                        || #
# ||      Based on earlier work from Lorenzo Baronti and Hampus Karlsson    || #
# ||                                                                        || #
# ||                           Katja Petzold Group                          || #
# ||                                May 2021                                || #
# ||                                                                        || #
# || ---------------------------------------------------------------------- || #
# ||                                                                        || #
# || GNU GPLv3 License                                                      || #
# ||                                                                        || #
# || Mutate-and-Chemical-Shift-Fingerprint (MCSF) Pipeline                  || #
# || Query local RNA BMRB database for chemical shift analysis.             || #
# || Copyright (c) 2021 Noah Hopkins & Magdalena Riad,                      || #
# || Katja Petzold Group, Karolinska Institutet                             || #
# ||                                                                        || #
# || This program is free software: you can redistribute it and/or modify   || #
# || it under the terms of the GNU General Public License as published by   || #
# || the Free Software Foundation, either version 3 of the License, or      || #
# || (at your option) any later version.                                    || #
# ||                                                                        || #
# || This program is distributed in the hope that it will be useful,        || #
# || but WITHOUT ANY WARRANTY; without even the implied warranty of         || #
# || MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          || #
# || GNU General Public License for more details.                           || #
# ||                                                                        || #
# || You should have received a copy of the GNU General Public License      || #
# || along with this program.  If not, see <https://www.gnu.org/licenses/>. || #
# ||                                                                        || #
# **************************************************************************** #

# ||--------
# || USAGE
# ||--------

# For a general overview, see: README.md
#
# For available input options and relevant details, see:
# Section 2, 'OPTIONS AND INPUT DATA' below.
#
# Additional documentation can be found in docstrings in
# section 3, 'FUNCTIONS'.
#
# Install required module 'base36' using:
#    pip install base36

# **********--------------------------------------------------------********** #
# |                           ~~ [1] IMPORTS ~~                              | #
# **********--------------------------------------------------------********** #

# Import built-ins
import sys
import os
import os.path
import decimal
import operator
import sqlite3

from contextlib import closing

if sys.version_info[0] is 2:
    from collections import OrderedDict

# Import third party modules / packages
import base36


# **********--------------------------------------------------------********** #
# |                    ~~ [2] OPTIONS AND INPUT DATA ~~                      | #
# **********--------------------------------------------------------********** #

OUTPUT_PATH = "out/"
# PATH TO OUTPUT FOLDER
# All output files will be generated into this directory.

DB_PATH = 'db/MCSF_DB_2021-05_v6.0.db'
# PATH TO SQLITE3 BMRB RNA DATABASE FILE
# Contains chemical shift data and sequence data, among other things.
# Data in MCSF_DB_2021-05_v6.0.db was extracted from NMR-STAR-files
# from the BMRB.

FN_MATCH_LIST = "in/adit_nmr_matched_pdb_bmrb_entry_ids_210419.csv"
# FILE NAME OF BMRB-PDB MATCH LIST
# Contains a mapping between BMRB entries and PDB entries.
# An updated match table can be found at:
#  https://bmrb.io/ftp/pub/bmrb/nmr_pdb_integrated_data/adit_nmr_matched_pdb_bmrb_entry_ids.csv

XY_PAIR_DATA = [
    # X   [Y]  nuc_A  nuc_B   X-Y pair list file name
    ('A', 'U', 'H3',  'N3', 'in/FRABASE_A[U]_pairs.txt'),
    ('G', 'U', 'H3',  'N3', 'in/FRABASE_G[U]_pairs.txt'),
    ('U', 'G', 'H1',  'N1', 'in/FRABASE_U[G]_pairs.txt'),
    ('C', 'G', 'H1',  'N1', 'in/FRABASE_C[G]_pairs.txt'),
]
# FRABASE PAIR DATA
#
# X and Y are paired residues forming the nt pair X-Y.
# N.B: Chemical shift values will be extracted from nuclei nuc_A and
# nuc_B of residue Y (no CS will be extracted from residue X!).
#
# A FRABASE pair list can be generated using the 'Base pair' query at:
#  https://rnafrabase.cs.put.poznan.pl/?act=Structural%20Elements
#
# '1st residue': 'X'
# '2nd residue': 'Y'
# 'Experimental method': 'NMR'
# 'Base pair classification': Up to you ('Any' used in the MCSF article)
# 'Base-base parameters': Up to you ('Any' used in the MCSF article)
# 'Include all models of the structure': Un-checked
#
# Multiple pairlists can be inputted at once following the preconfigured
# example above (A-U and G-U with H3-N3, U-G and C-G with H1-N1).

AVERAGE_CS = False
# PERFORM OPTIONAL AVERAGING OF CS VALUES.
# Set to True to average CS readings of each nuclei. I.e. in cases
# where the same nuclei is probed under varying experimental
# conditions (pH, salts, time, etc) those shifts are combined into one
# average value. Set to 'False' in the MCSF article.

VERBOSE = True
# PRINT INFO TO HOST.
# Set to True to print information to console. Might increase run time
# slightly. Prints detailed info about CS averagings if AVERAGE_CS is
# set to True.


# **********--------------------------------------------------------********** #
# |                           ~~ [3] FUNCTIONS ~~                            | #
# **********--------------------------------------------------------********** #


def extract_data(file_path, sep=" ", clean=None, ignore_header=True):
    """Extract data from file as list of lists
    
    Extracts data from file, return data as a list of lines. Each line
    is split into a sublist according to a separator. A cleaning
    function can be passed to clean each line element.
    
    Parameters
    ----------
    file_path : str
        Path to file.
        
    sep : str
        Separator to separate line(s).
        
    clean : function
        A function to perform optional cleaning of line elements.
        
    ignore_header : bool
        If True, skip header comments (lines starting with '#') in
        input file.
        
    Returns
    -------
    d : list[list[str]]
    """
    d = list()
    if clean:
        with open(normpath(file_path)) as f:
            for line in f.readlines():
                if sep is " ":
                    line_arr = line.split()
                else:
                    line_arr = line.split(sep)
                line_arr = [clean(e) for e in line_arr]
                if ignore_header:
                    if line_arr[0].startswith('#'):
                        continue
                d.append(line_arr)
    else:
        with open(normpath(file_path)) as f:
            for line in f.readlines():
                if sep is " ":
                    line_arr = line.split()
                else:
                    line_arr = line.split(sep)
                if ignore_header:
                    if line_arr[0].startswith('#'):
                        continue
                d.append(line_arr)
    return d


def extract_data_lines(file_path, ignore_header=True):
    """Extract data from file as list of lines (strings)
    
    Extracts data from file, return data as a list of lines. Each line
    is extracted as a string.
    
    Parameters
    ----------
    file_path : str
        Path to file.
        
    ignore_header : bool
        If True, skip header comments (lines starting with '#') in
        input file.
        
    Returns
    -------
    d : list[str]
    """
    d = list()
    with open(normpath(file_path)) as f:
        for line in f.readlines():
            if ignore_header:
                if line.startswith('#'):
                    continue
            line = line.split(" ", 1)[1]
            d.append(line)
    return d


def write_data(out_path, data, sep=" "):
    """Write data to file
    
    Takes in a list of lists. Each element of outer list is writen to
    file as a separate line, with sublist elements separated by sep.
    
    Parameters
    ----------
    out_path : str
        Path to output file.
        
    data : iterable[iterable]
        Data to write to file.
        
    sep : str
        Separator to separate line elements. Default is space.
        
    Returns
    -------
    None : NoneType
    """
    with open(normpath(out_path), "w") as f:
        lines = (sep.join(line) + '\n' for line in data)
        f.writelines(lines)


def max_pdb_id(pdb_ids):
    """Choose maximum PDB ID from iterable collection of PDB ID strings
    
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


def match_frabase_with_bmrb(bmrb_pdb_match_list, frabase_pair_data):
    """Match FRABASE pairs with BMRB entries using BMRB-PDB match table
    
    In cases where one BMRB entry has multiple associated PDB IDs or in
    cases where one PDB ID has multiple BMRB IDs, reduce to a 1:1
    mapping where the newest deposited PDB and BMRB entries are used.
    
    WARNING: This function assumes that the newest BMRB and PDB entries
    have the highest ID number. This assumption is correct in most
    cases, but is not guaranteed to be true in all cases. For more info
    see docs of subroutines reduce_match_list() and max_pdb_id().
    
    Flow and connections
    --------------------
    A) BMRB <--(match_list)--> PDB <-- FRABASE (pair_data; has PDB ID)
         |
       match_frabase_with_bmrb(match_list, pair_data)
         ↓
    B) matched_pair_data (FRABASE data, BMRB ID, and PDB ID)
                                  where BMRB <-(1:1)-> PDB
        
    Parameters
    ----------
    bmrb_pdb_match_list : list[list[str, str]]
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
    
    frabase_pair_data : list
        FRABASE base pair list.
        Can be generated using the FRABASE 'Base pair' query at:
            rnafrabase.cs.put.poznan.pl/?act=Structural%20Elements
        
    Returns
    -------
    matched_frabase_pair_data : list
        The rows in frabase_pair_data which have a matching BMRB entry.
    """
    # Create a new match list with 1:1 mapping between PDB and BMRB IDs
    reduced_bmrb_pdb_match_dict = reduce_match_list(bmrb_pdb_match_list)
    reduced_pdb_bmrb_match_dict = {v: k for k, v  # reverse dict
                                   in reduced_bmrb_pdb_match_dict.items()}
    
    # Extract PDB IDs from reduced match list
    pdb_ids = reduced_bmrb_pdb_match_dict.values()
    
    # Extract FRABASE pair data if stated PDB ID in reduced match list
    matched_frabase_pair_data = []
    for frabase_pair in frabase_pair_data:
        pair_pdb_id = frabase_pair[1]
        if pair_pdb_id in pdb_ids:
            matched_frabase_pair_data.append(
                frabase_pair + [reduced_pdb_bmrb_match_dict[pair_pdb_id]])
    
    return matched_frabase_pair_data


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
    bmrb_pdb_match_list : list[list[str, str]]
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
    
    # Create new match list with unique PDB IDs (pick newest ones)
    pdb_reduced_bmrb_pdb_match_dict = reduce_pdb_ids(bmrb_pdb_match_list)

    # Create new match list with unique BMRB IDs (pick newest ones)
    bmrb_reduced_pdb_bmrb_match_dict = reduce_bmrb_ids(bmrb_pdb_match_list)
    
    # Combine into a 1:1 mapping
    reduced_match_dict = {}
    for pdb, bmrb in bmrb_reduced_pdb_bmrb_match_dict.items():
        if pdb in pdb_reduced_bmrb_pdb_match_dict.values():
            reduced_match_dict[bmrb] = pdb
    
    return reduced_match_dict


def reduce_pdb_ids(match_list):
    """If more than one PDB ID per BMRB ID, choose most recent PDB ID
    """
    # Collect PDB IDs for each BMRB ID
    bmrb_pdb_match_dict = dict()
    for _bmrb_id, _pdb_id in match_list:
        if _bmrb_id in bmrb_pdb_match_dict:
            bmrb_pdb_match_dict[_bmrb_id].append(_pdb_id)
        else:
            bmrb_pdb_match_dict[_bmrb_id] = [_pdb_id]
    
    # Create new match list with unique PDB IDs (pick highest ID)
    bmrb_pdb_match_dict = {
        bmrb_id: (max_pdb_id(pdb_ids) if (len(pdb_ids) > 1) else pdb_ids[0])
        for bmrb_id, pdb_ids in bmrb_pdb_match_dict.items()
    }
    
    return bmrb_pdb_match_dict


def reduce_bmrb_ids(match_list):
    """If more than one BMRB ID per PDB ID, choose most recent BMRB ID
    """
    # Collect BMRB IDs for each PDB ID
    pdb_bmrb_match_dict = dict()
    for _bmrb_id, _pdb_id in match_list:
        if _pdb_id in pdb_bmrb_match_dict:
            pdb_bmrb_match_dict[_pdb_id].append(_bmrb_id)
        else:
            pdb_bmrb_match_dict[_pdb_id] = [_bmrb_id]

    # Create new match list with unique BMRB IDs (pick highest ID)
    pdb_bmrb_match_dict = {
        bmrb_id: (max(pdb_ids) if (len(pdb_ids) > 1) else pdb_ids[0])
        for bmrb_id, pdb_ids in pdb_bmrb_match_dict.items()
    }
    
    return pdb_bmrb_match_dict


def average_shift_table(shift_table, verbose=False):
    """Averages CS values for each identical nuclei in a BMRB entry
    
    Takes in a list of tuples (i.e. a list of 'Shift' table rows)
    created from a "SELECT * FROM Shift WHERE entry_id=(?)" SQLite DB
    query (or more stringent; all rows needs to be from the same BMRB
    entry). Output is a list of tuples where each individual nuclei
    has its CS values and errors averaged.
    
    Two nuclei in the list are deemed the same if the shift_list_id is
    different but all of the following attributes are identical:
      bmrb_id
      assembly_id
      entity_id
      entity_assembly_id
      residue_id
      nt_type
      nuc_type

    If identical nuclei are found, their CS values and errors are
    averaged (significant digits accounted for).
    
    WARNING: If ambiguity codes in the set are not all the same, an
    ambiguity code of '-1' will be given.

    Parameters
    ----------
    shift_table : list[tuple[int, int, int, int, int,
                           int, str, str, int, int, int]]
        A list of tuples (rows) created from an SQL query to the
        'Shift' table.
        
    verbose : bool
        Set to True to print detailed information about averagings.
    
    Returns
    -------
    cs_table_averaged : list[tuple[int, int, int, int, int,
                           int, str, str, int, int, int]]
        Averaged list of 'Shift' table rows.
    """
    if not shift_table or len(shift_table) < 2:
        return shift_table
    
    # Sort together "identical" nuclei
    shift_table = sorted(shift_table,
                         key=operator.itemgetter(0, 2, 3, 4, 5, 6, 7))
    
    if verbose:
        global AVG_COUNTER
        AVG_COUNTER += 1
        print('\ni: {}, BMRB: {} \nBefore:'
              ''.format(AVG_COUNTER, shift_table[0][0]))
        print(shift_table)
    
    # Group up identical nuclei and perform averaging
    last_nuc = 0
    cs_table_averaged = []
    cs_data_to_avg = []
    for cs_row in shift_table:
        current_nuc = [
            cs_row[0],  # bmrb_id
            cs_row[2],  # assembly_id
            cs_row[3],  # entity_id
            cs_row[4],  # entity_assembly_id
            cs_row[5],  # residue_id
            cs_row[6],  # nt_type
            cs_row[7],  # nuc_type
        ]
        
        # If nuc identical to last nuc, add to list
        if current_nuc == last_nuc:
            cs_data_to_avg.append(cs_row)
        else:
            # If not identical
            if cs_data_to_avg:
                # and list with previous identical nuclei is not empty,
                # perform averaging on previous nuclei:
                cs_table_averaged.append(avg_cs(cs_data_to_avg))
                # Then begin to collect identical nuclei again again
                cs_data_to_avg = []
            # Add current nuclei to new list
            cs_data_to_avg.append(cs_row)
        
        last_nuc = current_nuc

    cs_table_averaged.append(avg_cs(cs_data_to_avg))
    
    # Re-sort output
    cs_table_averaged = sorted(cs_table_averaged,
                               key=operator.itemgetter(0, 1, 6, 7, 2, 3, 4, 5))
    # Print info
    if verbose:
        print('After:')
        print(cs_table_averaged)
    
    return cs_table_averaged


def avg_cs(cs_to_avg):
    """Average CS values and errors
    
    WARNING: If ambiguity codes in the set are not all the same, an
    ambiguity code of '-1' will be given.

    Parameters
    ----------
    cs_to_avg : list[tuple[int, int, int, int, int,
                           int, str, str, int, int, int]]
        Rows from Shift table to perform averaging on.
        
    Returns
    -------
    averaged : tuple[int, int, int, int, int,
                     int, str, str, int, int, int]
        Shift table row with averaged CS value and error.
    """
    shifts = []
    errs = []
    cs_decimals = []
    err_decimals = []
    ambiguity = []
    for cs_res_row in cs_to_avg:
        cs_val = cs_res_row[8]
        cs_val_err = cs_res_row[9]
        cs_amb_code = cs_res_row[10]
        
        # Extract number of decimals from data
        cs_val_dec = abs(
            decimal.Decimal(str(cs_val)).as_tuple().exponent
        )
        cs_val_err_dec = abs(
            decimal.Decimal(str(cs_val_err)).as_tuple().exponent
        )
        
        shifts.append(cs_val)
        errs.append(cs_val_err)
        cs_decimals.append(cs_val_dec)
        err_decimals.append(cs_val_err_dec)
        ambiguity.append(cs_amb_code)
    
    # Calculate number of decimals in resulting average
    cs_dec = min(cs_decimals)
    err_dec = min(err_decimals)
    
    # Calculate averages
    _avg_cs = round(sum(shifts) / len(shifts), cs_dec)
    _avg_err = round(sum(errs) / len(errs), err_dec)
    
    # Handle ambiguity codes
    ambiguity = list(set(ambiguity))
    if len(ambiguity) is not 1:
        ambiguity = -1
    else:
        ambiguity = ambiguity[0]

    # Create result vector
    averaged = list(cs_to_avg[0][:8]) + [_avg_cs, _avg_err, ambiguity]
    
    return averaged


def sqlite_query(cur, sql, args):
    """Perform SQLite3 query and return results
    """
    cur.execute(sql, args)
    return cur.fetchall()
    
    
def extract_nuclei_data_from_db(path_db, matched_pair_list,
                                nuc_type_A, nuc_type_B,
                                out_A_path, out_B_path,
                                average_cs, verbose):
    """Extract data from local BMRB RNA DB using matched pair list
    
    Extracts data on nuc_A and nuc_B in nt Y from local BMRB RNA DB
    using the matched (BMRB-PDB-FRABASE) pair list. Writes output to
    files.
    
    Two files are generated:
        1) 'out_X[y]_A.txt'
            Data on nuc_A (in nt Y) from local RNA BMRB DB.
        2) 'out_X[y]_B.txt'
            Data on nuc_B (in nt Y) from local RNA BMRB DB.
            
    Parameters
    ----------
    path_db : str
        Path to local BMRB RNA DB-file.
        
    matched_pair_list : list[list[str, ...]]
        Matched (BMRB-PDB-FRABASE) pair list.
    
    nuc_type_A : str
        Nucleus type (e.g. 'H1' or 'N3') of nucleus A in nt Y.
        
    nuc_type_B : str
        Nucleus type (e.g. 'H1' or 'N3') of nucleus B in nt Y.
        
    out_A_path : str
        Output path for nuc_A data-file.
    
    out_B_path : str
        Output path for nuc_B data-file.
        
    average_cs : bool
        Set to 'True' to average CS readings of each nuclei. I.e. in
        cases where the same nuclei is probed under varying experimental
        conditions (pH, salts, time, etc) those shifts are combined into
        one average value. Was set to 'False' in the MCSF article.
        
    verbose : bool
        Set to 'True' to print information to console. Prints detailed
        info about CS averagings if AVERAGE_CS is set to True.
    
    Returns
    -------
    None : NoneType
    """
    # Create output files and connect to DB file
    with open(out_A_path, "w+") as out_file_nuc_A, \
            open(out_B_path, "w+") as out_file_nuc_B:
        with closing(sqlite3.connect(path_db)) as conn:
            
            # Search through matchlist for pairs which have data in DB
            i_A = 1
            i_B = 1
            for XY_pair in matched_pair_list:
                
                # Assemble identifying info of nuc_A and nuc_B
                nt_Y_A_info = (
                    XY_pair[26],  # entry_id
                    XY_pair[12],  # residue_id
                    XY_pair[9],   # nt_type
                    nuc_type_A    # nuc_type
                )
                nt_Y_B_info = nt_Y_A_info[:3] + (nuc_type_B,)
                
                # Extract nuc_A and nuc_B data from DB if it exists
                sql = ("SELECT * "
                       "FROM Shift "
                       "WHERE ("
                       "    entry_id = ? "
                       "    AND residue_id = ? "
                       "    AND nt_type = ? "
                       "    AND nuc_type = ?"
                       ")")
                
                with closing(conn.cursor()) as cur:
                    db_rows_Y_A = sqlite_query(cur, sql, nt_Y_A_info)
                    db_rows_Y_B = sqlite_query(cur, sql, nt_Y_B_info)

                # Perform optional averaging of CS values
                if average_cs:
                    db_rows_Y_A = average_shift_table(db_rows_Y_A,
                                                      verbose=verbose)
                    db_rows_Y_B = average_shift_table(db_rows_Y_B,
                                                      verbose=verbose)
                
                # Write nuc_A data to file "out_X[y]_A.txt"
                for row_A in db_rows_Y_A:
                    
                    row_A = [str(e).strip() for e in row_A]
                    
                    data_XY_A = ' '.join([
                        str(i_A),     # idx
                        XY_pair[1],   # pdb_id
                        row_A[0],     # entry_id
                        row_A[1],     # shift_list_id
                        row_A[2],     # assembly_id
                        row_A[3],     # entity_id
                        row_A[4],     # entity_assembly_id
                        
                        # nt data
                        row_A[5],     # residue_id
                        row_A[6],     # nt_type
                        XY_pair[10],  # dbn
                        
                        # nt bp partner data
                        XY_pair[7],   # residue_id_bp_partner
                        XY_pair[4],   # nt_type_bp_partner
                        XY_pair[5],   # dbn_bp_partner
                        
                        # nuc data
                        row_A[7],     # nuc_type
                        row_A[8],     # cs_val
                        row_A[9],     # cs_val_err
                        row_A[10],    # ambiguity_code
                    ])
                    
                    out_file_nuc_A.write(data_XY_A)
                    out_file_nuc_A.write('\n')
                    i_A += 1
                    
                # Write nuc_A data to file "out_X[y]_A.txt"
                for row_B in db_rows_Y_B:
                    row_B = [str(e).strip() for e in row_B]
    
                    data_XY_B = " ".join([
                        str(i_B),     # idx
                        XY_pair[1],   # pdb_id
                        row_B[0],     # entry_id
                        row_B[1],     # shift_list_id
                        row_B[2],     # assembly_id
                        row_B[3],     # entity_id
                        row_B[4],     # entity_assembly_id
    
                        # nt data
                        row_B[5],     # residue_id
                        row_B[6],     # nt_type
                        XY_pair[10],  # dbn
    
                        # nt bp partner data
                        XY_pair[7],   # residue_id_bp_partner
                        XY_pair[4],   # nt_type_bp_partner
                        XY_pair[5],   # dbn_bp_partner
    
                        # nuc data
                        row_B[7],     # nuc_type
                        row_B[8],     # cs_val
                        row_B[9],     # cs_val_err
                        row_B[10],    # ambiguity_code
                    ])
                    
                    out_file_nuc_B.write(data_XY_B)
                    out_file_nuc_B.write('\n')
                    i_B += 1


def remove_duplicate_cs_vals(merged_AB):
    """Remove duplicate CS values from list of Y (XY bp) A-B nuc data
    
    Since it is possible for a nucleotide to have more than one
    base-pairing partner we might get the same CS value multiple times
    in our output.

    To prevent this from happening in our out_X[Y]_A-B.txt-file, we
    explicitly look for and remove data points that are duplicates
    in terms of CS values.

    WARNING:
    It is therefore possible for there to exist more bp partners
    than just the one specified in the output file (i.e. in cases
    where one nt has more than one pb partners).
    As it stands, the bp partner specified in the output file is
    the one that occurs first in the FRABASE pair list.
    
    Parameters
    ----------
    merged_AB : Iterable[Tuple[str]]
        A list of Y (XY bp) A-B nuc data.
    
    Returns
    -------
    merged_AB : Iterable[Tuple[str]]
        A list of Y (XY bp) A-B nuc data, with guaranteed unique CS vals
        for each entry.
    """
    merged_AB_copy = merged_AB[:]
    for AB_link_i in merged_AB:
        
        # Identifying information of nt Y
        Y_nt_data = [
            AB_link_i[1],  # pdb_id
            AB_link_i[2],  # entry_id
            AB_link_i[3],  # shift_list_id
            AB_link_i[4],  # assembly_id
            AB_link_i[5],  # entity_id
            AB_link_i[6],  # entity_assembly_id
            AB_link_i[7],  # residue_id
            AB_link_i[8],  # nt_type
        ]
        
        # List of all nt Y nuc_A-nuc_B data that matches Y_nt_data
        matching_AB_links = [AB_link for AB_link in merged_AB
                             if list(AB_link[1:9]) == Y_nt_data]
        
        # If CS values are identical, keep only first one
        sentinel = True
        for AB_link_j in matching_AB_links:
            cs_A, cs_B = AB_link_j[14], AB_link_j[18]
            if (AB_link_i[14] == cs_A) or (AB_link_i[18] == cs_B):
                if sentinel:
                    sentinel = False
                    continue
                try:
                    merged_AB_copy.remove(AB_link_j)
                except ValueError:
                    pass
                
    merged_AB = merged_AB_copy
    
    return merged_AB


def extract_from_pair_list(frabase_pair_list_data, path_db, path_output,
                           path_match_list, average_cs, verbose):
    """Main subroutine of program
    
    (a) Reduces the supplied BMRB-PDB match list into a 1-to-1 mapping.
    (b) Extracts relevant data from a supplied FRABASE pair list using
        the BMRB-PDB match list (i.e. extract data from pair list if
        PDB ID in pair list is also found in BMRB-PDB match list).
    (c) Condense data and make a connection between FRABASE structural
        data from pair list, PDB ID and BMRB ID.
    (d) Use this data to extract relevant data (CS values and more) from
        supplied local SQLite3 BMRB RNA DB.
    (e) Write resulting data to files. These files will later be merged
        into a single file that can be used as input for the plotting
        script '2_plot_cs_hsqc.py' data visualization purposes.
        
    Parameters
    ----------
    frabase_pair_list_data : Tuple[str, str, str, str, str]
        FRABASE pair list data including:
        (a) the pairing nucleotides (X-Y),
        (b) the nuclei of interest (nuc_A, nuc_B) in nt Y, as well as
        (c) the file name of the FRABASE pair list.
        
        Example:
            ('A', 'U', 'H3', 'N3', 'in/FRABASE_A[U]_pairs.txt')
    
    path_db : str
        Path to the local SQLite3 BMRB RNA DB.
        
    path_output : str
        Path to output folder. Output files will be written to this
        location.
        3 output files will be generated:
            1) 'out_X[Y]_A.txt'
                FRABASE, BMRB and PDB data for nuc_A in nt Y in pair XY.
            2) 'out_X[Y]_B.txt'
                FRABASE, BMRB and PDB data for nuc_B in nt Y in pair XY.
            3) 'out_X[Y]_A-B.txt'
                Combined data for nuc_A and nuc_B.
                Note: some data in files 1 and 2 might not be present in
                this combined file, since some nuc_A data might not have
                corresponding nuc_B data and vice versa.
        
    path_match_list : str
        Path to BMRB-PDB match list.
        
    average_cs : bool
        Set to 'True' to average CS readings of each nuclei. I.e. in
        cases where the same nuclei is probed under varying experimental
        conditions (pH, salts, time, etc) those shifts are combined into
        one average value. Was set to 'False' in the MCSF article.
        
    verbose : bool
        Set to 'True' to print information to console. Prints detailed
        info about CS averagings if 'average_cs' is set to 'True'.
        
    Returns
    -------
    None : NoneType
    """
    nt_X, nt_Y, nuc_A, nuc_B, path_pair_list = frabase_pair_list_data
    
    # Platform normalize paths
    path_output = normpath(path_output)
    path_db = normpath(path_db)
    path_pair_list = normpath(path_pair_list)
    path_match_list = normpath(path_match_list)
    
    # Construct file names and paths
    fn_A = 'out_{}[{}]_{}.txt'.format(nt_X, nt_Y, nuc_A)
    fn_B = 'out_{}[{}]_{}.txt'.format(nt_X, nt_Y, nuc_B)
    fn_all = 'out_{}[{}]_{}-{}.txt'.format(nt_X, nt_Y, nuc_A, nuc_B)
    
    out_nuc_A_path = join(path_output, fn_A)
    out_nuc_B_path = join(path_output, fn_B)
    out_nuc_AB_path = join(path_output, fn_all)
    
    # Load pair list and match list
    pair_list = extract_data(path_pair_list,
                             clean=lambda x: x.strip(),
                             sep=';')
    
    match_list = extract_data(path_match_list,
                              clean=lambda x: x.strip(),
                              sep=',')
    
    # Create a pair list with a 1:1 mapping between BMRB and PDB IDs
    matched_pair_list = match_frabase_with_bmrb(match_list,
                                                pair_list)

    # Go through matched pair list and extract data from DB if it exists
    extract_nuclei_data_from_db(path_db, matched_pair_list,
                                nuc_A, nuc_B,
                                out_nuc_A_path, out_nuc_B_path,
                                average_cs, verbose)

    # Extract data as list of lines (strings)
    nuc_A_data_arr = extract_data_lines(out_nuc_A_path)
    nuc_B_data_arr = extract_data_lines(out_nuc_B_path)
    
    # Remove duplicate lines (and keep order)
    if sys.version_info[0] is 2:
        nuc_A_data_arr = list(OrderedDict.fromkeys(nuc_A_data_arr))
        nuc_B_data_arr = list(OrderedDict.fromkeys(nuc_B_data_arr))
    else:
        nuc_A_data_arr = list(dict.fromkeys(nuc_A_data_arr))
        nuc_B_data_arr = list(dict.fromkeys(nuc_B_data_arr))
        
    # Note: Duplicates are a result of multi-paired bases as specified
    # in the FRABASE pair list. I.e. sometimes one nucleotide
    # interacts/pairs to multiple other nucleotides.

    # Cleaning of data
    nuc_A_data_arr = [A_d.strip('\n').split(" ")
                      for A_d in nuc_A_data_arr]
    nuc_B_data_arr = [B_d.strip('\n').split(" ")
                      for B_d in nuc_B_data_arr]
    
    # Combine data from nuc_A and nuc_b (of nt Y) into one output file
    merged_AB = []
    i = 1
    for nuc_A_data in nuc_A_data_arr:
        
        A_nt = (
            nuc_A_data[0],  # pdb_id
            nuc_A_data[1],  # entry_id
            nuc_A_data[2],  # shift_list_id
            nuc_A_data[3],  # assembly_id
            nuc_A_data[4],  # entity_id
            nuc_A_data[5],  # entity_assembly_id
            nuc_A_data[6],  # residue_id
            nuc_A_data[7],  # nt_type
        )
        
        for nuc_B_data in nuc_B_data_arr:
            B_nt = (
                nuc_B_data[0],  # pdb_id
                nuc_B_data[1],  # entry_id
                nuc_B_data[2],  # shift_list_id
                nuc_B_data[3],  # assembly_id
                nuc_B_data[4],  # entity_id
                nuc_B_data[5],  # entity_assembly_id
                nuc_B_data[6],  # residue_id
                nuc_B_data[7],  # nt_type
            )
            
            if A_nt == B_nt:
                matched_AB_data = (
                    # bp X-Y data:
                    str(i),               # i
                    str(nuc_A_data[0]),   # pdb_id
                    str(nuc_A_data[1]),   # entry_id
                    str(nuc_A_data[2]),   # shift_list_id
                    str(nuc_A_data[3]),   # assembly_id
                    str(nuc_A_data[4]),   # entity_id
                    str(nuc_A_data[5]),   # entity_assembly_id
                    str(nuc_A_data[6]),   # residue_id
                    str(nuc_A_data[7]),   # nt_type
                    str(nuc_A_data[8]),   # nt_dbn
                    str(nuc_A_data[9]),   # residue_id_bp_partner
                    str(nuc_A_data[10]),  # nt_type_bp_partner
                    str(nuc_A_data[11]),  # nt_dbn_bp_partner

                    # nt Y, nuc A data:
                    str(nuc_A_data[12]),  # nuc_type_h
                    str(nuc_A_data[13]),  # cs_data_h
                    str(nuc_A_data[14]),  # cs_err_data_h
                    str(nuc_A_data[15]),  # cs_ambiguity_code_h

                    # nt Y, nuc B data:
                    str(nuc_B_data[12]),  # nuc_type_n
                    str(nuc_B_data[13]),  # cs_data_n
                    str(nuc_B_data[14]),  # cs_err_data_n
                    str(nuc_B_data[15]),  # cs_ambiguity_code_n
                )
                merged_AB.append(matched_AB_data)
                
                i += 1
    
    # Remove duplicate CS values due to multi-paired primary nts
    merged_AB = remove_duplicate_cs_vals(merged_AB)
    
    write_data(out_nuc_AB_path, merged_AB, sep=' ')


def merge_files(pair_data, output_path):
    """Merge 'out_X[Y]_A-B.txt' output files into one
    """
    names = []
    files_to_combine = []
    for pair_list in pair_data:
        nt_X, nt_Y, nuc_A, nuc_B, path_pair_list = pair_list
        name = "{}[{}]_{}-{}".format(nt_X, nt_Y, nuc_A, nuc_B)
        files_to_combine.append("out_{}.txt".format(name))
        names.append(name)
        
    data_from_files = []
    for file_name in files_to_combine:
        file_path = join(output_path, file_name)
        data_from_files.append(
            extract_data(file_path, sep=" ")
        )

    data_all = []
    for data in data_from_files:
        data_all += data

    output_file_name = "out_ALL_({}).txt".format(")_(".join(names))
    output_file_path = join(output_path, output_file_name)
    write_data(output_file_path, data_all, sep=" ")
    
    
def normpath(path):
    """Normalize path
    """
    return os.path.normpath(path)


def join(root, branch):
    """Join and normalize paths
    """
    return os.path.normpath(os.path.join(os.path.normpath(root),
                                         os.path.normpath(branch)))


# **********--------------------------------------------------------********** #
# |                              ~~ [4] main ~~                              | #
# **********--------------------------------------------------------********** #

def main(db_path, output_path, fn_match_list, pair_data, average_cs, verbose):

    if verbose:
        import time
        start_time = time.time()
    else:
        time = None
        start_time = 0
    
    # Run main subroutine for each FRABASE pair list
    for _pair_list_data in pair_data:
        extract_from_pair_list(_pair_list_data, db_path, output_path,
                               fn_match_list, average_cs, verbose)
    # Combine all out_X[Y]_A-B.txt files
    merge_files(pair_data, output_path)
    
    if verbose:
        print("\nDone!\n")
        print("--- %s seconds ---\n" % (time.time() - start_time))
        

# **********--------------------------------------------------------********** #
# |                              ~~ [5] RUN ~~                               | #
# **********--------------------------------------------------------********** #

if VERBOSE:
    if AVERAGE_CS:
        # AVG_COUNTER counts the number of CS averagings.
        # For more info, see function average_cs_table().
        AVG_COUNTER = 0

main(db_path=DB_PATH, output_path=OUTPUT_PATH, fn_match_list=FN_MATCH_LIST,
     pair_data=XY_PAIR_DATA, average_cs=AVERAGE_CS, verbose=VERBOSE)
