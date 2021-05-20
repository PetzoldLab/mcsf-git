# MCSF - BMRB chemical shift query

The MCSF scripts constitutes a pipeline to query a local BMRB SQLite3 RNA database for 
chemical shifts (CS) and plot them in a HSQC-type of plot, on the example of imino nitrogens and 
protons.

## Description
See the documentation within the script 1_extract_from_pair_list.py for detailed information.

### Summary:
The script `1_extract_from_pair_list.py` takes as input:

 - (a) a BMRB-PDB match list (which can be downloaded from the BMRB), 

- (b) a FRABASE pair list (which can be retrieved using a FRABASE `Structural elements > Base pair` 
query), 

- (c) the local BMRB CS database (supplied, `MCSF_DB_v6.0.db`). 

In some cases
one BMRB ID is linked to multiple PDB IDs and in some cases one BMRB ID is linked to multiple PDB 
IDs. In any of these cases, the script chooses the entry with the highest ID number which is assumed
the most recent one.
Both BMRB and PDB entries are generally sequentially chosen. There are some minor exceptions to 
this, however: in the past, authors were allowed to choose their own PDB IDs. In addition, 
in some edge cases BMRB IDs might not be perfectly sequential. 

The script takes the BMRB-PDB match list and reduces it to a 1-to-1 mapping between BMRB and PDB 
entries before proceeding with data extraction.

For each nucleotide (nt)-pair in the FRABASE pair list, the PDB ID found is used to fish out a BMRB ID 
(using the BMRB-PDB match list) which is then appended to the pair table.

The script then queries the local SQLite3 BMRB RNA database (generated with scripts 
written by Noah Hopkins) to extract chemical shift data of the imino proton and nitrogen of 
nucleotides of interest. 
The CSs of the proton and nitrogen within one nucleotide are then linked. 
To be linked they need to have identical: `pdb_id`, `entry_id`, `shift_list_id`, `assembly_id`, 
`entity_id`, `entity_assembly_id`, `residue_id`, and `nt_type`. 

**Note 1:** It is possible that one nucleotide base-pairs to multiple other nucleotides. If this is 
the case, it will be stated in the pair table from FRABASE. Since each row in the output file 
contains data for one nt in a specific base pair and its corresponding proton and nitrogen CS readings, having 
both base pair `X-Y` and `X-Z` in the output file will result in duplicate CS values (for the `X` nt). 
Therefore, to prevent duplicate CS readings in the output, we keep only the first occurring base pair as
stated in the FRABASE table. It is therefore possible for there to exist more bp partners than 
just the one specified in the output file.

**Note 2:** The NMR-star files, which are the basis of the data included in the local bmrb database, 
sometimes contain multiple chemical shift lists. This can be due to different sample conditions (pH,
temperature, MgCl2 concentration). In the article, we treat the data from the same nucleus, but from
different chemical shift lists, as separate data points. An option to average all CS values for each
nucleus, for entries with multiple shift lists, have been added in case it is useful to you. 

## Usage

1. Download the BMRB-PDB matched list from: 
   
   https://bmrb.io/ftp/pub/bmrb/nmr_pdb_integrated_data/adit_nmr_matched_pdb_bmrb_entry_ids.csv
   
	1. Convert the `.csv`-file to `.txt` and name it for example: `BMRB_match_PDB_<DATE>.txt`
	2. Put the text file in the lists folder

2. Generate a FRABASE pair list using the `Structural elements > Base Pair` query at: 
   
   http://rnafrabase.cs.put.poznan.pl/?act=Structural%20Elements 

	1. For imino information of H1-N1 of G in a G-C base-pair: Input `C` as the residue of the 1st 
	   strand and `G` as the residue of the 2nd strand. Experimental method:`NMR`. Include all 
	   base-pair types by leaving the Base pair classification to `Any`. Download the pair list.
	2. Repeat step 2.i for A-U, G-U and U-G base-pairs. Note that the scripts will extract chemical 
	   shift information of the base of the 2nd nt (Y in X-Y pair), hence the need for both a G-U 
	   and a U-G list.

3. `1_extract_from_pair_list.py`
	1. Run the script `1_extract_from_pair_list.py` with the text files generated in step 1 and 2 as
	   inputs
	2. Output: `out_ALL_(A[U]_H3-N3)_(G[U]_H3-N3)_(U[G]_H1-N1)_(C[G]_H1-N1).txt`

4. `2_plot_cs_hsqc.py`
	1. Input: `out_ALL_(A[U]_H3-N3)_(G[U]_H3-N3)_(U[G]_H1-N1)_(C[G]_H1-N1).txt`
	2. Output: chemical shift plot (in `.eps`-format)
	
# Authors and acknowledgement

Noah Hopkins has developed the scripts generating the local chemical shift database. Noah Hopkins 
has together with Magdalena Riad developed the script to extract chemical shift data from the local
database based on scripts written by Lorenzo Baronti and Hampus Karlsson. The work was done within 
the group of Katja Petzold at Karolinska Institutet.
