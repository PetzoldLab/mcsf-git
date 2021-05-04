# mcsf - plotting imino chemical shifts

1. Download the BMRB-PDB matched list from: https://bmrb.io/ftp/pub/bmrb/nmr_pdb_integrated_data/adit_nmr_matched_pdb_bmrb_entry_ids.csv
	1. Convert the .csv file to .txt and name it for example: "BMRB_match_PDB_210419.txt"
	2. Put the text file in the lists folder

2. 1_retrieve_sec_from_fra.py 
	1. Input: “BMRB_match_PDB_.txt”
	2. Use 1A_retrieve_dbn_from_fra.py  to retrieve PDB IDs, sequences and dot bracket structures from FRABASE, for the PDB and BMRB matched IDs.
	3. Manually remove structures containing: "-" signs indicating structures containing nucleotides lacking structural data, e.g 	 -------(.(((......))).)-----
	4. Put the structures in a text file called “structures_wo_dash_210419.txt” in the “lists” folder

3. 2A_extract_from_dbn.py
	1. Install the python module “forgi”, any version
	2. In case the forgi module gives errors such as: “ImportError: No module called builtins found”, then comment out the built-ins in the following scripts:
		Residue.py
		Mcannotate.py
		stuff.py
		bulge_graph.py
	3. Input: “structures_wo_dash_210419.txt”, “BMRB_match_PDB_.txt”
	4. Output: nuclei_imino_21.txt

4. 3AB_reorganize.py
	1. Input: nuclei_imino_21.txt
	2. Output: imino_nuc_sparse_21.txt

5. 4AB_plot_cs_hsqc.py
	1. Input: imino_nuc_sparse_21.txt
	2. Output: chemical shift plot
