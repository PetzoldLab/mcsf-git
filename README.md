# mcsf 

1. 1_bmrb2pdb.py
	1. Download BMRB NMR star-file for hairpins containing nitrogen CS from https://bmrb.io/search/query_grid/query_1_20.html. An example is called bmrb_entries_15N_RNA_20
	2. Put bmrb_entries_15N_RNA_20 folder in the “main” folder
	3. Use 1_bmrb2pdb.py script to fish out PDB-IDs from bmrb_entries_15N_RNA_20
		1. This script prints out PDB and BMRB matched IDs which are to be put into a text file named for example “bmrb_pdb.txt”
		2. Put that text file in the “main” folder

2. 2_retrieve_sec_from_fra.py 
	1. Input: “bmrb_pdb.txt”
	2. Use 2_retrieve_sec_from_fra.py  to retrieve PDB ID, sequence and dot bracket structures from FRABASE, for the PDB and BMRB matched IDs.
	3. Manually remove structures containing: [, }, -, signs indicating pseudoknot structures, e.g 	((((((((.....[[[[[[.)))))))).........((((....))))..]]]]]]...
	4. Put the structures in a text file called “structures.txt” in the “lists” folder

3. 3_extract.py
	1. Install the python module “forgi” 1.1
	2. In case the forgi module gives errors such as: “ImportError: No module called builtins found”, then comment out the built-ins in the following scripts:
		Residue.py
		Mcannotate.py
		stuff.py
		bulge_graph.py
	3. Input: structures.txt, bmrb_pdb.txt
	4. Output: nuclei.txt

4. 4_reorganize_n_calib_2020_imino.py
	1. Input: nuclei.txt
	2. Output: nuc_sparse.txt

5. 5_plot_bmrb_cs_imino.py
	1. Input: nuc_sparse.txt
	2. Output: chemical shift plot
