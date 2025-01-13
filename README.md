# PV-BIOGRID EVIDENCE MAPPING
## Python virtual env

	python3 -m venv pv

	source pv/bin/env

	export DATA_DIR=<PATH_TO_DATA DIR>

Install necessary python modules

	pip3 install matplotlib matplotlib_venn psycopg2-binary

## Requirements:

1. *Required, must exist*: af2-models-split-ifresid_.tsv: 310582 rows, pdockq scores but no iptm scores
2. *Required, must exist*: af2_iptm_pdockq.tsv: 103969 entries from the AF2 file with pdockq scores and iptm scores
3. *Required, must exist*: suppl_ppi_models_.tsv: 486100 entries
4. *Generated from Biogrid data processing*: updated_bg_source.tsv: 976325 entries of human interactions (final output from stage 3 of Biogrid data processing)

## BIOGRID DATA PROCESSING: 
At each step, an input file is read and output file is generated, which then acts as the input for the next stage.
All paths must be changed to whatever your system is using. Fully-qualified paths are preferred. If you are copy/paste-ing the examples from this README, ensure you set DATA_DIR appropriately!

## Step 1. Download data from Biogrid and extract human related interactions and format header.
Comment: output file created named "biogrid_human_interactions.txt" which has no spaces in header
```
./biogrid_fetch_hformat.sh
```

### Step 2. BG data processing (Stage 1)
Read step 1 output file, and extract uniprot and refseq ids as fields/columns

```
./biogrid_processing.py -s1 \
        -i $DATA_DIR/biogrid_human_interactions.txt \
        -o $DATA_DIR/updated_bg_human_interactions.tsv \
        --verbose
```
*Stage 1 completed in 32 seconds*

### Step 3. BG data processing (Stage 2)
Read step 2 output file, and get the missing ids uniprot ids from ProtVar Database

Files required:

*pv_db.ini*: config file containing the PV DB conn details for stage 2. This can be either specified on the command line with `--config-file` or loaded from the current directory by default.

**NOTE**: Connection to the DB is required for this stage to work

```
./biogrid_processing.py -s2 \
        --config-file ../pv_db.ini
        -i $DATA_DIR/updated_bg_human_interactions.tsv \
        -o $DATA_DIR/updated_bg_human_interactions_PVDB.tsv \
      	--verbose
```
*Stage 2 completed in 321.92 seconds*

### Step 4. BG data processing (Stage 3)
Read step 3 output file, and perform "merging of column values" for a given interaction id, since biogrid data has 1:many relationship for a given interaction_id.

```
    ./biogrid_processing.py -s3 \
        -i $DATA_DIR/updated_bg_human_interactions_PVDB.tsv \
        -o $DATA_DIR/updated_bg_source.tsv \
        --verbose
```
*Stage 3 completed in 28.50 seconds*

## AF2, SUPPL_PPI, BG DATA Merging: 
'''pv_data_merging.py''' is a script that takes as input file listing the filenames to merge.
Output file is generated which is used in the subsequent steps before finally generating a tsv 
file which contains biogrid data mapped to af2 and suppl_ppi data

Files required:

- *af2-models-split-ifresid_.tsv*
- *suppl_ppi_models_.tsv*
- *af2_iptm_pdockq.tsv*

### Recap: comparison of AF2 and suppl_ppi files
![af2_suppl](https://github.com/user-attachments/assets/232b490f-da2c-409b-8e85-ed0e5e20f31e)

### Step 1) Concatenate files: Files 2 and 3

![supp_af2_receieved](https://github.com/user-attachments/assets/e798a598-dd3d-4d46-a1b6-e20a080bcc2b)

Generate an updated suppl_ppi file which includes the missing af2 entries.

Required Data:

`af2_iptm_pdockq.tsv`: data containing 103969 entries that were previously missing
`suppl_ppi_models_.tsv`: file from PV

**Note** `--input_file_list` must be a textfile listing files to concatenate/merge. For this step, the list must be the two files mentioned above. DO NOT provide '--common_col' command-line argument.

Example:

```
/full/path/suppl_ppi_models_.tsv
/full/path/af2_iptm_pdockq.tsv
```

Running the script.

Example:
```
./pv_data_merging.py --input_file_list $DATA_DIR/input_files_suppl_mAF2.txt --outfile $DATA_DIR/updated_suppl_ppi.tsv --verbose 
```
*ELAPSED TIME: 4.60 seconds*

### Step 2) Merge files on "interaction_id": File 1 and combined result of (Files 2 and 3)
![af2_updated_suppl](https://github.com/user-attachments/assets/4116549a-4085-48be-beed-048b381c4363)

Note that, as above, `--input_file_list` must be a file containing a list of files (with full paths) to concatenate/merge.

		## af2-models-split-ifresid_.tsv file from PV
		## updated_suppl_ppi.tsv output of Step 1
		These two files, one on each line == INPUT for script 

	  Filenames:
	/full/path/
	/full/path/

```
 	./pv_data_merging.py --input_file_list $DATA_DIR/input_files_AF2_updatedSuppl.txt --common_col interaction_id --outfile $DATA_DIR/af2_suppl_ppi_combined.tsv --verbose 
```
*ELAPSED TIME: 62.67 seconds*

### Step 3) FINAL MERGE:  Merge files on "interaction_id": (combined result of Files 2, 3 and 1) and File 4
![af2_bg_combined](https://github.com/user-attachments/assets/e162e3de-300b-432c-bfd3-58557aed29fe)

Note that, as above, `--input_file_list` must be a file containing a list of files (with full paths) to concatenate/merge.

        ## af2_suppl_ppi_combined.tsv output of Step 2
	## updated_bg_source.tsv output stage 3 of biogrid processing (biogrid_processing.py -s3)
	These two files, one on each line == INPUT for script
 
	Filenames:
	/full/path/af2_suppl_ppi_combined.tsv
	/full/path/updated_bg_source.tsv
  
```
./pv_data_merging.py --input_file_list $DATA_DIR/input_files_updatedAF2_biogrid.txt --common_col interaction_id --outfile $DATA_DIR/af2_suppl_ppi_biogrid_combined.tsv --verbose 
```
*ELAPSED TIME: 32.69 seconds*
