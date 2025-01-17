# PV-BIOGRID EVIDENCE MAPPING
## Python virtual env

	python3 -m venv pv

	source pv/bin/activate

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
        --config ../pv_db.ini \
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
    	--write_counts \
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
![venn_af2_suppl_ppi](https://github.com/user-attachments/assets/3a91b7c6-bdc0-4686-9e87-6ee39fcbb599)

### Step 1) Concatenate files: Files 2 and 3
![venn_mAF2_suppl_ppi](https://github.com/user-attachments/assets/44db58f1-0148-4068-a80d-019b23de7c84)

Generate an updated suppl_ppi file which includes the missing af2 entries.

Required Data:

`af2_iptm_pdockq.tsv`: data containing 103969 entries that were previously missing
`suppl_ppi_models_.tsv`: file from PV

**Note:** `--input_file_list` must be a textfile listing files to concatenate/merge. For this step, the list must be the two files mentioned above. **DO NOT** provide `--common_col` command-line argument.

Filenames:

```
/full/path/suppl_ppi_models_.tsv
/full/path/af2_iptm_pdockq.tsv
```

Running script:

```
./pv_data_merging.py --input_file_list $DATA_DIR/input_files_suppl_mAF2.txt \
     --outfile $DATA_DIR/updated_suppl_ppi.tsv \
     --verbose 
```
*ELAPSED TIME: 4.60 seconds*

### Step 2) Merge files on "interaction_id": File 1 and combined result of (Files 2 and 3)
![af2_updated_suppl](https://github.com/user-attachments/assets/4116549a-4085-48be-beed-048b381c4363)
![venn_updatedSuppl_AF2](https://github.com/user-attachments/assets/3f756c50-85c5-4927-8fa3-4bce3ba4126f)

**Note:** As above, `--input_file_list` must be a file containing a list of files (with full paths) to concatenate/merge.

	  	Filenames:
		/full/path/af2-models-split-ifresid_.tsv
		/full/path/updated_suppl_ppi.tsv
		
		## af2-models-split-ifresid_.tsv file from PV
        	## updated_suppl_ppi.tsv output of Step 1
        	These two files, one on each line == INPUT for script
        	
Running script:

```
./pv_data_merging.py --input_file_list $DATA_DIR/input_files_AF2_updatedSuppl.txt \
 	 --common_col interaction_id \
 	 --outfile $DATA_DIR/af2_suppl_ppi_combined.tsv \
 	 --verbose
```
*ELAPSED TIME: 62.67 seconds*

### Step 3) FINAL MERGE: Merge files on "interaction_id": (combined result of Files 2, 3 and 1) and File 4
![af2_bg_combined](https://github.com/user-attachments/assets/e162e3de-300b-432c-bfd3-58557aed29fe)

**Note:** Similarly, as above, `--input_file_list` must be a file containing a list of files (with full paths) to concatenate/merge.
 
		Filenames:
		/full/path/af2_suppl_ppi_combined.tsv
		/full/path/updated_bg_source.tsv	
		
        ## af2_suppl_ppi_combined.tsv output of Step 2
		## updated_bg_source.tsv output stage 3 of biogrid processing (biogrid_processing.py -s3)
		These two files, one on each line == INPUT for script
  
Running script:

```
./pv_data_merging.py \
     --input_file_list $DATA_DIR/input_files_updatedAF2_biogrid.txt \
     --common_col interaction_id \
     --outfile $DATA_DIR/af2_suppl_ppi_biogrid_combined.tsv \
     --verbose
```
*ELAPSED TIME: 32.69 seconds*

## PLOT VENN DIAGRAMS

This is sed to generate Venn Diagrms for *AF2, Supp_ppi, and Biogrid datasets*. These are based on *'interaction_id'* column found in all the datasets used.

### Compare AF2 and Suppl Data
```
./plot_venn.py --files $DATA_DIR/af2-models-split-ifresid_.tsv $DATA_DIR/suppl_ppi_models_.tsv \
    --key_column interaction_id \
    --labels AF2 Suppl_ppi \
    --plot_title "Comparison of AF2 and Suppl_ppi" \
    --output_file $DATA_DIR/venn_af2_suppl_ppi.png \
    --verbose
```

### Compare Suppl_ppi and Missing AF2

```
./plot_venn.py --files $DATA_DIR/af2_iptm_pdockq.tsv $DATA_DIR/suppl_ppi_models_.tsv \
    --key_column interaction_id \
    --labels M_AF2 Suppl_ppi \
    --plot_title "Comparison of M_AF2 and Suppl_ppi" \
    --output_file $DATA_DIR/venn_mAF2_suppl_ppi.png \
    --verbose
```
    
### Compare Updated Suppl_ppi and AF2

```
./plot_venn.py --files $DATA_DIR/updated_suppl_ppi.tsv $DATA_DIR/af2-models-split-ifresid_.tsv \
    --key_column interaction_id \
    --labels Updated_Suppl AF2  \
    --plot_title "Comparison of updated_Suppl and AF2" \
    --output_file $DATA_DIR/venn_updatedSuppl_AF2.png \
    --verbose
```
**Note:**
*updated_suppl_ppi.tsv*: result of combining the following files using *pv_data_merging.py*
- suppl_ppi_models_.tsv
- af2_iptm_pdockq.tsv

### Compare combined AF2_suppl_ppi and Biogrid

```
./plot_venn.py --files $DATA_DIR/af2_suppl_ppi_combined.tsv $DATA_DIR/updated_bg_source.tsv \
    --key_column interaction_id \
    --labels updated_AF2 BG \
    --plot_title "Comparison of updated_AF2 and Biogrid"" \
    --output_file $DATA_DIR/venn_updatedAF2_BG.png \
    --verbose
```
**Note:**
*af2_suppl_ppi_combined.tsv*: result of combining the following files using *pv_data_merging.py*
- af2-models-split-ifresid_.tsv
- updated_suppl_ppi.tsv

### Triple comparison: AF2, updated suppl_ppi, and Biogrid

```
./plot_venn.py --files $DATA_DIR/updated_suppl_ppi.tsv $DATA_DIR/af2-models-split-ifresid_.tsv $DATA_DIR/updated_bg_source.tsv \
    --key_column interaction_id \
    --labels updatedSuppl AF2 BG  \
    --plot_title "Comparison of updated_Suppl_ppi, AF2 and Biogrid" \
    --output_file $DATA_DIR/venn3_af2_suppl_bg.png \
    --verbose
```
