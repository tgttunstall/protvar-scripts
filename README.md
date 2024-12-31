# BIOGRID DATA PROCESSING: 

At each step, an input file is read and output file is generated, which then acts as the input for the next stage.
Files required:
	pv_db.ini: config file containing the PV DB conn details for stage 2

## Step1. Download data from Biogrid and extract human related interactions and format header.
# comment: output file created named "biogrid_human_interactions.txt" which has no spaces in header
    ./biogrid_fetch_hformat.sh 

## Step 2. BG data processing (Stage 1): Read step 1 output file, and extract uniprot and refseq ids as fields/columns
# script: 
    ./biogrid_processing.py -s1 
        -i /home/pub/Work/data_arise_proteome/protvar/biogrid/biogrid_human_interactions.txt
        -o /home/pub/Work/data_arise_proteome/protvar/biogrid/updated_bg_human_interactions.tsv
        --verbose

## Step 3. BG data processing (Stage 2): Read step 2 output file, and get the missing ids uniprot ids from ProtVar Database      
# script: 
    ./biogrid_processing.py -s2 
        -i /home/pub/Work/data_arise_proteome/protvar/biogrid/updated_bg_human_interactions.tsv 
        -o /home/pub/Work/data_arise_proteome/protvar/biogrid/updated_bg_human_interactions_PVDB.tsv
      	-- verbose
# ==> Stage 2 completed in 321.92 seconds
     	
Step 4. BG data processing (Stage 3): Read step 3 output file, and perform "merging of column values" for a given interaction id
since biogrid data has 1:many relationship for a given interaction_id.
# script: 
    ./biogrid_processing.py -s3 
        -i /home/pub/Work/data_arise_proteome/protvar/biogrid/updated_bg_human_interactions_PVDB.tsv
        -o /home/pub/Work/data_arise_proteome/protvar/biogrid/updated_bg_source.tsv
        --verbose
# ==> Stage 3 completed in 28.50 seconds

##########################
# AF2, SUPPL_PPI, BG DATA Merging: 
pv_data_merging.py is a script that takes as input file listing the filenames to merge.
output file is generated which is used in  the subsequent steps before final generating a tsv 
file which is biogrid data mapped to af2 and suppl_ppi data
##########################
Files required:

- suppl_ppi_models_.tsv
- af2-models-split-ifresid_.tsv
- af2_iptm_pdockq.tsv



MISSION: Combine files
	1. af2-models-split-ifresid_.tsv: 310582 rows, pdockq scores but no iptm
	2. af2_iptm_pdockq.tsv: 103969 entries from the AF2 file ^^ that David provdied pdockq and iptm scores for.
	3. suppl_ppi_models_.tsv: 486100 entries
	4. updated_bg_source.tsv: 976325 entries of human interactions (final output from stage 3 of Biogrid data processing


# Step 1) Concatenate files: Files 2 and 3
	# Generate an updated suppl_ppi file which includes the missing af2 entries.
		## af2_iptm_pdockq.tsv: data received from David which are 103969 entries that were previously missing
		## suppl_ppi_models_.tsv: file from PV
		These two files, one on each line == INPUT for script 
		Output: 590069 rows
# DO NOT provide '--common_col' cmd arg so it does a simple concatenate
script: 
	./pv_data_merging.py --input_file_list /home/pub/Work/data_arise_proteome/protvar/input_files_suppl_mAF2.txt --outfile /home/pub/Work/data_arise_proteome/protvar/output/updated_suppl_ppi.tsv --verbose 

# ==> ELAPSED TIME: 4.60 seconds.

# Step 2) Merge files on "interaction_id": (Files 2 and 3) and File 1
		## af2-models-split-ifresid_.tsv: file from PV
		## updated_suppl_ppi.tsv: output of Step 1
		These two files, one on each line == INPUT for script 
script: 
	./pv_data_merging.py --input_file_list /home/pub/Work/data_arise_proteome/protvar/input_files_AF2_updatedSuppl.txt --common_col interaction_id --outfile /home/pub/Work/data_arise_proteome/protvar/output/af2_suppl_ppi_combined.tsv --verbose 

#==> ELAPSED TIME: 62.67 seconds.

# Step 3) FINAL MERGE:  Merge files on "interaction_id": (Files 2 and 3 and 1) and File 4
		## af2_suppl_ppi_combined: output of Step 2
		## updated_bg_source.tsv: output stage 3 of biogrid processing (biogrid_processing.py -s3)
		These two files, one on each line == INPUT for script 
script: 
	./pv_data_merging.py --input_file_list /home/pub/Work/data_arise_proteome/protvar/input_files_updatedAF2_biogrid.txt --common_col interaction_id --outfile /home/pub/Work/data_arise_proteome/protvar/output/af2_suppl_ppi_biogrid_combined.tsv --verbose 

# ==> ELAPSED TIME: 32.69 seconds.


