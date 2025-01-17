#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 14:10:38 2024

@author: tanu
"""
###############################################################################
# Process interaction data from the BioGRID database through multiple stages,
# Each stage of the script serves a specific purpose:

# Stage 1: Extracts and processes UniProt and RefSeq IDs from the given input file, updating the data accordingly
# and writing the results to an output file specified by the user.

# Stage 2: Connects to a PostgreSQL database to fetch missing UniProt IDs based on RefSeq IDs present in the input file.
# The script updates these IDs in the data, optionally generating a tracking dict to capture details of the data modifications.

# Stage 3: Merges columns of data based on user-specified criteria, facilitating the consolidation of interaction data
# across different datasets. This stage also supports the generation of tracking information regarding the merge operations.

# The script is controlled via command-line arguments, enabling users to specify the stage, input and output files,
# database configurations, and other options such as verbose logging and tracking data output.

# Usage of this script requires a proper setup of the Python environment with necessary dependencies installed,
# and appropriate database access configurations for fetching UniProt IDs.
# Custom functions imported are clearly indicated in the imports

# Example usage:
#     python this_script.py -s 1 -i 'input_file.tsv' -o 'output_file.tsv' --verbose
###############################################################################
import os
import sys
import argparse
import configparser
import time

# my custom functions
#from my_pv import *
from my_pv.pv_functions import read_file_to_list_of_dicts, extract_ids, process_biogrid_file, dict_value_merge
from my_pv.pv_functions import write_list_of_dicts, write_nested_dict
###############################################################################
def setup_arguments():
    parser = argparse.ArgumentParser(description="Process BioGRID data in multiple stages.")
    
    parser.add_argument('-s', '--stage', 
                        type = int, choices = [1, 2, 3], 
                        required = True, 
                        help = "Stage to execute (1, 2, or 3)")
    
    parser.add_argument('-i', '--input', 
                        required = True, 
                        help = "Input file path for the stage")
    
    parser.add_argument('-o', '--output', 
                        required = True, 
                        help = "Output file path for the stage")
    
    parser.add_argument('--config', 
                        default = '/home/tanu/git/protvar/scripts/pv_db.ini', 
                        help = "Path to the configuration file")

    parser.add_argument('-v', 
                        '--verbose', 
                        action = 'store_true',
                        help = "Enable verbose output")
    
    # New argument to control whether to write the tracking data
    parser.add_argument('--write_counts', 
                        action = 'store_true',
                        help ="Enable or disable output of tracking counts for merged data in Stage 3")

    return parser

def load_configuration(config_path):
    """
    Load database configuration from a file.
    """
    config = configparser.ConfigParser()
    config.read(config_path)
    return dict(config['db'])

def run_stage1(input_file, output_file, verbose):
    print("Starting Stage 1: Uniprot and RefSeq ID extraction")
    start_time = time.time()
    
    data_list = read_file_to_list_of_dicts(input_file)
    updated_data_list = extract_ids(data_list)
    #write_list_of_dicts_to_tsv(updated_data_list, output_file)
    write_list_of_dicts(updated_data_list, output_file, delimiter = "\t")
    
    if verbose:
        print(f"Input file: {input_file}")
        print(f"Output file: {output_file}")

        elapsed = time.time() - start_time
        print(f"Stage 1 completed in {elapsed:.2f} seconds. \nData is processed and written to {output_file}")

def run_stage2(input_file, output_file, config_path, verbose):
    print("Starting Stage 2: Fetching missing UniProt IDs")

    # Load the configuration from the specified path
    db_params = load_configuration(config_path)

    start_time = time.time()
    final_data, tracking_dict = process_biogrid_file(input_file, db_params)
    #write_list_of_dicts_to_tsv(final_data, output_file)  
    
    # Write data with missing uniprot IDs retrieved from DB 
    write_list_of_dicts(final_data, output_file, delimiter = "\t" )  

    #OPTIONAL: you can write out the tracking dict
    dirname, fname = os.path.split(output_file)
    tracking_filename = os.path.join(dirname, os.path.splitext(fname)[0] + '_tracking_info_db.tsv')
    print("Writing tracking dict t: {tracking_filename}")
          
    write_nested_dict(data = tracking_dict, 
                       output_file = tracking_filename,
                       key_column_name = "Index", 
                       export_keys = True,
                       delimiter = '\t', 
                       merged_value_delimiter = ',')

    if verbose:
        print(f"Input file: {input_file}")
        print(f"Output file: {output_file}")
        elapsed = time.time() - start_time
        print(f"Stage 2 completed in {elapsed:.2f} seconds. Data is updated and written to {output_file}")

def run_stage3(input_file, output_file, verbose, write_counts):
    start_time = time.time()
    print("Starting Stage 3: Merging data values of specified columns")

    cols_for_value_merge = [
        'Publication_Identifiers',
        'Taxid_Interactor_A',
        'Taxid_Interactor_B',
        'Interaction_Types',
        'Source_Database',
        'Interaction_Identifiers',
        'Confidence_Values'
    ]

    # Merging data based on specified columns and tracking counts for cross check
    value_merged_bg_data, value_merged_counts = dict_value_merge(input_file = input_file,
                                                                 key_column = 'interaction_id',
                                                                 selected_columns = cols_for_value_merge)

    # write_nested_dict(data = value_merged_bg_data,
    #                   output_file = output_file,
    #                   delimiter = '\t',
    #                   merged_value_delimiter = ',')
    
    write_nested_dict(data = value_merged_bg_data, 
                      output_file = output_file, 
                      key_column_name = "interaction_id", 
                      export_keys = True,
                      delimiter = '\t', 
                      merged_value_delimiter = ',')

    
    print(f"Data with merged values written to: {output_file}")

    # OPTIONAL: write tracking information if the flag is set  
    if write_counts:
        print("Writing tracking dict with counts as well...")
        dirname, fname = os.path.split(output_file)
        counts_filename = os.path.join(dirname, os.path.splitext(fname)[0] + '_COUNTS.tsv')
        
        print(f"Tracking counts written to: {counts_filename}")
        write_nested_dict(data = value_merged_counts, 
                          output_file = counts_filename, 
                          key_column_name = "interaction_id", 
                          export_keys = True,
                          delimiter = '\t', 
                          merged_value_delimiter = ',')

    if verbose:
        print(f"Input file: {input_file}")
        print(f"Output file: {output_file}")
        print(f"Using the following columns to merge values: {cols_for_value_merge}")
        elapsed = time.time() - start_time
        print(f"Stage 3 completed in {elapsed:.2f} seconds. \nMerged data is written to {output_file}")


def main():
    parser = setup_arguments()
    args = parser.parse_args()
  
    # Execute the appropriate stage
    if args.stage == 1:
        run_stage1(args.input, args.output, args.verbose)
    elif args.stage == 2:
        run_stage2(args.input, args.output, args.config, args.verbose)
    elif args.stage == 3:
        run_stage3(args.input, args.output, args.verbose, args.write_counts)
if __name__ == "__main__":
    main()
###############################################################################