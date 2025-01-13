#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 14:10:38 2024

@author: tanu
"""
###############################################################################
import os
import sys
import argparse
import configparser
import time
from my_pv import * # my custom functions

###############################################################################
def setup_arguments():
    parser = argparse.ArgumentParser(description="Process BioGRID data in multiple stages.")
    
    parser.add_argument('-s', '--stage', type = int, choices = [1, 2, 3], required = True, 
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
                        help = "Enable output of tracking counts for merged data in Stage 3")

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
    write_list_of_dicts_to_tsv(updated_data_list, output_file)
    
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
    write_list_of_dicts_to_tsv(final_data, output_file)
    # OPTIONAL: you can also write out the tracking dict for debug

    if verbose:
        print(f"Input file: {input_file}")
        print(f"Output file: {output_file}")
        elapsed = time.time() - start_time
        print(f"Stage 2 completed in {elapsed:.2f} seconds. Data is updated and written to {output_file}")

def run_stage3(input_file, output_file, verbose, write_counts = True):
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

    # Optionally write tracking information if the flag is set  
    if write_counts:
        print("Writing tracking dict with counts as well...")
        dirname, fname = os.path.split(output_file)
        counts_filename = os.path.join(dirname, os.path.splitext(fname)[0] + '_COUNTS.tsv')
        write_nested_dict(data = value_merged_counts,
                          output_file = counts_filename,
                          delimiter = '\t',
                          merged_value_delimiter = ',')
        if verbose:
            print(f"Tracking counts written to {counts_filename}")

    if verbose:
        print(f"Input file: {input_file}")
        print(f"Output file: {output_file}")
        print(f"Using the following columns to merge values: {cols_for_value_merge}")
        elapsed = time.time() - start_time
        print(f"Stage 3 completed in {elapsed:.2f} seconds. Merged data is written to {output_file}")


def main():
    parser = setup_arguments()
    args = parser.parse_args()

    # Configuration and setup
    #db_params = load_configuration(args.config)
    
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