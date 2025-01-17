#!/usr/bin/env python3

# -*- coding: utf-8 -*-
"""
Created on Mon Dec 30 14:25:10 2024

@author: tanu
"""
###############################################################################
# Merges or concatenates two data files based cmd option. 
# It provides the flexibility to either merge the data based on a common key column (performing an inner join) 
# or concatenate them vertically (similar to an `rbind` in R).

# Example Usage: Merging Data with a Common Column 
#    pv_data_merging.py --input_file_list file_list.txt --common_col interaction_id --outfile merged_data.tsv --verbose

# Example Usage: Concatenating Data Without a Common Column
#    pv_data_merging.py --input_file_list file_list.txt --outfile concatenated_output.tsv --verbose

# Arguments:
    
# - `--input_file_list`: Specifies a file containing paths to two input files to be merged, one per line. Max of 2 files
# - `--common_col`: Specifies the column to be used as a key for merging.
#                   If not provided, the script concatenates the files.
# - `--outfile`: Specifies the path for the merged output file.
# - `--verbose`: Enables verbose output to provide detailed processing information including time.

# Features:

# - Merge Based on Common Column: If a common column is specified, the script checks for its presence in both files and merges the data based on this column.
# - Concatenate Data: If no common column is specified, the script concatenates the data from both files row-wise.
# - Verbose Output: Optional verbose output to track the script’s processing details and timings.
# - Flexible Output: Allows specification of the output file path.

# Custom functions imported are clearly indicated in the imports

###############################################################################
import argparse
import os
import sys
import csv
import time
#from time import time

# my custom functions
#from my_pv import *
from my_pv.pv_functions import read_file_to_dict, merge_dicts, write_nested_dict, subset_dict_of_dicts
###############################################################################
############
# Function: 
############
def read_paths_from_file(file_path):
    """
    Reads a file where each line is a file path and returns a list of these paths.
    """
    with open(file_path, 'r') as file:
        return [line.strip() for line in file if line.strip()]
############
# Function: 
############
def check_column_in_file(file_path, column_name, delimiter = '\t'):
    """
    Checks if a given column name exists in the file header.
    """
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter = delimiter)  # Adjust delimiter if necessary
        headers = next(reader)
        return column_name in headers
###############################################################################
def main():
    parser = argparse.ArgumentParser(description="Merge two data files with optional control over merging behavior.")
    
    parser.add_argument('--input_file_list', 
                        help = 'File containing list of input files to merge, one per line. Max of 2 files')
    parser.add_argument('--common_col', 
                        default = None, # 'interaction_id
                        help ='Optional common column to use for merging. If not provided, data files are concatenated.')
    parser.add_argument('--outfile', 
                        required = True, 
                        help = 'File path for the merged output.')
    parser.add_argument('--verbose', 
                        action = 'store_true', 
                        help = 'Enable verbose output.')

    args = parser.parse_args()
    
    if args.verbose:
        start_time = time.time()
        print(f"BEGUN time: {start_time:.2f} seconds")

    # Read file paths from the input file list
    file_paths = read_paths_from_file(args.input_file_list)
    print(f"Files to be merged/combined: {file_paths}")
    
    if len(file_paths) != 2:
        raise ValueError("Exactly two file paths are required.")

    ###############
    # MERGING: common_col provided as cmd argument
    ###############             
    if args.common_col:
        # Verify the presence of the common column in each file
        for file_path in file_paths:
            if not check_column_in_file(file_path, args.common_col):
                raise ValueError(f"The specified common column '{args.common_col}' is not found in the file {file_path}")
            else:
                print(f"Common column '{args.common_col}' found in {file_path}")
        
        # Read files into dictionaries using the common column as a unique key
        dict_list = [read_file_to_dict(file, unique_key = args.common_col) for file in file_paths]
        print(f"\nMerging files based on the common column provided: {args.common_col} which is merging on dict keys")
        
        if args.verbose:
            print(f"Files imported as dicts...\nNo. of keys in Dict 1: {len(dict_list[0])}")
            #print(f"\nFirst few entries in Dict 1:")
            #s1 = dict(list(dict_list[0].items())[0: 2]) 
            #pp.pprint(s1)
            print(f"\nNo. of keys in Dict 2: {len(dict_list[1])}")
            #print(f"\nFirst few entries in Dict 2:")
            #s2 = dict(list(dict_list[1].items())[0: 2]) 
            #pp.pprint(s2)
        
        #================
        # Inner join to keep all cols from primary dict
        #================
        primary_dict = dict_list[0]
        for supplemental_dict in dict_list[1:]:
            primary_dict, _ = merge_dicts(primary_dict, 
                                          supplemental_dict, 
                                          join_type = 'inner')
            
         # write_dict_to_tsv(data = primary_dict, 
         #                   filename = args.outfile, 
         #                   include_keys = True)
        
        write_nested_dict(data = primary_dict, 
                           output_file = args.outfile, 
                           key_column_name = "interaction_id", 
                           export_keys = False,
                           delimiter = '\t', 
                           merged_value_delimiter = ',')
        
        if args.verbose:
            print(f"combined dict length: {len(primary_dict)}")
            elapsed = time.time() - start_time
            print(f"END time: {time.time()}")
            print(f"ELAPSED TIME: {elapsed:.2f} seconds.\nData written to {args.outfile}")
        print(f"Merging complete")

    ###############
    # CONCATENATE: common_col NOT provided as cmd argument
    ###############           
    else:
        print(f"\nArgument common_col is empty. \nSo, concatenating files (like rbind)")
        #Read files into dictionaries using interaction_id as unique key so the dict structure serves with interaction_id as key
        dict_list = [read_file_to_dict(file, unique_key = 'interaction_id') for file in file_paths]
        if args.verbose:
            print(f"Files imported as dicts...\nNo. of keys in Dict 1: {len(dict_list[0])}")
            #print(f"\nFirst few entries in Dict 1:")
            #s1 = dict(list(dict_list[0].items())[0: 2]) 
            #pp.pprint(s1)
            print(f"\nNo. of keys in Dict 2: {len(dict_list[1])}")
            #print(f"\nFirst few entries in Dict 2:")
            #s2 = dict(list(dict_list[1].items())[0: 2]) 
            #pp.pprint(s2)

        try:
            #================
            # Dict comprehension
            #================
            result_dict = {**dict_list[0], **dict_list[1]} if len(dict_list) == 2 else {}
            
            # Filter result_dict to keep only specific keys as these are common between two files
            # otherwise directly writing the dict to file will result in lots of empty column values
            keys_to_keep = ['interaction_id', 'pdb', 'iptm', 'pdockq', 'pdockq_fd', 'sources']
            
            final_dict = subset_dict_of_dicts(result_dict, keys_to_keep)
            if args.verbose:
                print(f"No. of output columns: {len(keys_to_keep)}, \nOutput columns are: {keys_to_keep}")
                print(f"Concatenated dict length: {len(final_dict)}")

               
             # write_dict_to_tsv(data = final_dict, 
             #                   filename = args.outfile, 
             #                   include_keys = True,
             #                   ordered_columns = None)
            
            write_nested_dict(data = final_dict, 
                               output_file = args.outfile, 
                               key_column_name = "interaction_id", 
                               export_keys = False,
                               delimiter = '\t', 
                               merged_value_delimiter = ',')
            
            print(f"Concatenating complete")

        except ValueError as e:
            print(f"Error: {e}")

    if args.verbose:
        elapsed = time.time() - start_time
        print(f"END time: {time.time()}")
        print(f"ELAPSED TIME: {elapsed:.2f} seconds.\nData written to {args.outfile}")

if __name__ == '__main__':
    main()
