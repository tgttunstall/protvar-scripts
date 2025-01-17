#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 16:25:14 2025

@author: tanu
"""
###############################################################################
"""
Generates Venn Diagrams for two or three data sets based on a specified key column which should be common across the datasets.
It accepts multiple file paths, extracts unique identifiers from each file based on a specified key column, 
and plots a Venn diagram to visually represent the overlap among these sets.
Optionally users can provide labels for each dataset and a title for the Venn diagram.

Usage:
    python plot_venn.py -f FILE1.tsv FILE2.tsv [--labels Label1 Label2] --key_column ColumnName [-t "Custom Title"] [--output_file "output.png"] [--verbose]

Arguments:
    -f, --files: List of paths to input files.
    --key_column: The column name to extract keys for Venn diagram
    --labels: Optional labels for each set in the Venn diagram. If omitted, default labels ("Set 1", "Set 2", etc.) will be generated.
    -t, --plot_title: Optional title for the Venn diagram. Default is 'Venn Diagram of Sets'.
    --output_file: Path to save the generated Venn diagram image. Default is 'output_venn.png'.
    --verbose: Optional Prints detailed logs about the script's execution, including timings.

Features:
    - Handles two or three input files to generate corresponding Venn diagrams.
    - Provides error handling for missing key columns across files.
    - Automatically handles delimiter detection for CSV and TSV file formats.
    - Supports optional verbose output for detailed operation insight.

This script utilizes the 'matplotlib_venn' library to plot Venn diagrams and requires CSV files to be formatted correctly with headers.

Example Command:
    python plot_venn.py -f data1.csv data2.csv data3.csv --key_column "ID" --labels "Data Set 1" "Data Set 2" "Data Set 3" -t "Comparison of IDs" --verbose
"""
###############################################################################
import argparse
import os
import sys
import time
#from time import time # does not work!
#from my_pv import common_columns
#from venn_functions import *
from my_pv.pv_functions import common_columns
from my_pv.venn_functions import read_keys_from_file, plot_venn
###############################################################################
def setup_arguments():
    """
    Sets up command line arguments for the script.
    
    Returns:
        ArgumentParser: The configured argument parser.
    """
    parser = argparse.ArgumentParser(description = "Generate Venn diagrams from provided CSV/TSV files.")
    parser.add_argument('-f', '--files', 
                        nargs = '+', 
                        required = True, 
                        help = 'Paths to the input files.')
    parser.add_argument('--key_column', 
                        required = True, 
                        help = 'Column name to extract keys for Venn diagram.')
    parser.add_argument('--labels', 
                        nargs = '*', 
                        help = 'Optional labels for the sets in the Venn diagram. If not provided, labels will be generated automatically.')
    parser.add_argument('-t', '--plot_title', 
                        default = 'Venn Diagram of Sets', 
                        help = 'Title for the Venn diagram.')
    parser.add_argument('--output_file', 
                        default = 'output_venn.png', 
                        help = "Output file for the Venn diagram.")
    parser.add_argument('--verbose', 
                        action = 'store_true', 
                        help = 'Enable verbose output.')
    return parser


def main():
    
    parser = setup_arguments()
    args = parser.parse_args()
    
    start_time = time.time()
    
    if args.verbose:
        print(f"\nReading files: '{args.files}'")
        
    # Ensure common key column is present in all files
    print(f"\nChecking if '{args.key_column}' is present in files provided...")
    
    common_cols = common_columns(args.files, delimiters = 'auto', force = False)
    
    if common_cols is None or args.key_column not in common_cols:
        sys.stderr.write(f"Error: Specified key column '{args.key_column}' is not present/common across all files.\n")
        sys.exit(1)
    else:
        print(f"'{args.key_column}' is present in files. Proceeding to generate plots...")

    # Extract keys from each file
    if args.verbose:
        print(f"\nStarting the process to read sets for plots.")
        
    sets = [read_keys_from_file(file, args.key_column) for file in args.files]
    
    if args.verbose:
        print(f"\nSets read successfully. Elapsed time: '{time.time() - start_time:.2f}' seconds.")
    
    # Labels
    if args.labels and len(args.labels) != len(sets):
        print("\nError: The number of labels provided does not match the number of files.")
        return
    elif not args.labels:
        args.labels = [f"Set {i+1}" for i in range(len(sets))]
        
    if args.verbose:
        print(f"\nGenerating Venn diagram with labels: {args.labels}")
        
    # Generate and save Venn diagram
    plot_venn(*sets, labels = args.labels, plot_title = args.plot_title, output_file = args.output_file)
    
    if args.output_file:
        print(f"\nVenn diagram saved to {args.output_file}")
    
    if args.verbose:
        print(f"\nVenn diagram successfully generated. Total elapsed time: {time() - start_time:.2f} seconds.")

if __name__ == '__main__':
    main()