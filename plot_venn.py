#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 16:25:14 2025

@author: tanu
"""
import argparse
import os
import sys
import time
#from time import time # does not work!
#from my_pv import common_columns
#from venn_functions import *
from my_pv.pv_functions import common_columns
from my_pv.venn_functions import read_keys_from_file, plot_venn

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