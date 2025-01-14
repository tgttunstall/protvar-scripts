#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 14 16:25:14 2025

@author: tanu
"""
import argparse
import os
import sys
from time import time
#from matplotlib_venn import venn2, venn3
#import matplotlib.pyplot as plt
from my_pv import * # my custom functions (already imports matplotlib)

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
    parser.add_argument('-k', '--key_column', 
                        required = True, 
                        help = 'Column name to extract keys for Venn diagram.')
    parser.add_argument('-l', '--labels', 
                        nargs = '+', 
                        help = 'Labels for each set in the Venn diagram.')
    parser.add_argument('-t', '--plot_title', 
                        default = 'Venn Diagram of Sets', 
                        help = 'Title for the Venn diagram.')
    parser.add_argument('--output_file', 
                        default = 'output_venn.png', 
                        help = "Output file for the Venn diagram.")
    return parser


def main():
    parser = setup_arguments()
    args = parser.parse_args()
    
    # Ensure common key column is present in all files
    if not common_columns(args.files, delimiters = 'auto', force = False):
        sys.stderr.write(f"Error: Specified key column '{args.key_column}' is not common across all files.\n")
        sys.exit(1)

    # Extract keys from each file
    sets = [read_keys_from_file(file, args.key_column) for file in args.files]
    
    if args.labels and len(args.labels) != len(sets):
        print("Error: The number of labels provided does not match the number of files.")
        return
    
    # Create labels if not provided
   # if not args.labels:
    #    args.labels = [f"Set {i+1}" for i in range(len(sets))]
    
    # Generate and save Venn diagram
    plot_venn(*sets, labels = args.labels, plot_title = args.plot_title, output_file = args.output_file)
    print(f"Venn diagram saved to {args.output_file}")

if __name__ == '__main__':
    main()