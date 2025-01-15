#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 12:33:41 2025

@author: tanu
"""

import os
import sys
import csv
from matplotlib_venn import venn2, venn3
import matplotlib.pyplot as plt
###############################################################################
############
# Function: 
############
def read_keys_from_file(file_path, column_name):
    """
    Reads a file and extracts a set of unique values from a specified column.
    Used for Venn Diagram plotting primarily.

    Args:
        file_path (str): Path to the file to read.
        column_name (str): Name of the column from which to extract unique values.

    Returns:
        set: A set of unique values from the specified column.
    """
    unique_values = set()
    try:
        with open(file_path, mode = 'r', newline = '') as file:
            reader = csv.DictReader(file, delimiter='\t')  # Adjust delimiter if necessary
            for row in reader:
                if column_name in row and row[column_name] != 'N/A':
                    unique_values.add(row[column_name])
    except FileNotFoundError:
        sys.exit(f"Error: File not found {file_path}")
    except KeyError:
        sys.exit(f"Error: Column {column_name} not found in file {file_path}")
    return unique_values

############
# Function: 
############
def plot_venn(*sets, labels = None, plot_title = 'Venn Diagram of Sets', output_file = None):
    """
    Plots a Venn diagram for two or three sets.
    
    Args:
        *sets: Variable number of set arguments.
        labels (list of str, optional): Labels for the sets in the Venn diagram. Must match the number of sets if provided.
        plot_title (str): Title of the plot.
        output_file (str, optional): Path to save the output file. If None, the plot is not saved.
    """
    
    if len(sets) not in [2, 3]:
        sys.exit("Error: Venn diagrams can only be plotted for two or three sets.")
    
    if labels and len(labels) != len(sets):
        sys.exit(f"Error: Number of labels {len(labels)} does not match number of sets {len(sets)}.")
    
    plt.figure(figsize=(8, 8))
    if len(sets) == 2:
        venn2(sets, labels)
    elif len(sets) == 3:
        venn3(sets, labels)
    plt.title(plot_title)
    
    if output_file:
        plt.savefig(output_file)  # Save the plot to the specified file
        print(f"Plot saved to {output_file}")
    
    plt.show()  # Display the plot

# Example usage
#set1 = {'a', 'b', 'c', 'd'}
#set2 = {'b', 'c', 'e'}
#set3 = {'c', 'd', 'e', 'f'}

#plot_venn(set1, set2, set3, 
#          labels = ['Group 1', 'Group 2', 'Group 3'], 
#          plot_title = 'Sample Venn Diagram')
###############################################################################