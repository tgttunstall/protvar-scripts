#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 14:07:18 2024

@author: tanu
"""
import sys
import time
import csv
import re
import pprint as pp
from collections import defaultdict, Counter
import psycopg2
#import matplotlib.pyplot as plt
#from matplotlib_venn import venn2, venn3
###############################################################################
############
# Function: 
############
def eprint(*args, **kwargs):
    """
    Prints provided arguments to standard error (stderr).

    Args:
        *args: Variable length argument list.
        **kwargs: Arbitrary keyword arguments.
    """
    print(*args, file = sys.stderr, **kwargs)
###############################################################################
############
# Function: 
############
def secs2time(secs):
    """
    Converts a time duration in seconds to a human-readable format.

    Args:
        secs (float): Time duration in seconds.

    Returns:
        str: A formatted string representing the time in "HHh MMm SSs" format.
    """
    minutes, seconds = divmod(secs, 60)
    hours, minutes = divmod(minutes, 60)
    return "{:02.0f}h {:02.0f}m {:02.0f}s".format(hours, minutes, seconds)
###############################################################################
#############
# Function: subset n entries from dict
#############
def subset_dict(input_dict, num_entries):
    """
    This function subsets a dictionary based on x number of entries.
    
    Parameters:
    - input_dict: The original dictionary to be subsetted.
    - num_entries: number of entries to be subsetted 
    
    Returns:
    - A new dictionary containing only the specified number of entries.
    """
    # Initialize an empty dictionary to store the subset
    subset = {}
 
    # Get the specified number of keys
    keys = list(input_dict.keys())[:num_entries]
    
    # Iterate over the keys and add them to the subset dictionary
    for key in keys:
        subset[key] = input_dict[key]
    
    return subset
###############################################################################
#############
# Function: subset dict based on specific keys
############
def subset_dict_by_keys(input_dict, keys):
    """
    This function subsets a dictionary based on specific keys.
    
    Parameters:
    - input_dict: The original dictionary to be subsetted.
    - keys: A list of keys to include in the subset.
    
    Returns:
    - A new dictionary containing only the specified keys.
    """
    # Initialize an empty dictionary to store the subset
    subset = {}

    # Iterate over the specified keys and add them to the subset dictionary if they exist in the input_dict
    for key in keys:
        if key in input_dict:
            subset[key] = input_dict[key]

    return subset

# Example usage
#a = subset_dict_by_keys(ex2D, [0,1])
###############################################################################
###########
# Function: 
###########   
def common_columns(files, delimiters = 'auto', force = False):
    """
    Determines the common columns between the headers of multiple files, using either automatically determined delimiters
    based on file extensions or user-specified delimiters. Can force processing despite warnings.

    Args:
        files (list of str): List of file paths.
        delimiters (list of str or 'auto'): A list of delimiters corresponding to each file or 'auto' to determine automatically.
        force (bool): If True, continues processing even if potential issues are detected with the delimiters.
            potential issue is if it returns a single column which likely indicates the delimiter is incorrect.
        

    Returns:
        set: Set of common column names, or None if an unrecoverable error occurs.
    """
    if delimiters == 'auto':
        delimiters = []
        for file in files:
            if file.endswith('.csv'):
                delimiters.append(',')
            elif file.endswith('.tsv') or file.endswith('.txt'):
                delimiters.append('\t')
            else:
                print(f"Error: File type of {file} not recognized. Processing terminated.")
                return None
    elif len(delimiters) != len(files):
        print("Error: The number of delimiters must match the number of files.")
        return None

    headers = []
    for file, delimiter in zip(files, delimiters):
        try:
            with open(file, 'r', newline='') as f:
                reader = csv.reader(f, delimiter = delimiter)
                header = next(reader)
                if len(header) == 1 and not force:
                    print(f"Warning: Check delimiter for {file}. Only one column detected, possibly incorrect delimiter. \nRun with 'force=True' if you wish to continue...")
                    if not force:
                        return None
                headers.append(set(header))
        except Exception as e:
            print(f"Error reading file {file}: {e}")
            return None

    # Compute the intersection of all header sets
    if headers:
        common_headers = set.intersection(*headers)
        return common_headers
    else:
        print("\nNo common columns found...")
        return set()

#Example usage:
#file_paths = ['/home/pub/Work/data_arise_proteome/protvar/file1.csv',
#              '/home/pub/Work/data_arise_proteome/protvar/file0.csv', 
#              '/home/pub/Work/data_arise_proteome/protvar/file2.csv', 
#              '/home/pub/Work/data_arise_proteome/protvar/file2_c.csv', 
#              '/home/pub/Work/data_arise_proteome/protvar/file3.tsv',
#              '/home/pub/Work/data_arise_proteome/protvar/file4.tsv']

#my_delimiters_wrong = ['\t', ',', '\t']
# FIXME, what if the delimiter is not correct?
#my_delimiters = [',',',',',','\t', '\t'] 

#(common_columns(file_paths, delimiters=my_delimiters_wrong, force = True))  # Uses the list of delimiters
#print(common_columns(file_paths, delimiters = my_delimiters, force = True))  # Uses the list of delimiters
#print(common_columns(file_paths, delimiters = 'auto'))  # Uses a single delimiter for all
#print(common_columns(file_paths, delimiters = 'auto', force = True))  # Uses a single delimiter for all

#with open('/home/pub/Work/data_arise_proteome/protvar/file3.tsv', newline = '') as f:
    #header = next(csv.reader(f))
###############################################################################
############
# Function: 
############
def read_file_to_dict(file_path, unique_key = None, delimiter = '\t'):
    """
    Reads a delimited file and converts it into a dictionary of dictionaries, ensuring all columns are represented.
    Missing data for any column is filled with 'N/A'. If no unique_key is provided, each row's index is used as the key.

    Args:
        file_path (str): Path to the file.
        unique_key (str, optional): Column name to use as the unique key for the dictionary.
        delimiter (str, optional): Delimiter character for splitting columns in the file.

    Returns:
        dict: A dictionary of dictionaries, where each key is the value from the unique_key column or the row index,
              and each dictionary entry contains all columns with missing values filled as 'N/A'.
    """
    try:
        with open(file_path, mode = 'r', newline = '') as file:
            #print(f"Reading...: {file_path}")
            reader = csv.DictReader(file, delimiter = delimiter)
            headers = reader.fieldnames
            default_dict = {header: 'N/A' for header in headers}
            result = {}
            for index, row in enumerate(reader):
                complete_row = default_dict.copy()
                complete_row.update(row)
                key = complete_row[unique_key] if unique_key else index
                if key in result:
                    sys.stderr.write(f"Duplicate key found: {key}\n")
                result[key] = complete_row
            return result
    except Exception as e:
        sys.stderr.write(f"Failed to read {file_path}: {e}\n")
        sys.exit(1)
        
###############################################################################
############
# Function: 
############
def merge_dicts(primaryD, supplementalD, join_type = 'inner'):
    """
    Merges entries from two dictionaries based on the specified SQL-style join type.

    Args:
        primaryD (dict): Primary dictionary where data will be merged into.
        supplementalD (dict): Supplemental dictionary from which data will be merged.
        join_type (str): Type of join to perform ('inner', 'left', 'right', 'full').

    Returns:
        dict: Dictionary containing merged data based on the specified join type.

    Description of Join Types:
        - 'inner': Merges only the keys that exist in both primaryD and supplementalD. 
                   The result will only include keys found in both dictionaries. If a key is found 
                   in both, values from supplementalD will overwrite values from primaryD for overlapping keys.

        - 'left': Merges all keys from the primaryD and only the keys from supplementalD that exist in primaryD.
                  The result includes all keys from the primaryD. For keys that exist in both dictionaries, 
                  values from supplementalD will overwrite values from primaryD for those keys.

        - 'right': Merges all keys from the supplementalD and only the keys from primaryD that exist in supplementalD.
                   The result includes all keys from the supplementalD. For keys that exist in both dictionaries, 
                   values from primaryD will overwrite values from supplementalD for those keys.

        - 'full': Merges all keys from both primaryD and supplementalD.
                  The result includes all keys from both dictionaries. For keys that exist in both, 
                  values from supplementalD will overwrite values from primaryD. For keys unique to each dictionary, 
                  they are included as-is from their respective source.

    Each join type handles the merging of keys and their associated values differently.
    """
    result_dict = {}
    overwritten_values = {}  # Dictionary to store overwritten values

    # Define the keys to process based on the join type
    keys = {
        'inner': set(primaryD.keys()) & set(supplementalD.keys()),
        'left': set(primaryD.keys()),
        'right': set(supplementalD.keys()),
        'full': set(primaryD.keys()) | set(supplementalD.keys())
    }.get(join_type, set())

    # Merge data based on the join type
    for key in keys:
        entry_from_primary = primaryD.get(key, {})
        entry_from_supplemental = supplementalD.get(key, {})

        # Prepare to merge by copying data from primary
        merged_entry = entry_from_primary.copy()

        # Check for and track overwritten values
        for k, v in entry_from_supplemental.items():
            if k in merged_entry:
                if merged_entry[k] != v:
                    if key not in overwritten_values:
                        overwritten_values[key] = {}
                    overwritten_values[key][k] = merged_entry[k]
            merged_entry[k] = v

        # Assign the merged entry to the result dictionary
        result_dict[key] = merged_entry

    return result_dict, overwritten_values

# Example usage:
#primaryD = {'a': {'x': 1}, 'b': {'x': 2}}
#supplementalD = {'b': {'y': 3}, 'c': {'y': 4}}

#print("Inner Join:", merge_dicts(primaryD, supplementalD, 'inner'))
#print("Left Join:", merge_dicts(primaryD, supplementalD, 'left'))
#print("Right Join:", merge_dicts(primaryD, supplementalD, 'right'))
#print("Full Join:", merge_dicts(primaryD, supplementalD, 'full'))

# Example usage
#primaryD = {'a': {'x': 1, 'z': 3}, 'b': {'x': 2, 'z': 5}}
#supplementalD = {'b': {'x': 3, 'z': 4}, 'c': {'y': 4}}
#pp.pprint(merge_dicts(primaryD, supplementalD, 'inner'))
#pp.pprint(merge_dicts(primaryD, supplementalD, 'left'))
#pp.pprint(merge_dicts(primaryD, supplementalD, 'right'))
#pp.pprint(merge_dicts(primaryD, supplementalD, 'full'))

#result, overwrites = merge_dicts(primaryD, supplementalD, 'inner')
#print("Merged Data:", result)
#print("Overwritten Values:", overwrites)
###############################################################################
############
# Function: 
############
def subset_dict_of_dicts(data, keys_to_keep):
    """
    Subsets a dictionary of dictionaries, keeping only specific keys in the inner dictionaries.

    Args:
        data (dict): The original dictionary of dictionaries.
        keys_to_keep (list or set): The keys to keep in each inner dictionary.

    Returns:
        dict: A new dictionary with each inner dictionary containing only the specified keys.
        
      """
    # Create a new dictionary with only the desired keys in the inner dictionaries
    return {outer_key: {inner_key: inner_val for inner_key, inner_val in inner_dict.items() if inner_key in keys_to_keep}
            for outer_key, inner_dict in data.items()}

# # Example usage
# data = {
#     'item1': {'key1': 'value1', 'key2': 'value2', 'key3': 'value3'},
#     'item2': {'key1': 'value4', 'key2': 'value5', 'key3': 'value6'},
#     'item3': {'key1': 'value7', 'key2': 'value8', 'key4': 'value9'},
#     'item4': {'key1': 'value7', 'key4': 'value11', }

# }

# keys_to_keep = ['key1', 'key3']
# subset_data = subset_dict_of_dicts(data, keys_to_keep)
# import pprint as pp
# pp.pprint(subset_data)

###############################################################################
############
# Function: 
############
def read_file_to_list_of_dicts(filename):
    """
    Reads a TSV (Tab-Separated Values) file and converts each row into a dictionary,
    with the column headers serving as keys. This function returns a list of dictionaries,
    each representing a single row of data.
    
    Args:
        filename (str): The path to the TSV file that will be read.
    
    Returns:
        list: A list of dictionaries where each dictionary corresponds to a row in the file.
              The keys in each dictionary are derived from the column headers of the TSV file.
    
    Example:
        File: 'data.tsv' with the following content:
            Name    Age    City
            Alice   30     New York
            Bob     25     Los Angeles
            Charlie 35     Chicago
    
        Calling this function like so:
            data_list = read_file_to_list_of_dicts('data.tsv')
    
        Would return:
            [
                {'Name': 'Alice', 'Age': '30', 'City': 'New York'},
                {'Name': 'Bob', 'Age': '25', 'City': 'Los Angeles'},
                {'Name': 'Charlie', 'Age': '35', 'City': 'Chicago'}
            ]
    """
    with open(filename, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')  # Assuming tab-delimited
        data_list = [row for row in reader]  # List comprehension to collect all data rows
    return data_list

# Example Usage
#filename = '/home/pub/Work/data_arise_proteome/protvar/biogrid/sample_biogrid.txt'
#data_list = read_file_to_list_of_dicts(filename)
#pp.pprint(data_list[0])  # Print the first record 
###############################################################################
############
# Function: 
############
#def extract_and_combine_ids(data_list):
def extract_ids(data_list):
    """
    Extracts UniProt and RefSeq IDs from a list of dictionaries that contain biological data,
    specifically from fields 'Alt_IDs_Interactor_A' and 'Alt_IDs_Interactor_B'. Each entry
    in the list is updated with new keys for the extracted IDs ('uniprot_id1', 'refseq_id1',
    'uniprot_id2', 'refseq_id2').
    
    Args:
        data_list (list of dict): A list of dictionaries where each dictionary represents a biological
                                  interaction record. Each record must have the fields 'Alt_IDs_Interactor_A'
                                  and 'Alt_IDs_Interactor_B' containing the data strings to be searched.
    
    Returns:
        list of dict: The updated list of dictionaries, with each dictionary now including four new keys:
                      'uniprot_id1', 'refseq_id1', 'uniprot_id2', and 'refseq_id2'. Each of these fields
                      will contain the extracted ID if found, or 'N/A' if not found.
    
    Example:
       
    
    This function is particularly useful in bioinformatics for preprocessing data to isolate identifier information
    for further analysis or database integration.
    """

    # Compile regex patterns for efficiency in a loop
    uniprot_id_pattern = re.compile(r'uniprot/swiss-prot:([^\|]+)')
    refseq_id_pattern = re.compile(r'refseq:([NXW][CGMRP]_\d+)')

    # Iterate through each dictionary in the list
    for entry in data_list:
        # Process Alt_IDs_Interactor_A
        alt_ids_a = entry.get('Alt_IDs_Interactor_A', '')
        uniprot_id_a = uniprot_id_pattern.search(alt_ids_a)
        refseq_id_a = refseq_id_pattern.search(alt_ids_a)

        # Assigning new keys to the dictionary
        entry['uniprot_id1'] = uniprot_id_a.group(1) if uniprot_id_a else 'N/A'
        entry['refseq_id1'] = refseq_id_a.group(1) if refseq_id_a else 'N/A'

        # Process Alt_IDs_Interactor_B
        alt_ids_b = entry.get('Alt_IDs_Interactor_B', '')
        uniprot_id_b = uniprot_id_pattern.search(alt_ids_b)
        refseq_id_b = refseq_id_pattern.search(alt_ids_b)

        # Assigning new keys to the dictionary
        entry['uniprot_id2'] = uniprot_id_b.group(1) if uniprot_id_b else 'N/A'
        entry['refseq_id2'] = refseq_id_b.group(1) if refseq_id_b else 'N/A'

        # Create an interaction ID based on uniprot IDs
        #entry['interaction_id'] = f"{entry['uniprot_id1']}_{entry['uniprot_id2']}"

    return data_list

# Example usage
# input_data = [
#      {'Alt_IDs_Interactor_A': 'uniprot/swiss-prot:P12345|refseq:NP_123456',
#       'Alt_IDs_Interactor_B': 'uniprot/swiss-prot:P23456|refseq:NP_234567'},
#      {'Alt_IDs_Interactor_A': 'uniprot/swiss-prot:P34567',
#       'Alt_IDs_Interactor_B': 'refseq:NP_345678'}
#  ]
# extracted_data = extract_ids(input_data)

# for entry in extracted_data:
#     print(entry['uniprot_id1'], entry['refseq_id1'], entry['uniprot_id2'], entry['refseq_id2'])
    
# Output:
#P12345 NP_123456 P23456 NP_234567
#P34567 N/A N/A NP_345678
###############################################################################
############
# Function: 
############
def write_list_of_dicts(data_list, output_filename, delimiter = "\t"):
    """
    Writes a list of dictionaries to a TSV file.

    Args:
        data_list (list of dict): List of dictionaries containing the data.
        output_filename (str): The filename where the TSV data will be saved.
    """
    if not data_list:
        print("No data to write.")
        return

    # Open the file for writing
    with open(output_filename, 'w', newline = '') as file:
        # Use the csv.DictWriter to handle the dictionary writing
        writer = csv.DictWriter(file, fieldnames = data_list[0].keys(), delimiter = delimiter)
        
        # Write the header automatically based on the keys of the first dictionary
        writer.writeheader()
        
        # Write the data
        writer.writerows(data_list)

# Example usage
#output_filename = 'biogrid_data.tsv'
#write_list_of_dicts(updated_data_list, output_filename, delimiter = "\t")

# Example usage
#input_file = '/home/pub/Work/data_arise_proteome/protvar/biogrid/updated_bg_human_interactions.tsv'
#input_file = '/home/pub/Work/data_arise_proteome/protvar/biogrid/sample_biogrid.txt'

# #result_data, merged_ids = process_biogrid_file(input_file, 
#                                    merge_fields = ['ID_Interactor_A',
#                                                  'ID_Interactor_B',
#                                                  'Taxid_Interactor_A',
#                                                  'Taxid_Interactor_B',
#                                                  'Interaction_Types',
#                                                  'Source_Database',
#                                                  'Interaction_Identifiers',
#                                                  'Confidence_Values'],
#                                    output_format = 'list')


#print("Merged Interaction IDs:", merged_ids)
###############################################################################
############
# Function: 
############
def dict_value_merge(input_file, key_column = 'interaction_id', selected_columns = None):
    """
    Merges values in a TSV file based on a specified key column and aggregates data into lists,
    while tracking the number of merges for each column under each key.

    Args:
        input_file (str): The path to the TSV file to process.
        key_column (str): The column name to use as the key for merging.
        selected_columns (list, optional): List of column names to merge. If None, all columns are merged.

    Returns:
        dict: A dictionary with keys as unique entries from the key column and values as dictionaries
              of data aggregated into lists.
        dict: A tracking dictionary that counts the number of entries merged for each column under each key.
    """
    merged_data = defaultdict(lambda: defaultdict(list))
    merge_counts = defaultdict(lambda: defaultdict(int))
    
    with open(input_file, 'r', newline = '') as file:
        reader = csv.DictReader(file, delimiter='\t')

        # Determine which columns to merge
        if selected_columns:
            all_columns = selected_columns
        else:
            all_columns = [col for col in reader.fieldnames if col != key_column]
                
        for row in reader:
            key_value = row[key_column]
            if key_value == 'N/A':
                continue  # Skip entries with invalid key values
    
            for col in all_columns:
                col_value = row[col]
                merged_data[key_value][col].append(col_value)
                merge_counts[key_value][col] += 1
                    
    # Convert defaultdict to dict for return
    return dict(merged_data), dict(merge_counts)

# Example usage:
#infile = '/home/pub/Work/data_arise_proteome/protvar/biogrid/updated_sample_biogrid_PVDB.tsv'
#key_column = 'interaction_id'
# my_cols = ['Publication_Identifiers',
# 'Taxid_Interactor_A',
# 'Taxid_Interactor_B',
# 'Interaction_Types',
# 'Source_Database',
# 'Interaction_Identifiers',
# 'Confidence_Values']

# data, counts = dict_value_merge(input_file = infile, key_column = 'interaction_id', selected_columns = my_cols)
# pp.pprint(data['P41143_Q5JY77'])
# a = subset_dict_by_keys(data, ['P41143_Q5JY77'])
# b = subset_dict_by_keys(counts, ['P41143_Q5JY77', 'A8MVS5_Q8ND76'])
    
# Asserts to check the correctness of data
#assert isinstance(data, dict)
#assert isinstance(counts, dict)
#assert 'P41143_Q5JY77' in data  # Check for expected keys
#assert len(data['P41143_Q5JY77']['Taxid_Interactor_A']) == counts['P41143_Q5JY77']['Taxid_Interactor_A']
#print("Test passed: Data and counts match expected values.")

###############################################################################
############
# Function: 
############
def write_nested_dict(data, 
                      output_file, 
                      key_column_name = "ID", 
                      export_keys = True, 
                      delimiter = "\t", 
                      merged_value_delimiter = ","):
    """
    Writes a nested dictionary to a TSV file.

    Args:
        data (dict): The nested dictionary where the outer key is an identifier (e.g., 'interaction_id'),
                     and each inner dictionary contains column names as keys with values as lists or single values.
        output_file (str): Path to the output TSV file.
        key_column_name (str): The name of the column that will store dictionary keys. Defaults to "ID".
        export_keys (bool): Whether to include the dictionary keys in the output file. Defaults to True.
        delimiter (str): Field delimiter for the TSV file. Defaults to '\t' for TSV format.
        merged_value_delimiter (str): Delimiter to join list values in the output. Defaults to ','.
    """
    if not data:
        raise ValueError("The input data dictionary is empty. No file will be written.")

    with open(output_file, 'w', newline='') as file:
        # Collect all possible fieldnames from the data to handle varying fields across entries
        all_fields = set()
        for columns in data.values():
            all_fields.update(columns.keys())

        # Ensure the key column is appropriately handled
        # If the dict of dicts has extra fields, then the output file will have all fieldnames with N/A where applicable
        fieldnames = [key_column_name] + sorted(all_fields) if export_keys else sorted(all_fields)
        writer = csv.DictWriter(file, fieldnames=fieldnames, delimiter=delimiter)
        writer.writeheader()

        for key, columns in data.items():
            row = {key_column_name: key} if export_keys else {}
            
            for field in fieldnames:
                if export_keys and field == key_column_name:
                    continue  # Skip as it's already handled
                values = columns.get(field, None)  # Use None to clearly indicate missing values
                if values is None:
                    row[field] = 'N/A'  # Explicitly set missing fields to 'N/A'
                elif isinstance(values, list):
                    row[field] = merged_value_delimiter.join(map(str, values))  # Handle list values
                else:
                    row[field] = str(values)  # Convert single values to string  
            writer.writerow(row)

# Example usage:
# test_dict = {'A8MVS5_Q8ND76': 
#               {'Confidence_Values': ['score:0.999690761'],
#               'Interaction_Identifiers': ['biogrid:3044379'],
#               'Interaction_Types': ['psi-mi:"MI:0915"(physical association)'],
#                                 'Publication_Identifiers': ['pubmed:33961781'],
#                                 'Source_Database': ['psi-mi:"MI:0463"(biogrid)'],
#                                 'Taxid_Interactor_A': ['taxid:9606'],
#                                 'Taxid_Interactor_B': ['taxid:9606']},
              
#               'P41143_Q5JY77': {'Confidence_Values': ['-','-'],
#                                 'Interaction_Identifiers': ['biogrid:11769',
#                                                             'biogrid:905123'],
#                                 'Interaction_Types': ['psi-mi:"MI:0407"(direct interaction)', 
#                                                       'psi-mi:"MI:0407"(direct interaction)'],
#                                 'Publication_Identifiers': ['pubmed:12142540',
#                                                             'pubmed:15086532'],
#                                 'Source_Database': ['psi-mi:"MI:0463"(biogrid)',
#                                                     'psi-mi:"MI:0463"(biogrid)'],
#                                 'Taxid_Interactor_A': ['taxid:9606',
#                                                       'taxid:9606'],
#                                 'Taxid_Interactor_B': ['taxid:9606',
#                                                       'taxid:9606']},
#               'A_B': 
#                {'Confidence_Values': 'score:0.999690761',
#                 'Interaction_Identifiers': 'biogrid:3044379',
#                 'Interaction_Types': 'psi-mi:"MI:0915"(physical association)',
#                 'Publication_Identifiers': 'pubmed:33961781',
#                 'Source_Database': 'psi-mi:"MI:0463"(biogrid)',
#                 'Extra_field': '123X'}
#               }
             
    
# write_nested_dict(data = test_dict, 
#                   output_file = "/home/pub/Work/data_ebi_protvar/test_dict_flattened.tsv", 
#                   key_column_name="interaction_id", 
#                   export_keys=True,
#                   delimiter='\t', 
#                   merged_value_delimiter=',')

###############################################################################
############
# Function: 
############
def fetch_uniprot_from_db(refseq_id, cursor):
    """
    Fetch the UniProt ID associated with a given RefSeq ID from the database.
    
    This function executes a SQL query to retrieve the UniProt accession number
    for a specified RefSeq ID. It handles multiple possible outcomes: no matching UniProt ID,
    a single UniProt ID, or multiple UniProt IDs. In cases where multiple IDs are found,
    it logs a message and returns the first one.
    
    Parameters:
    - refseq_id (str): The RefSeq ID for which the UniProt ID is to be fetched.
    - cursor (Cursor): A database cursor object used to execute the query.
    
    Returns:
    - tuple:
        - str: The UniProt ID if found, 'N/A*' if no match is found, 'error' if an exception occurs.
        - bool: True if multiple UniProt IDs were found for the RefSeq ID, False otherwise.
    
    Exceptions:
    - Raises an exception and prints an error message if the database query fails.
    
    Usage Example:
    >>> cursor = db_connection.cursor()
    >>> fetch_uniprot_from_db('NP_00001', cursor)
    ('P12345', False)  # Example of a single match
    
    It is primarly used by the process_biogrid_file() to get the ids
    """
    
    try:
        #query = f"SELECT uniprot_acc FROM uniprot_refseq WHERE refseq_acc = '{refseq_id}'"
        query = f"SELECT uniprot_acc FROM uniprot_refseq WHERE refseq_acc LIKE '{refseq_id}%'"
        cursor.execute(query)
        results = cursor.fetchall()
        if not results:
            print(f"No UniProt ID found for RefSeq ID: {refseq_id}")
            #return 'not_in_pv_db', False
            return 'N/A*', False
        elif len(results) > 1:
            print(f"Multiple UniProt IDs found for RefSeq ID: {refseq_id}, selecting the first one.")
            return results[0][0], True  # Return the first result but indicate multiple were found.
        return results[0][0], False
    except Exception as e:
        print(f"Error fetching UniProt ID for RefSeq ID {refseq_id}: {e}")
        return 'error', False  # Return 'error' and False if an exception occurs.
###############################################################################
############
# Function: 
############
def process_biogrid_file(input_file, db_params):
    """Processes the BioGRID TSV file to update missing UniProt IDs using a PostgreSQL database with debug info."""
    
    """
    Processes a BioGRID TSV file to update missing UniProt IDs based on corresponding RefSeq IDs
    using a PostgreSQL database. It also collects debug information for each row processed.

    This function reads a TSV file containing biogrid interaction data, checks for missing UniProt IDs,
    and attempts to fetch and update these IDs from a database using RefSeq IDs. It keeps track of
    any rows where updates occur and whether multiple potential UniProt IDs were found for single RefSeq IDs.

    Parameters:
    - input_file (str): The path to the BioGRID TSV file to be processed.
    - db_params (dict): A dictionary containing database connection parameters such as host, database name,
                        user, and password. This is used to establish a connection to the PostgreSQL database.

    Returns:
    - tuple:
        - list: A list of dictionaries, each representing a row from the input file with updated UniProt IDs.
        - dict: A dictionary where keys are row indices (1-indexed) and values are lists of strings describing
                the updates or actions taken for each row, e.g., 'multiple_uniprot_for_refseq1'.

    Processes the file by:
    - Connecting to the PostgreSQL database using the provided db_params.
    - Iterating over each row in the TSV file.
    - For each row, checking if the UniProt IDs (uniprot_id1 and uniprot_id2) are missing and, if so,
      attempting to fetch corresponding UniProt IDs using RefSeq IDs.
    - Collecting debug information about which rows had multiple potential matches or were updated.
    - Closing the database connection and returning the processed data along with debug info.

    Usage:
    >>> db_params = {'dbname': 'bioinfo_db', 'user': 'dbuser', 'password': 'secret', 'host': 'localhost'}
    >>> data, debug_info = process_biogrid_file('biogrid_data.tsv', db_params)
    """
    
    conn = psycopg2.connect(**db_params)
    cursor = conn.cursor()

    data = []
    debug_info = defaultdict(list)

    with open(input_file, 'r', newline = '') as file:
        reader = csv.DictReader(file, delimiter = '\t')
        for idx, row in enumerate(reader):
            #updated = False
            #print(f"Processing row {idx + 1}")
            
            if row['uniprot_id1'] == 'N/A' and row['refseq_id1'] != 'N/A':
                new_id, multiple = fetch_uniprot_from_db(row['refseq_id1'], cursor)
                if multiple:
                    debug_info[idx + 1].append('multiple_uniprot_for_refseq1')
                row['uniprot_id1'] = new_id
                #updated = True
                debug_info[idx + 1].append('up1_updated')


            if row['uniprot_id2'] == 'N/A' and row['refseq_id2'] != 'N/A':
                new_id, multiple = fetch_uniprot_from_db(row['refseq_id2'], cursor)
                if multiple:
                    debug_info[idx + 1].append('multiple_uniprot_for_refseq2')
                row['uniprot_id2'] = new_id
                #updated = True
                debug_info[idx + 1].append('up1_updated')

            # if updated:
            #     print(f"Updated UniProt IDs in row {idx + 1}")
            #     debug_info[idx + 1].append('updated')
            # row['interaction_id'] = f"{row['uniprot_id1']}_{row['uniprot_id2']}"
            # data.append(row)
            row['interaction_id'] = f"{row['uniprot_id1']}_{row['uniprot_id2']}"
            data.append(row)

    cursor.close()
    conn.close()
    return data, dict(debug_info)
###############################################################################