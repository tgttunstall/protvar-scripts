#!/bin/bash

########################################################################
# Script to process BioGRID data file for human interactions
    # 1. Downloads the current release from: https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-4.4.240/BIOGRID-ALL-4.4.240.mitab.zip

    # 2. Unzips the downloaded file
    # 3. Modifies the header to remove '#' and replace whitespace with '_'
    # 4. Extracts data related to human interactions i.e. taxid: 9606 from fields 10 and 11
        # so extracting either fields with taxid:9606
    # 5. Cleans up temporary files based on specified options

# Override defaults using args
# Usage: ./script_name.sh -u <URL> -f <FILE> -o <OUTPUT_FILE> [--clean[=ext1,ext2,...]] [--force] [--time]

# Options:
#   -u <URL>       : URL to download the BioGRID data file (optional, default provided)
#   -f <FILE>      : Name of the input file (optional, default provided)
#   -o <OUTPUT_FILE>: Name of the output file (required)
#   --clean=[ext1,ext2,...]: Clean up input files. If extensions are specified, only those will be removed. This gives the option to clean the zip file, txt file, etc. 
#   --force        : Force re-download of the data file even if it exists
#   --time         : Time the execution of major operations
#
# Example usage:
#   ./script_name.sh -u https://example.com/data.zip -f data.zip -o output.txt --clean=zip,txt --force --time

########################################################################
URL='https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.240'
FILE='BIOGRID-ALL-4.4.240.mitab.zip'
OUTPUT_FILE='biogrid_human_interactions.txt'

CLEAN=false
FORCE=false
TIME=false
CLEAN_EXTENSIONS=()

while getopts u:f:o:-: flag
do
    case "${flag}" in
        u) URL=${OPTARG};;
        f) FILE=${OPTARG};;
        o) OUTPUT_FILE=${OPTARG};;
        -) case "${OPTARG}" in
               clean) CLEAN=true;;
               clean=*) CLEAN=true; IFS=',' read -ra CLEAN_EXTENSIONS <<< "${OPTARG#*=}";;
               force) FORCE=true;;
               time) TIME=true;;
               *) echo "Invalid option: --${OPTARG}" >&2; exit 1;;
           esac;;
        *) echo "Invalid option: -${flag}" >&2; exit 1;;
    esac
done

# Print the options script is run with whether default or custom
echo "Script running with the following:"
echo "URL: $URL"
echo "Input file: $FILE"
echo "Output file: $OUTPUT_FILE"


DATA_FILE=$(basename -s .zip ${FILE}).txt

time_cmd() {
    if [ "$TIME" = true ]; then
        time "$@"
    else
        "$@"
    fi
}

echo "=================================================================="

# Check if zip and unzipped files exist, then start processing...
if [ "$FORCE" = true ]; then
    echo "Force option used. Downloading and processing..."
    time_cmd wget -c "${URL}/${FILE}" || { echo "Download failed" >&2; exit 1; }
    time_cmd unzip -o "${FILE}" || { echo "Unzip failed" >&2; exit 1; }
elif [ -f "${FILE}" ] && [ -f "${DATA_FILE}" ]; then
    echo "Both zip and unzipped files exist. Starting processing..."
elif [ -f "${FILE}" ]; then
    echo "Zip file exists. Unzipping and starting processing..."
    time_cmd unzip -o "${FILE}" || { echo "Unzip failed" >&2; exit 1; }
elif [ -f "${DATA_FILE}" ]; then
    echo "Unzipped file exists. Starting processing..."
else
    echo "Neither zip nor unzipped file exists. Downloading..."
    time_cmd wget -c "${URL}/${FILE}" || { echo "Download failed" >&2; exit 1; }
    time_cmd unzip -o "${FILE}" || { echo "Unzip failed" >&2; exit 1; }
fi

# Processing Header
echo "No. of lines in unzipped file: $(wc -l ${DATA_FILE})"
echo "Modifying header: removing '#' and replacing 'whitespace' with '_'"
time_cmd sed -E '1s/^#//; 1s/ /_/g; 1q' "${DATA_FILE}" > "${OUTPUT_FILE}" || { echo "Header modification failed" >&2; exit 1; }

# Extracting human interactions using taxid:9606 from fields 10 and 11
echo "Extracting human interactions i.e. taxid:9606..."
time_cmd awk -F'\t' '($10 ~ /taxid:9606/ || $11 ~ /taxid:9606/)' "${DATA_FILE}" >> "${OUTPUT_FILE}" || { echo "Extraction of human interactions failed" >&2; exit 1; }

echo "Results saved in ${OUTPUT_FILE}"
echo "No. of lines post extraction: $(wc -l ${OUTPUT_FILE})"
echo "=================================================================="
# Clean up 
# Define allowed extensions
ALLOWED_EXTENSIONS=("zip" "txt" "mitab")

echo "Cleaning up..."
if [ "$CLEAN" = true ]; then
    if [ ${#CLEAN_EXTENSIONS[@]} -eq 0 ]; then
        rm "${DATA_FILE}" "${FILE}" || echo "Failed to remove some files" >&2
        echo "Removed input files"
    else
        for ext in "${CLEAN_EXTENSIONS[@]}"; do
            if [[ " ${ALLOWED_EXTENSIONS[@]} " =~ " ${ext} " ]]; then
                if [ -f "${FILE%.*}.${ext}" ]; then
                    rm "${FILE%.*}.${ext}" && echo "Removed ${FILE%.*}.${ext}"
                else
                    echo "File ${FILE%.*}.${ext} not found, skipping."
                fi
            else
                echo "Warning: Extension '${ext}' is not in the allowed list and will be ignored." 
                echo "Valid extensions for clean up are: ${ALLOWED_EXTENSIONS[*]}"
            fi
        done
    fi
else
    echo "No cleanup performed. Use --clean option to remove input files."
fi
