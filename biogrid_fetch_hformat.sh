#!/bin/bash

########################################################################
# Script to process BioGRID data file for human interactions
    # 1. Downloads the current release from: https://downloads.thebiogrid.org/File/BioGRID/Release-Archive/BIOGRID-4.4.240/BIOGRID-ALL-4.4.240.mitab.zip

    # 2. Modifies the header to replace whitespace with underscore
    # 3. Extracts data related to human interactions i.e. taxid: 9606

# Override defaults using args
# Usage: ./script_name.sh -u <URL> -z <FILE> -o <OUTPUT_FILE>
########################################################################
URL='https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-4.4.240'
FILE='BIOGRID-ALL-4.4.240.mitab.zip'
OUTPUT_FILE='biogrid_human_interactions.txt'


while getopts u:z:o: flag
do
    case "${flag}" in
        u) URL=${OPTARG};;
        z) FILE=${OPTARG};;
        o) OUTPUT_FILE=${OPTARG};;
        *) echo "Invalid option"; exit 1;;
    esac
done

DATA_FILE=$(basename -s .zip ${FILE}).txt

if [ -z "${URL}" ] || [ -z "${FILE}" ] || [ -z "${OUTPUT_FILE}" ]; then
    echo "Usage: $0 -u <URL> -z <FILE> -o <OUTPUT_FILE>"
    exit 1
fi

echo "Downloading BioGRID data from ${URL}/${FILE}"
wget -c "${URL}/${FILE}"

echo "Unzipping ${FILE}..."
unzip -o "${FILE}"

echo "Modifying header and extracting human interactions..."
head -1 "${DATA_FILE}" | sed -E 's/ /_/g; s/#//g' > "${OUTPUT_FILE}"
awk -F'\t' '($10 ~ /taxid:9606/ || $11 ~ /taxid:9606/)' "${DATA_FILE}" >> "${OUTPUT_FILE}"

echo "Results saved in ${OUTPUT_FILE}"
echo "Cleaning up..."
rm "${FILE}"
echo "Done."
