#!/usr/bin/env bash

# The URL to the tar file
RNA_DATAURL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179994/suppl/GSE179994_all.Tcell.rawCounts.rds.gz"
TCR_DATAURL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE179nnn/GSE179994/suppl/GSE179994_all.scTCR.tsv.gz"
DATADIR="prepared-data"
RNA_FILE="$DATADIR/GSE179994_all.Tcell.rawCounts.rds.gz"
TCR_FILE="$DATADIR/GSE179994_all.scTCR.tsv.gz"

echo "- Make the directory for prepared data ..."
mkdir -p prepared-data

echo "- Download the RNA data if needed ..."
if [ ! -e $RNA_FILE ]; then
    wget -O $RNA_FILE $RNA_DATAURL
fi

echo "- Download the TCR data if needed ..."
if [ ! -e $TCR_FILE ]; then
    wget -O $TCR_FILE $TCR_DATAURL
fi

echo "- Prepare the RNA data ..."
Rscript prepare-rna.R

echo "- Prepare the TCR data ..."
Rscript prepare-tcr.R

echo "- Done!"
