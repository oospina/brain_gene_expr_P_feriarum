# This script takes a list of sample IDs (eg. 'I' numbers) and returns sums of raw counts from files
# containing two columns: transcript ID and raw count. eg. the result of cut -f1,13 on coverageBed output files, or from Alan's java script
# The script must be in the same folder with files to be sumarized, and must be ran from there
# Use as python summarizeRawCounts.py [file-with-sample-IDs]

import os
import sys
import pandas as pd

# Takes the name of the file with sample IDs
arg = sys.argv[1]
# Get the current directory and sets it as working dir 
path = os.getcwd()
os.chdir(path)

# open sample IDs file and reads its lines
samples = open(arg, "r")
ids = samples.read().splitlines()

fname = ids[0] + '_RawCountsFromBAM_genes.csv'
df_all  = pd.read_csv(fname, names=['RefTrans', str(ids[0])])

for i in range(1, len(ids)):
	fname = ids[i] + '_RawCountsFromBAM_genes.csv' 
	df = pd.read_csv(fname, names=['RefTrans', str(ids[i])])
	df_all = pd.merge(df_all, df, on="RefTrans", how="outer")

df_all.to_csv("sumRawCounts.csv", index=False, float_format='%.0f')
