import os
import sys
import pandas as pd


def read_params(pfile):
# read_params has no csv option, because some paramsfiles have comma separated values within the params columns
    if '.tsv' in pfile:
        separator = '\t'
        try:
            paramsDF = pd.read_csv(pfile, sep=separator, index_col=0, compression= "infer")
        except Exception as e:
            sys.stderr.write(f"\n\tError: {pfile} file is not properly formatted. Please fix.\n\n")
            print(e)
    elif '.xls' in pfile:
        try:
            paramsDF = pd.read_excel(pfile, index_col=0, compression= "infer")
        except Exception as e:
            sys.stderr.write(f"\n\tError: {pfile} file is not properly formatted. Please fix.\n\n")
            print(e)
    return paramsDF

def read_counts_or_params(countfile):
    # allow csv, tsv, or xls/xlsx
    if '.tsv' in countfile or '.csv' in countfile:
        separator = '\t'
        if '.csv' in countfile:
            separator = ','
        try:
            counts = pd.read_csv(countfile, sep=separator, index_col=0, compression= "infer")
        except Exception as e:
            sys.stderr.write(f"\n\tError: {countfile} file is not properly formatted. Please fix.\n\n")
            print(e)
    elif '.xls' in countfile:
        try:
            counts = pd.read_excel(countfile, index_col=0, compression= "infer")
        except Exception as e:
            sys.stderr.write(f"\n\tError: {countfile} file is not properly formatted. Please fix.\n\n")
            print(e)
    return counts
