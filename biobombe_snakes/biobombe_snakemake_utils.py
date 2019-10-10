import os
import sys
import pandas as pd

def read_params(pfile):
    if '.tsv' in pfile:
        separator = '\t'
        try:
            paramsDF = pd.read_csv(pfile, dtype=str, sep=separator, index_col=0)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {pfile} file is not properly formatted. Please fix.\n\n")
            print(e)
    elif '.xls' in countfile:
        try:
            paramsDF = pd.read_excel(pfile, dtype=str, index_col=0)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {pfile} file is not properly formatted. Please fix.\n\n")
            print(e)
    return paramsDF
