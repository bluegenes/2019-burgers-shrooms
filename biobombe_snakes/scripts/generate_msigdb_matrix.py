# # Processing MSigDB Gene Sets into Binary Matrix
# # Modified from https://github.com/greenelab/BioBombe/blob/master/3.build-hetnets/generate_msigdb_matrix.ipynb
#
# This script loads the full MSigDB gene set `.gmt` file (version 6.1) and outputs a binary, gene by gene set matrix indicating gene membership in the given gene set.
#
# **Note that we exclude gene sets with restrictive licences (KEGG, Biocarta, and The AAAS/STKE Cell Signaling Database)**

import os
import sys
import csv
import argparse
import numpy as np
import pandas as pd


def make_template_matrix(msigdb_file, checkblacklist=True): #blacklist = ['KEGG', 'BIOCARTA', 'ST_'], checkblacklist=True):
    """
    Retrieve all genes and pathways from given msigdb .gmt file

    Output:
    sorted gene by pathways pandas dataframe. Entries indicate membership
    """
    all_db_pathways = []
    all_db_genes = []
    # Resources with restrictive licenses
    blacklist = ('KEGG', 'BIOCARTA', 'ST_')

    # Get a set of all genes and all pathways in MSigDB (not blacklisted)
    with open(msigdb_file, 'r') as msigdb_fh:
        msigdb_reader = csv.reader(msigdb_fh, delimiter='\t')

        for row in msigdb_reader:
            signature_name = row[0]
            signature_genes = row[2:]

            if checkblacklist:
                if signature_name.startswith(blacklist):
                    continue

            all_db_pathways.append(signature_name)
            all_db_genes += signature_genes

    big_msigdb_df = pd.DataFrame(0, index=set(all_db_genes), columns=all_db_pathways)
    big_msigdb_df = big_msigdb_df.sort_index()
    big_msigdb_df = big_msigdb_df.T.sort_index().T

    # Loop through file again to populate dataframe. This is a fast implementation
    with open(msigdb_file, 'r') as msigdb_fh:
        msigdb_reader = csv.reader(msigdb_fh, delimiter='\t')
        for row in msigdb_reader:
            signature_name = row[0]
            signature_genes = row[2:]
            if checkblacklist:
                if signature_name.startswith(blacklist):
                    continue

            for gene in signature_genes:
                big_msigdb_df.at[gene, signature_name] = 1

    return big_msigdb_df



# Store .gmt files
#full_msigdb_file = os.path.join('data', 'msigdb.v6.1.entrez.gmt')

# ## Process MSigDB gmt files into large matrix

def write_msigdb_file(msigdb_file, out_binary_matrix, checkblacklist=True):
    full_msigdb_df = make_template_matrix(msigdb_file, checkblacklist=True)
    print(full_msigdb_df.shape)
    full_msigdb_df.to_csv(out_binary_matrix, sep='\t', compression='bz2')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--input_msigdb_file', help="input_full_msigdb_file", default= "msigdb/msigdb.v6.1.entrez.gmt")
    p.add_argument('--out_matrix', help="output_matrix", default= "msigdb/full_msigdb_binary_matrix.tsv.bz2")
    p.add_argument('--checkblacklist', action="store_true", help='check the blacklist?', default=False)
    args = p.parse_args()
    sys.exit(write_msigdb_matrix(args.input_msigdb_file, args.out_matrix, args.checkblacklist))

