"""
extended/modified from biobombe scripts https://github.com/greenelab/BioBombe, by N.T.Pierce
and https://github.com/greenelab/nf1_inactivation/blob/master/scripts/process_rnaseq.py
"""

import os
import sys
import random
import argparse
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn import preprocessing

random.seed(42)

#wget https://osf.io/ek9nu/download -O haptophyta_orthogroup.quant.tsv

def read_counts(countfile):
    if '.tsv' in countfile or '.csv' in countfile:
        separator = '\t'
        if '.csv' in countfile:
            separator = ','
        try:
            counts = pd.read_csv(countfile, dtype=str, sep=separator, index_col=0)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {countfile} file is not properly formatted. Please fix.\n\n")
            print(e)
    elif '.xls' in countfile:
        try:
            counts = pd.read_excel(countfile, dtype=str, index_col=0)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {countfile} file is not properly formatted. Please fix.\n\n")
            print(e)
    return counts

def preprocess_data(countfile, out_file=None, scale = "min_max", percTest = 0.1, mad= False, num_mad_genes = 5000, out_folder = ""):
    """
    Zero-one scale the expression data

    Output:
    Writes normalized matrix (if output=True) and mad genes to file
    (if mad=True); returns the normalized matrix if output=False
    """

    #################
    # Preprocessing #
    #################

    # load data to pd dataframe
    #expr_data = pd.read_csv(countfile, dtype='str', sep='\t', index_col=0)
    expr_data = read_counts(countfile)

    # Drop all row names with unidentified gene name
    #expr_data = data[-data.index.str.contains('?', regex=False)]

    expr_data = expr_data.sort_index() # sort by gene name

    if scale == "min_max":
       # Zero-one normalize
        min_max_scaler = preprocessing.MinMaxScaler()
        expr_scaled = min_max_scaler.fit_transform(expr_data.T)
    elif scale == 'zscore':
        expr_scaled = preprocessing.scale(expr_data.T, axis=0)

    expr_norm = pd.DataFrame(expr_scaled, index=expr_data.columns, columns=expr_data.index).T  # transform back

    # check that output folder exists, make it if not
    if out_folder:
        if not os.path.exists(out_folder):
            os.mkdirs(out_folder)

    # write scaled output

    if not out_file:
        out_file = os.path.join(out_folder, os.path.basename(countfile).rsplit('.')[0] + '.scaled.tsv')

    expr_norm.to_csv(out_file, sep='\t', header=True, index=True)
    #print(min_max_scaler.data_max_)

    ############################
    # Split Test, Training sets
    ############################

    trainDF, testDF = train_test_split(expr_norm, test_size=percTest, random_state =42) #, shuffle=False

    # write training set to file
    train_file = os.path.join(out_folder, os.path.basename(countfile).rsplit('.')[0] + '.processed.train.tsv.gz')
    trainDF.to_csv(train_file, sep='\t', compression='gzip', float_format='%.3g')

    # write test set to file
    test_file = os.path.join(out_folder, os.path.basename(countfile).rsplit('.')[0] + '.processed.test.tsv.gz')
    testDF.to_csv(test_file, sep='\t', compression='gzip', float_format='%.3g')

    ###########################################
    # Sort on Median Absolute Deviation "mad" #
    ###########################################

    ### THIS CURRENTLY ONLY DOES MAD FOR TRAINING SET, NOT ALL GENES, as per BioBombe notebooks
    if mad:

        mad_file = os.path.join(out_folder, os.path.basename(countfile).rsplit('.')[0] + '.processed.train.mad.' + str(num_mad_genes) + '.tsv.gz')
        # sort
        mad_genes_df = pd.DataFrame(trainDF.mad(axis=0).sort_values(ascending=False)).reset_index()
        mad_genes_df.columns = ['gene_id', 'median_absolute_deviation']
        mad_genes_df.to_csv(mad_file, sep='\t', index=False)
        #subset to only top genes
        #top_mad_genes = mad_genes_df.iloc[0:num_mad_genes, ].index
        # subset original dataset
        #expr_subset_df = expr_norm.loc[:, top_mad_genes]
        # Write to file
        #expr_subset.to_csv(mad_file, sep='\t', index=False)


if __name__ == '__main__':
    #countfile = "haptophyta_orthogroup.quant.tsv"

    p = argparse.ArgumentParser()
    p.add_argument("countfile", help = "csv,tsv,or xls rnaseq gene expression table")
    p.add_argument("-s", "--scale",help="method to scale expression matrix", default='min_max')
    p.add_argument("--percTest", help="percent of data to set aside as test set. Training set will be 1-percTest", default = 0.1, type=float)
    p.add_argument("--mad", help="subset mad genes", action="store_true", default=False)
    p.add_argument("--output_folder", default="")
    p.add_argument("--output_filename", default=None)
    p.add_argument("--num_mad_genes", help="number of highest median absolute deviation genes to output", type=int, default=5000)
    args = p.parse_args()
    sys.exit(preprocess_data(args.countfile, args.output_filename, args.scale, args.percTest, args.mad, args.num_mad_genes, args.output_folder))
