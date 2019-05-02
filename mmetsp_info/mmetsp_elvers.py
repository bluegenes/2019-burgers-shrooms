#! /usr/bin/env python
"""
Auto generate mmetsp elvers files from csv
"""
import os
import sys
import argparse
import pandas as pd



def build_fq2(row):
    base_link = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"
    SRR = row['SRR']
    if row['LibraryLayout'] == "PAIRED":
        fq2 = base_link + SRR[0:6] + '/00' + SRR[-1] + '/' + SRR  + '/' + SRR + '_2.fastq.gz'
    return fq2


def build_elvers_samples(samples_file, out_csv, subset_list=None):
    if '.tsv' in samples_file or '.csv' in samples_file:
        separator = '\t'
        if '.csv' in samples_file:
            separator = ','
        try:
            samples = pd.read_csv(samples_file, dtype=str, sep=separator)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {samples_file} file is not properly formatted. Please fix.\n\n")
            print(e)
    elif '.xls' in samples_file:
        try:
            samples = pd.read_excel(samples_file, dtype=str)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {samples_file} file is not properly formatted. Please fix.\n\n")
            print(e)

    base_link = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"
    if "SRR" in samples.columns and "LibraryLayout" in samples.columns:
        samples['fq1'] = samples['SRR'].apply(lambda x : base_link + x[0:6] + '/00' + x[-1] + '/' + x + '/' + x + '_1.fastq.gz')
        samples['fq2'] = samples.apply(lambda row : build_fq2(row), axis=1)

    elvers_samples = samples.loc[:, ['SampleName','SRR','fq1','fq2']]
    elvers_samples.rename(columns={"SampleName": "sample", "SRR": "unit"}, inplace=True)

    # for now, just keep this as a large csv file. We can do all preprocessing using a single file.
    elvers_samples.to_csv(out_csv,  sep=',', index=False)


#def per_sample_elvers():
# write a per-sample csv and yaml file!


#def subset_samples_file(samples_df, subset_names)
    #if subset_list:
    #    sub_samples = subset_samples_file(samples, subset_list)
    #    sub_samples.to_csv(out_csv, sep=', ')
    #else:




if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('samples_file')
    p.add_argument('-o', '--out_csv')
    p.add_argument('--subset_list')
    args = p.parse_args()
    sys.exit(build_elvers_samples(args.samples_file, args.out_csv, args.subset_list))
