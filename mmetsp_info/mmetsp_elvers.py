#! /usr/bin/env python
"""
generate mmetsp elvers file(s) from csv
"""
import os
import sys
import argparse
import pandas as pd



def build_fq2(row):
    base_link = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/"
    #SRR = row['SRR']
    SRR = row['Run.x']
    if row['LibraryLayout'] == "PAIRED":
        fq2 = base_link + SRR[0:6] + '/00' + SRR[-1] + '/' + SRR  + '/' + SRR + '_2.fastq.gz'
    return fq2


def build_elvers_samples(samples_file, out_csv, subset_list=None, exclude_list=None):
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
    #if "SRR" in samples.columns and "LibraryLayout" in samples.columns:
    if "Run.x" in samples.columns and "LibraryLayout" in samples.columns:
        #samples['fq1'] = samples['SRR'].apply(lambda x : base_link + x[0:6] + '/00' + x[-1] + '/' + x + '/' + x + '_1.fastq.gz')
        samples['fq1'] = samples['Run.x'].apply(lambda x : base_link + x[0:6] + '/00' + x[-1] + '/' + x + '/' + x + '_1.fastq.gz')
        samples['fq2'] = samples.apply(lambda row : build_fq2(row), axis=1)
        samples['reference'] = samples['SampleName'].apply(lambda x :  'MMETSP_assemblies_figshare/' + x + '.trinity_out_2.2.0.Trinity.fasta.renamed.fasta')

    elvers_samples = samples.loc[:, ['SampleName','Run.x','fq1','fq2', 'reference']]
    elvers_samples.rename(columns={"SampleName": "sample", "Run.x": "unit"}, inplace=True)

    subset, exclude = [],[]
    if subset_list:
        with open(subset_list) as subsetFile:
            subset = [line.strip() for line in subsetFile]
            elvers_samples = elvers_samples[elvers_samples['sample'].isin(subset)]
    if exclude_list:
        with open(exclude_list) as excludeFile:
            exclude = [line.strip() for line in excludeFile]
            elvers_samples = elvers_samples[~elvers_samples['sample'].isin(exclude)]

    # for now, just keep this as a large csv file. We can do all preprocessing using a single file.
    elvers_samples.to_csv(out_csv,  sep=',', index=False)




if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('samples_file')
    p.add_argument('-o', '--out_csv')
    p.add_argument('--subset_list')
    p.add_argument('--exclude_list')
    args = p.parse_args()
    sys.exit(build_elvers_samples(args.samples_file, args.out_csv, args.subset_list, args.exclude_list))
