import os
import sys
import argparse
import pandas as pd
import glob


# requires python >= 3.6
# run: python copy_samplefiles.py haptophyta.txt --source ../../pep/ --destination ../../haptophyta_pep

#quant files are haptophyta_out/quant/MMETSP0006_SRR1296917_x_haptophyta_MMETSP0006.salmon/quant.sf
#namemaps are namemaps/MMETSP0503.trinity_out_2.2.0.Trinity.fasta.dammit.namemap.csv
# orthogroups file is haptophyta_orthofinder_results/Results_june_2019/Orthogroups/Orthogroups.tsv
# in header, <samplename>.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep --> remove trailing info, not needed

# because we need to read in each namemap, it'd be better to work with a single sample at a time. so iterate through the COLUMNS of the
# orthogroup table, rather than the rows (each orthogroup)

def build_counttable(orthogroups, quant_dir, og_counttable, samplelist=None, namemap_dir=None):
    # read samples -> sample list
    # reading entire orthogroup table -->  pandas
    if samplelist:
        with open(samplelist) as f:
            samples = [line.strip() for line in f]
    orthogroups = pd.read_csv(orthogroups, sep = '\t', header= 0, dtype=str ) #index_col=0) # read orthogroup samples
    orthogroups.columns =  orthogroups.columns.str.split('.').str[0]

    for sample in orthogroups.columns:
        if samplelist:
            if sample not in samples:
                continue
        orthogroups[sample] = orthogroups[sample].str.split('|').str[0]
        sample_og_dict = dict(zip(orthogroups[sample], orthogroups['Orthogroup']))
        # dammit: orthogroup Dictionary
    # need to build trinity name: orthogroup dictionary for each sample
        if namemap_dir:
            namemap = glob.glob(os.path.join(namemap_dir, sample + '*.csv'))
            dammit2trin = pd.read_csv(namemap[0], sep=',', header=0, dtype=str)
            dammit2trin['original'] = dammit2trin['original'].str.split(' ').str[0]
            #trinity to dammit map!
            dam2trinD = dict(zip(dammit2trin['renamed'], dammit2trin['original']))
            trinOG = {}
            for dam, trin in dam2trinD.items():
                og = sample_og_dict.get(dam, None)
                if og:
                    trinOG[trin] = og
        else:
            # if no namemaps, assume we had trinity names to start!
            trinOG = sample_og_dict

           # trinOG = trinity to orthogroup map! --> now let's grab quant files!
        try:
            quant_file = glob.glob(os.path.join(quant_dir, sample + '*', "quant.sf"))[0]
        except Exception as e:
             sys.stderr.write(f"\n\tWarning: cannot find quant file for sample {sample}, dropping from count table. \n\n")
             orthogroups.drop(columns =sample, inplace=True)
             continue

        quant = pd.read_csv(quant_file, header =0, sep = '\t')
        quant['Name'] = quant['Name'].str.rsplit('-').str[2] # get trinity name
        trinCounts = dict(zip(quant['Name'], quant['TPM']))
        #trinCounts = dict(zip(quant['Name'], quant['NumReads']))
        #if not countD:
        countD = {}
        for trin, OG in trinOG.items():
            countD[OG] = trinCounts.get(trin, 0)

        #just overwrite old orthogroups dictionary
        orthogroups[sample] = orthogroups['Orthogroup'].map(countD)


    orthogroups.to_csv(path_or_buf=og_counttable, sep = '\t', index=False)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('orthogroups', default = "Orthogroups.tsv")
    p.add_argument('--sample_list')
    p.add_argument('--namemap_directory')
    p.add_argument('--quant_directory', required=True)
    p.add_argument('--orthogroup_counttable', required=True)
    args = p.parse_args()
    sys.exit(build_counttable(args.orthogroups, args.quant_directory, args.orthogroup_counttable, args.sample_list, args.namemap_directory))
