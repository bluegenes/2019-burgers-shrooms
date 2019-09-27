import os
import sys
import argparse
import pandas as pd
import glob


# requires python >= 3.6
#quant files are haptophyta_out/quant/MMETSP0006_SRR1296917_x_haptophyta_MMETSP0006.salmon/quant.sf
#namemaps are namemaps/MMETSP0503.trinity_out_2.2.0.Trinity.fasta.dammit.namemap.csv
# orthogroups file is haptophyta_orthofinder_results/Results_june_2019/Orthogroups/Orthogroups.tsv
# in header, <samplename>.trinity_out_2.2.0.Trinity.fasta.transdecoder.pep --> remove trailing info, not needed


def build_gtmap(orthogroups_file, gtmap_dir, samplelist=None, namemap_dir=None):
    # reading entire orthogroup table -->  pandas
    if not os.path.exists(gtmap_dir):
        os.mkdir(gtmap_dir)

    if samplelist:
        # optional sample list to care about -> sample list
        with open(samplelist) as f:
            samples = [line.strip() for line in f]
    #orthogroups = pd.read_csv(orthogroups, sep = '\t', header= 0, dtype=str ) #index_col=0) # read orthogroup samples
    orthogroups = pd.read_csv(orthogroups_file, sep = '\t', header= 0, dtype=str, index_col=0)
    orthogroups.columns =  orthogroups.columns.str.split('.').str[0] # split off the excess, leaving just the MMETSP id

    for sample in orthogroups.columns:
        print(sample)
        if samplelist: # if we have a subset to go through, check that we care about this sample
            if sample not in samples:
                continue

        ### let's not be fancy, just iterate through, build dictionaries along the way
        temp = orthogroups[sample].dropna()
        # dammit: orthogroup Dictionary
        sample_og_dict = {}
        for og, tx_str in temp.items():
            tx_list = tx_str.split(" ")
            for tx in tx_list:
                if all([tx.startswith("Transcript_"), "+" not in tx, "-" not in tx]):
                     tx = tx.split("|")[0] # remove orf identifiers
                     sample_og_dict[tx] = og # dictionary of transcript: orthogroup

        ### OPTIONAL - keep track of all the sample_og_dicts
        #orthogroup_dicts[sample] = sample_og_dict # dict of dicts to save the temp d. But we can probably just iterate through this sample rn...
        # grab the dataframe for this sample, expand multiple transcripts
        #orthogroups[sample] = orthogroups[sample].str.split(' ')
        #orthogroups[sample] = orthogroups[sample].str.split('|').str[0] # this is for a single transcript. Instead, there are many. Need to do differently
        #sample_og_dict = dict(zip(orthogroups[sample], orthogroups['Orthogroup']))

    # need to build trinity name: orthogroup dictionary for each sample

        # for tximport, we want gene /t trans aka orthogroup /t trinity_name
        if namemap_dir:
            namemap = glob.glob(os.path.join(namemap_dir, sample + '*.csv'))
            dammit2trin = pd.read_csv(namemap[0], sep=',', header=0, dtype=str)
            dammit2trin['original'] = dammit2trin['original'].str.split(' ').str[0]
            #trinity to dammit map!
            dam2trinD = dict(zip(dammit2trin['renamed'], dammit2trin['original']))
            orthogroup_gene_trans_map_file = os.path.join(gtmap_dir, sample + "_orthogroup_gtmap.tsv")
            with open(orthogroup_gene_trans_map_file, 'w') as og_map:
            #for dam, trin in dam2trinD.items():
                for dam, trin in dam2trinD.items():
                    og = sample_og_dict.get(dam, None)
                    if og:
           #        trinOG[trin] = og
                        og_map.write(f"{og}\t{trin}\n")
            print('done!')
        else:
            # if no namemaps, assume we had trinity names to start!
            #trinOG = sample_og_dict
            orthogroup_gene_trans_map_file = os.path.join(gtmap_dir, sample + "_orthogroup_gtmap.tsv")
            with open(orthogroup_gene_trans_map_file, 'w') as og_map:
                for trin, og in sample_og_dict.items():
                    og_map.write(f"{og}\t{trin}\n")

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('orthogroups', default = "Orthogroups.tsv")
    p.add_argument('--sample_list')
    p.add_argument('--namemap_directory')
    p.add_argument('--gtmap_directory', default= "orthogroup_gtmaps")
    args = p.parse_args()
    sys.exit(build_gtmap(args.orthogroups, args.gtmap_directory, args.sample_list, args.namemap_directory))
