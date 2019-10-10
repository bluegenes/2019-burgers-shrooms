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

def build_counttable(orthogroups_file, gtmap_dir, samplelist=None, namemap_dir=None, quant_dir=None, verbose=False):
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
        if verbose:
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

    # if we started with dammit names, need to build trinity name: orthogroup dictionary for each sample
        trinOG = {}
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
                        trinOG[trin] = og
                        # for tximport, we want gene /t trans aka orthogroup /t trinity_name
                        og_map.write(f"{og}\t{trin}\n")
        else:
            # if no namemaps, assume we had trinity names to start!
            #trinOG = sample_og_dict
            orthogroup_gene_trans_map_file = os.path.join(gtmap_dir, sample + "_orthogroup_gtmap.tsv")
            with open(orthogroup_gene_trans_map_file, 'w') as og_map:
                for trin, og in sample_og_dict.items():
                    trinOG[trin] = og
                    og_map.write(f"{og}\t{trin}\n")

        # trinOG = trinity to orthogroup map! --> now let's grab quant files!

        if quant_dir:
            try:
                quant_file = glob.glob(os.path.join(quant_dir, sample + '*', "quant.sf"))[0]
                quant = pd.read_csv(quant_file, header =0, sep = '\t')
            except Exception as e:
                sys.stderr.write(f"\n\tWarning: cannot find quant file for sample {sample}, dropping from count table. \n\n")
                orthogroups.drop(columns =sample, inplace=True)
                continue
            quant['Name'] = quant['Name'].str.rsplit('-').str[2] # get trinity name
            trinCounts = dict(zip(quant['Name'], quant['TPM']))
            #trinCounts = dict(zip(quant['Name'], quant['NumReads']))
            #if not countD:
            countD = {}
            for trin, OG in trinOG.items():
                countD[OG] = trinCounts.get(trin, 0)

            #just overwrite old orthogroups dictionary
            #orthogroups[sample] = orthogroups['Orthogroup'].map(countD)

            orthogroups[sample] = orthogroups.index.map(countD)

    og_quanttable_file = orthogroups_file.rsplit('.', 1)[0] + '.quant.tsv'
    #orthogroups.to_csv(path_or_buf=og_quanttable_file, sep = '\t', index=True)
    orthogroups.to_csv(path_or_buf=og_quanttable_file, sep = '\t', index=True, na_rep='0')

if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('orthogroups', default = "Orthogroups.tsv")
    p.add_argument('--sample_list')
    p.add_argument('--namemap_directory')
    p.add_argument('--quant_directory')
    p.add_argument('--gtmap_directory', default= "orthogroup_gtmaps")
    p.add_argument('-v','--verbose', action='store_true', default=False)
    args = p.parse_args()
    sys.exit(build_counttable(args.orthogroups, args.gtmap_directory, args.sample_list, args.namemap_directory, args.quant_directory, args.verbose))
