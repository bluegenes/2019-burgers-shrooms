###################################################################
"""
Rename fasta contig names for orthofinder
"""
###################################################################

import sys
import os
import glob
import argparse
import screed


def rename_single_fasta(speciesNum, pepfile, seqmap, pep_out):
    with open(pep_out, 'w') as out:
        with screed.open(pepfile) as seqs:
            for read in seqs:
                try:
                    newname = speciesNum + '_' + seqmap[read.name]
                except:
                    print(read.name)
                    import pdb;pdb.set_trace()
                out.write('>' + newname + '\n')
                out.write(read.sequence + '\n')

def rename_fasta(speciesID_file, seqID_file, pep_dir, out_dir):
    # build dictionary of speciesID: pepfile
    speciesD = {}
    with open(speciesID_file, 'r') as f:
        for line in f:
            num, pepfile = line.strip().split(": ")
            speciesD[num] = pepfile.rsplit(".fasta", 1)[0]
    
    with open(seqID_file, 'r') as f:
        seqD = {}
        species = 0
        prev_species = 0
        for line in f:
            info, longID = line.strip().split(": ", 1)
            species, contigNum = info.split("_")
            if int(species) > prev_species:
                sp = str(prev_species)
                pep_fasta = os.path.join(pep_dir, speciesD[sp])
                pep_fasta_out = os.path.join(out_dir, f"Species{sp}.fa")
                rename_single_fasta(sp, pep_fasta, seqD, pep_fasta_out)
                # clear previous, dictionary
                prev_species = int(species)
                print(species)
                seqD = {}
            seqD[longID] = contigNum
        # do last species
        pep_fasta = os.path.join(pep_dir, speciesD[species])
        pep_fasta_out = os.path.join(out_dir, f"Species{species}.fa")
        rename_single_fasta(pep_fasta, seqD, pep_fasta_out)

# Sequence_IDs.txt file:
# 0_0: Transcript_0|m.1 Transcript_0|g.1  ORF Transcript_0|g.1 Transcript_0|m.1 type:3prime_partial len:404 (-) Transcript_0:1-1209(-)

if __name__ == '__main__':
    """Function: Take in a list of trinity gene names, a dammit namemap,
and a dammit fasta file. Output fasta of matching dammit contigs.
"""
    p = argparse.ArgumentParser()
    p.add_argument('--speciesIDs')
    p.add_argument('--seqIDs')
    p.add_argument('--pep_dir')
    p.add_argument('-o', '--out_dir', default=os.getcwd())
    args = p.parse_args()
    rename_fasta(args.speciesIDs, args.seqIDs, args.pep_dir, args.out_dir)
