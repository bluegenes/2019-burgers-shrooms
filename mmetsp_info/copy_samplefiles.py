import os
import sys
import argparse
import glob
from shutil import copyfile


# requires python >= 3.6
# run: python copy_samplefiles.py haptophyta.txt --source ../../pep/ --destination ../../haptophyta_pep


def copy_assemblyfiles(samplelist, source, destination):
    if not os.path.exists(destination):
        try: 
            os.mkdir(destination)
        except Exception as e:
            sys.stderr.write(f"\n\tError: cannot make {destination} destination folder. Please fix.\n\n")
            sys.exit(-1)
            
    with open(samplelist, 'r') as f:
        for line in f:
            sample = line.strip()
            samplefile = glob.glob(os.path.join(source, sample + '*'))
            if len(samplefile) > 1:
                sys.stderr.write(f"\n\tError: only expecting a single match per sample name. Fix samplenames (this sample: {sample}) or edit this script :) \n\n")
                sys.exit(-1)
            sampleF = samplefile[0]
            dest = os.path.join(destination, os.path.basename(sampleF))
            if dest.endswith('.pep'):
                #dest = dest.split('.pep')[0] + '.fasta'
                dest = dest + '.fasta'
            copyfile(sampleF, dest)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('sample_list')
    p.add_argument('--source', default = os.getcwd())
    p.add_argument('--destination', required=True)
    args = p.parse_args()
    sys.exit(copy_assemblyfiles(args.sample_list, args.source, args.destination))

