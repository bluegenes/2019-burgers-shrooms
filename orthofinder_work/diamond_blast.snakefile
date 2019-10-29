"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  diamond_blast.snakefile --use-conda # use -n for dry run
"""

import os
import sys
import itertools
from itertools import product


#ortho_dir = "/pylon5/mc5phkp/ntpierce/mmetsp_orthofinder_results/Results_oct_2019_1/WorkingDirectory/"
pep_dir = "/home/ntpierce/2019-burgers-shrooms/mmetsp_info/mmetsp_pep"
ortho_dir = "/home/ntpierce/2019-burgers-shrooms/orthofinder_work/orthfinder_blast_results"
species_ids = "/home/ntpierce/2019-burgers-shrooms/orthofinder_work/SpeciesIDs.txt"
sequence_ids = "/home/ntpierce/2019-burgers-shrooms/orthofinder_work/SequenceIDs.txt"

# build dictionary of speciesID: pepfile
speciesD = {}
with open(species_ids, 'r') as f:
    for line in f:
        num, pepfile = line.split(": ")
        speciesD[num] = pepfile

pepfiles = [os.path.join(ortho_dir, f"Species{x}.fa") for x in speciesD.keys()]

def build_outputs(num_species=658):
    filenames = []
    for i, j in itertools.product(range(num_species), repeat=2):
        if i==j:
            continue
        filenames.append(os.path.join(ortho_dir, f"Blast{i}_{j}.txt.gz"))
    return filenames

num_species =  int(4) #int(658)  # can count this from *fa files in directory
output_blast_filenames = build_outputs(num_species)

rule all:
    input: 
        output_blast_filenames

rule rename_fasta_files_and_contigs:
    input: 
        speciesIDfile = species_ids,
        seqIDfile = sequence_ids,
    output: 
         pepfiles
    params:
        peptide_dir = pep_dir,
        out_dir = ortho_dir,
    conda: 
        "orthofinder_diamond.yml"
    shell:
        """
        python rename_fasta.py --speciesIDs {input.speciesIDfile} --seqIDs {input.seqIDfile} --pep_dir {params.peptide_dir} -o {params.out_dir}
        """

rule diamond_makedb:
    input: 
        #pep = lambda wildcards: os.path.join(ortho_dir, speciesD[wildcards.s1])
        pep = os.path.join(ortho_dir, "Species{s1}.fa"),
    output: 
        os.path.join(ortho_dir, "diamondDBSpecies{s1}.dmnd")
    conda: 
        "orthofinder_diamond.yml"
    params:
        prefix = "diamondDBSpecies"
    shell:
        """
        diamond makedb --in {input} -d {params.prefix}
        """

# for orthofinder, output needs to be of form: Blast{species}_{species}.txt.gz
rule diamond_blastx:
    input: 
        pep = os.path.join(ortho_dir,"Species{s1}.fa"),
        db =  os.path.join(ortho_dir,"diamondDBSpecies{s2}.dmnd")
    output:
        os.path.join(ortho_dir,"Blast{s1}_{s2}.txt.gz")
    conda: 
        "orthofinder_diamond.yml"
    shell:
        """
        diamond blastx -d {input.db} -q {input.pep} -o {output}
        """

