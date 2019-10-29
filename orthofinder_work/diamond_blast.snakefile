"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  diamond_blast.snakefile --use-conda # use -n for dry run
"""

import os
import sys
import itertools
from itertools import product


#ortho_dir = "/pylon5/mc5phkp/ntpierce/mmetsp_orthofinder_results/Results_oct_2019_1/WorkingDirectory/"
ortho_dir =  "./" #"/pylon5/mc5phkp/ntpierce/mmetsp_pep"
pep_dir =  "../mmetsp_info/mmetsp_pep" #"/pylon5/mc5phkp/ntpierce/mmetsp_pep"
species_ids = "SpeciesIDs.txt"

# build dictionary of speciesID: pepfile
speciesD = {}
with open(species_ids, 'r') as f:
    for line in f:
        num, pepfile = line.split(": ")
        speciesD[num] = pepfile

def build_outputs(num_species=658):
    filenames = []
    for i, j in itertools.product(range(num_species), repeat=2):
        if i==j:
            continue
        filenames.append(os.path.join(ortho_dir, f"Blast{i}_{j}.txt.gz"))
    return filenames

num_species = int(658)  # can count this from *fa files in directory
output_blast_filenames = build_outputs(num_species)

rule all:
    input: 
        output_blast_filenames


rule diamond_makedb:
    input: 
        pep = os.path.join(pep_dir, speciesD[{s1}])
        #pep = "Species{s1}.fa",
    output: 
        "diamondDBSpecies{s1}.dmnd"
    conda: 
        "diamond_environment.yml"
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
        db = os.path.join(ortho_dir,"diamondDBSpecies{s2}.dmnd")
    output:
        os.path.join(ortho_dir,"Blast{s1}_{s2}.txt.gz")
    conda: 
        "diamond_environment.yml"
    shell:
        """
        diamond blastx -d {input.db} -q {input.pep} -o {output}
        """

