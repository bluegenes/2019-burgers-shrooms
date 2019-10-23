"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  diamond_blast.snakefile --use-conda # use -n for dry run
"""

import os
import sys
import itertools
from itertools import product

def build_outputs(num_species=658):
    filenames = []
    for i, j in itertools.product(range(num_species), repeat=2):
        if i==j:
            continue
        filenames.append(f"Blast{i}_{j}.txt.gz")
    return filenames

num_species = int(658)  # can count this from *fa files in directory
output_blast_filenames = build_outputs(num_species)

rule all:
    input: 
        output_blast_filenames

rule diamond_makedb:
    input: 
        pep = "{s1}.fa",
    output: 
        "diamondDBSpecies{s1}.dmnd"
    conda: 
        "diamond_environment.yml"
    params:
        prefix = "diamondDBSpecies397"
    shell:
        """
        diamond makedb --in {input} -d {output}
        """

# for orthofinder, output needs to be of form: Blast{species}_{species}.txt.gz
rule diamond_blastx:
    input: 
        pep = "{s1}.fa",
        db = "diamondDBSpecies{s2}.dmnd"
    output:
        "Blast{s1}_{s2}.txt.gz"
    conda: 
        "diamond_environment.yml"
    shell:
        """
        diamond blastx -d {input.db} -q {input.pep} -o {output}
        """

