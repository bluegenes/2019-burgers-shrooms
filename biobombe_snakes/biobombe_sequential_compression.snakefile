"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  biobombe_sequential_compression.snakefile --use-conda -n
"""

import os
import sys
import pandas as pd
from biobombe_snakemake_utils import read_params

zsweep_params = read_params("config/z_parameter_sweep_MMETSP.tsv")
zsweep_paramsD = zsweep_params.to_dict()

# this dictionary is  {component: {variable: value, var2: val2}}


SAMPLE = 'haptophyta_orthogroup'

rule all:
    input: 
#        f"results/dae/{SAMPLE}_paramsweep_summary.txt",f"results/vae/{SAMPLE}_paramsweep_summary.txt"
        f"figures/viz_results/{SAMPLE}_z_parameter_final_loss_dae.png", f"figures/viz_results/{SAMPLE}_z_parameter_final_loss_vae.png"


def generate_dae_combos(w):
    dae_combos = []
    for component, zparams in zsweep_paramsD.items():
        learning_rate = zparams["dae_lr"]
        batch_size = zparams["dae_batch_size"]
        epochs = zparams["dae_epochs"]
        sparsity = zparams["dae_sparsity"]
        noise = zparams["dae_noise"]
        dae_combos.append(f"results/dae/{w.sample}_comp{component}_lr{learning_rate}_bs{batch_size}_e{epochs}_sp{sparsity}_ns{noise}.tsv")
    return dae_combos

rule run_dae:
    input: 
        "data/{sample}.scaled.tsv"
    output: 
        "results/dae/{sample}_comp{component}_lr{learning_rate}_bs{batch_size}_e{epochs}_sp{sparsity}_ns{noise}.tsv"
    conda: 'environment.yml'
    shell:
       """
       python dae.py    --input_data {input}
                        --learning_rate {wildcards.learning_rate}
                        --batch_size {wildcards.batch_size}
                        --epochs {wildcards.epochs}
                        --sparsity {wildcards.sparsity}
                        --noise {wildcards.noise}
                        --output_filename {output}
                        --subset_mad_genes
        """

rule summarize_paramsweep_dae:
    input: generate_dae_combos
    output: "results/dae/{sample}_paramsweep_summary.txt"
    params: 
        results_dir = directory("results/dae")
    shell:
        """ 
        python biobombe_scripts/summarize_paramsweep.py -r {params.results_dir} -f {output}
        """

def generate_vae_combos(w):
    vae_combos = []
    for component, zparams in zsweep_paramsD.items():
        learning_rate = zparams["vae_lr"]
        batch_size = zparams["vae_batch_size"]
        epochs = zparams["vae_epochs"]
        kappa = zparams["vae_kappa"]
        vae_combos.append(f"results/vae/{w.sample}_comp{component}_lr{learning_rate}_bs{batch_size}_e{epochs}_k{kappa}.tsv")
    return vae_combos

rule run_vae:
    input: 
        "data/{sample}.scaled.tsv"
    output:
        "results/vae/{sample}_comp{component}_lr{learning_rate}_bs{batch_size}_e{epochs}_k{kappa}.tsv"
    conda: 'environment.yml'
    shell:
       """
       python vae.py    --input_data {input}
                        --learning_rate {wildcards.learning_rate}
                        --batch_size {wildcards.batch_size}
                        --epochs {wildcards.epochs}
                        --kappa {wildcards.kappa}
                        --output_filename {output}
                        --subset_mad_genes
        """

rule summarize_paramsweep_vae:
    input: generate_vae_combos
    params:
        results_dir = directory("results/vae")
    output: 
        "results/vae/{sample}_paramsweep_summary.txt"
    shell:
        """ 
        python biobombe_scripts/summarize_paramsweep.py -r {params.results_dir} -f {output}
        """


rule visualize_tybalt_paramsweep:
    input: 
        dae = "results/dae/{sample}_paramsweep_summary.txt", 
        vae = "results/vae/{sample}_paramsweep_summary.txt",
    params:
        dataset_name = "haptophyta_orthogroups",
    output: 
        "figures/viz_results/{sample}_z_parameter_final_loss_dae.png",
        "figures/viz_results/{sample}_z_parameter_final_loss_vae.png"
    script:
        "visualize-parameter-sweep.r"


