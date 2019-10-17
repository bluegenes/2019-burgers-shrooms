"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  biobombe_sequential_compression.snakefile --use-conda -n
"""

import os
import sys
import pandas as pd
from biobombe_snakemake_utils import read_params


paramsfile = "config/z_parameter_sweep_MMETSP.tsv"
zsweep_params = read_params(paramsfile)
zsweep_paramsD = zsweep_params.to_dict()
# this dictionary is  {zdim: {variable: value, var2: val2}}


SAMPLE = 'haptophyta_orthogroup'
ZDIMS = zsweep_paramsD.keys()
EXTS = ["tybalt_training_hist.tsv","adage_training_hist.tsv", "gene_corr.tsv.gz", "sample_corr.tsv.gz", "reconstruction.tsv"]

rule all:
    input: 
         expand("model_results/ensemble_z_results/{zdim}_components/{sample}_{zdim}_components_{ext}", sample=SAMPLE, zdim=ZDIMS, ext=EXTS),
         expand("model_results/ensemble_z_results/{zdim}_components_shuffled/{sample}_{zdim}_shuffled_components_{ext}", sample=SAMPLE, zdim=ZDIMS, ext=EXTS)

rule preprocess_data:
    input: 
        "data/{sample}.quant.tsv"
    output:
        train = "data/{sample}.processed.train.tsv.gz",
        test = "data/{sample}.processed.test.tsv.gz",
        mad = "data/{sample}.processed.mad.tsv.gz"
    params:
        outdir = "data"
    conda: 
        "environment.yml"
    shell:
        """
        python process_expression_data.py {input} --mad --output_folder {params.outdir}
        """

rule train_models:
    input: 
        train = "data/{sample}.processed.train.tsv.gz",
        test = "data/{sample}.processed.test.tsv.gz",
        mad = "data/{sample}.processed.mad.tsv.gz"
    output: 
        "model_results/ensemble_z_results/{zdim}_components/{sample}_{zdim}_components_tybalt_training_hist.tsv",
        "model_results/ensemble_z_results/{zdim}_components/{sample}_{zdim}_components_adage_training_hist.tsv",
        "model_results/ensemble_z_results/{zdim}_components/{sample}_{zdim}_components_gene_corr.tsv.gz",
        "model_results/ensemble_z_results/{zdim}_components/{sample}_{zdim}_components_sample_corr.tsv.gz",
        "model_results/ensemble_z_results/{zdim}_components/{sample}_{zdim}_components_reconstruction.tsv",
        directory("model_results/ensemble_z_matrices/{sample}_components_{zdim}")
    params:
        out_dir = "model_results",
        paramsF = paramsfile
    conda: 
        'environment.yml'
    shell:
        """
        python train_models_single_zdim.py   {input.train} {input.test}
                                             --input_mad {input.mad}
                                             --basename {wildcards.sample}
                                             --paramsfile {params.paramsF}
                                             --zdim {wildcards.zdim}
                                             --outdir {params.out_dir}
        """

rule train_models_shuffle:
    input:
        train = "data/{sample}.processed.train.tsv.gz",
        test = "data/{sample}.processed.test.tsv.gz",
        mad = "data/{sample}.processed.mad.tsv.gz"
    output:
        "model_results/ensemble_z_results/{zdim}_components_shuffled/{sample}_{zdim}_shuffled_components_tybalt_training_hist.tsv",
        "model_results/ensemble_z_results/{zdim}_components_shuffled/{sample}_{zdim}_shuffled_components_adage_training_hist.tsv",
        "model_results/ensemble_z_results/{zdim}_components_shuffled/{sample}_{zdim}_shuffled_components_gene_corr.tsv.gz",
        "model_results/ensemble_z_results/{zdim}_components_shuffled/{sample}_{zdim}_shuffled_components_sample_corr.tsv.gz",
        "model_results/ensemble_z_results/{zdim}_components_shuffled/{sample}_{zdim}_shuffled_components_reconstruction.tsv",
        directory("model_results/ensemble_z_matrices/{sample}_components_{zdim}")
    params:
        out_dir = "model_results",
        paramsF = paramsfile
    conda:
        'environment.yml'
    shell:
        """
        python train_models_single_zdim.py  {input.train} {input.test}
                                             --input_mad {input.mad}
                                             --basename {wildcards.sample}
                                             --paramsfile {params.paramsF}
                                             --zdim {wildcards.zdim}
                                             --outdir {params.out_dir}
                                             --shuffle
        """

rule reconstruct_results:
    input:
    output:
        recon_cost = "model_results/reconstruction_cost_{sample}.tsv",
        recon_cost_fig = "model_results/figures/reconstruction_cost_{sample}.tsv",
        vae_recon_cost_fig = "model_results/figures/vae_training_reconstruction_{sample}.tsv"
    params:
        results_dir = "model_results",
        figures_dir = "model_results/figures"
        #basename = {wildcards.sample}, # would this work?
    conda:
        "environment.yml"
    shell:
        """
        "visualize_reconstruction.r"
        """

