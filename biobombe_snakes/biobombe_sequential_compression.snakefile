"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  biobombe_sequential_compression.snakefile --use-conda # use -n for dry run
# note: currently working on ubuntu, NOT working on osx
"""

import os
import sys
import pandas as pd
from biobombe_snakemake_utils import read_params


paramsfile = "config/z_parameter_sweep_MMETSP.tsv"
#paramsfile = "config/z_parameter_sweep_MMETSP_test.tsv"
zsweep_params = read_params(paramsfile)
zsweep_paramsD = zsweep_params.to_dict()
# this dictionary is  {zdim: {variable: value, var2: val2}}


SAMPLE = 'haptophyta_orthogroup'
ZDIMS = zsweep_paramsD.keys()
EXTS = ["tybalt_training_hist.tsv","adage_training_hist.tsv", "sample_corr.tsv.gz", "reconstruction.tsv"] #"gene_corr.tsv.gz"

rule all:
    input: 
         expand("model_results/ensemble_z_results/{zdim}_components/{sample}_{zdim}_components_{ext}", sample=SAMPLE, zdim=ZDIMS, ext=EXTS),
         expand("model_results/ensemble_z_results/{zdim}_components_shuffled/{sample}_{zdim}_components_shuffled_{ext}", sample=SAMPLE, zdim=ZDIMS, ext=EXTS),
         expand("model_results_mad/ensemble_z_results/{zdim}_components/{sample}_{zdim}_components_{ext}", sample=SAMPLE, zdim=ZDIMS, ext=EXTS),
         expand("model_results_mad/ensemble_z_results/{zdim}_components_shuffled/{sample}_{zdim}_components_shuffled_{ext}", sample=SAMPLE, zdim=ZDIMS, ext=EXTS),
         #"model_results/reconstruction_cost_{sample}.tsv", "model_results_mad/reconstruction_cost_{sample}.tsv",


rule preprocess_data:
    input: 
        "data/{sample}.quant.tsv"
    output:
        train = "data/{sample}.train.processed.tsv.gz",
        test = "data/{sample}.test.processed.tsv.gz",
        mad = "data/{sample}.mad.processed.tsv.gz"
    params:
        outdir = "data"
    conda: 
        "environment.yml"
    shell:
        """
        export KERAS_BACKEND=tensorflow
        python process_expression_data.py {input} --mad --output_folder {params.outdir}
        """

rule preprocess_data_scale:
    input: 
        "data/{sample}.quant.tsv"
    output:
        train = "data/{sample}.train.processed.zeroone.tsv.gz",
        test = "data/{sample}.test.processed.zeroone.tsv.gz",
        mad = "data/{sample}.mad.processed.zeroone.tsv.gz"
    params:
        outdir = "data"
    conda: 
        "environment.yml"
    shell:
        """
        export KERAS_BACKEND=tensorflow
        python process_expression_data.py {input} --mad --output_folder {params.outdir} --scale --scale_method "min_max"
        """

rule train_models:
    input: 
        train = "data/{sample}.train.processed.tsv.gz",
        test = "data/{sample}.test.processed.tsv.gz",
        mad_train = "data/{sample}.mad.train90.processed.tsv.gz",
        mad_test = "data/{sample}.mad.test10.processed.tsv.gz"
    output: 
        "model_results/ensemble_z_results/{zdim}_components/{sample}_{zdim}_components_tybalt_training_hist.tsv",
        "model_results/ensemble_z_results/{zdim}_components/{sample}_{zdim}_components_adage_training_hist.tsv",
        #"model_results/ensemble_z_results/{zdim}_components/{sample}_{zdim}_components_gene_corr.tsv.gz",
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
        export KERAS_BACKEND=tensorflow
        python train_models_single_zdim.py  {input.train} {input.test} --basename {wildcards.sample} --paramsfile {params.paramsF} --zdim {wildcards.zdim} --outdir {params.out_dir}
        """

rule train_models_shuffle:
    input:
        train = "data/{sample}.train.processed.tsv.gz",
        test = "data/{sample}.test.processed.tsv.gz",
        mad = "data/{sample}.mad.processed.tsv.gz",
        mad_train = "data/{sample}.mad.train90.processed.tsv.gz",
        mad_test = "data/{sample}.mad.test10.processed.tsv.gz"
    output:
        "model_results/ensemble_z_results/{zdim}_components_shuffled/{sample}_{zdim}_components_shuffled_tybalt_training_hist.tsv",
        "model_results/ensemble_z_results/{zdim}_components_shuffled/{sample}_{zdim}_components_shuffled_adage_training_hist.tsv",
        #"model_results/ensemble_z_results/{zdim}_components_shuffled/{sample}_{zdim}_components_shuffled_gene_corr.tsv.gz",
        "model_results/ensemble_z_results/{zdim}_components_shuffled/{sample}_{zdim}_components_shuffled_sample_corr.tsv.gz",
        "model_results/ensemble_z_results/{zdim}_components_shuffled/{sample}_{zdim}_components_shuffled_reconstruction.tsv",
        directory("model_results/ensemble_z_matrices/{sample}_components_{zdim}")
    params:
        out_dir = "model_results",
        paramsF = paramsfile
    conda:
        'environment.yml'
    shell:
        """
        export KERAS_BACKEND=tensorflow
        python train_models_single_zdim.py   {input.train} {input.test} --basename {wildcards.sample} --paramsfile {params.paramsF} --zdim {wildcards.zdim} --outdir {params.out_dir} --shuffle
        """

rule train_models_mad:
    input:
        train = "data/{sample}.train.processed.tsv.gz",
        test = "data/{sample}.test.processed.tsv.gz",
        mad = "data/{sample}.mad.processed.tsv.gz",
        mad_train = "data/{sample}.mad.train90.processed.tsv.gz",
        mad_test = "data/{sample}.mad.test10.processed.tsv.gz"
    output:
        "model_results_mad/ensemble_z_results/{zdim}_components/{sample}_{zdim}_components_tybalt_training_hist.tsv",
        "model_results_mad/ensemble_z_results/{zdim}_components/{sample}_{zdim}_components_adage_training_hist.tsv",
        #"model_results_mad/ensemble_z_results/{zdim}_components/{sample}_{zdim}_components_gene_corr.tsv.gz",
        "model_results_mad/ensemble_z_results/{zdim}_components/{sample}_{zdim}_components_sample_corr.tsv.gz",
        "model_results_mad/ensemble_z_results/{zdim}_components/{sample}_{zdim}_components_reconstruction.tsv",
        directory("model_results_mad/ensemble_z_matrices/{sample}_components_{zdim}")
    params:
        out_dir = "model_results_mad",
        paramsF = paramsfile
    conda:
        'environment.yml'
    shell:
        """
        export KERAS_BACKEND=tensorflow
        python train_models_single_zdim.py   {input.train} {input.test} --mad_train {input.mad_train} --mad_test {input.mad_test} --basename {wildcards.sample} --paramsfile {params.paramsF} --zdim {wildcards.zdim} --outdir {params.out_dir}
        """

rule train_models_shuffle_mad:
    input:
        train = "data/{sample}.train.processed.tsv.gz",
        test = "data/{sample}.test.processed.tsv.gz",
        mad = "data/{sample}.mad.processed.tsv.gz",
        mad_train = "data/{sample}.mad.train90.processed.tsv.gz",
        mad_test = "data/{sample}.mad.test10.processed.tsv.gz"
    output:
        "model_results_mad/ensemble_z_results/{zdim}_components_shuffled/{sample}_{zdim}_components_shuffled_tybalt_training_hist.tsv",
        "model_results_mad/ensemble_z_results/{zdim}_components_shuffled/{sample}_{zdim}_components_shuffled_adage_training_hist.tsv",
        #"model_results_mad/ensemble_z_results/{zdim}_components_shuffled/{sample}_{zdim}_components_shuffled_gene_corr.tsv.gz",
        "model_results_mad/ensemble_z_results/{zdim}_components_shuffled/{sample}_{zdim}_components_shuffled_sample_corr.tsv.gz",
        "model_results_mad/ensemble_z_results/{zdim}_components_shuffled/{sample}_{zdim}_components_shuffled_reconstruction.tsv",
        directory("model_results_mad/ensemble_z_matrices/{sample}_components_{zdim}")
    params:
        out_dir = "model_results_mad",
        paramsF = paramsfile
    conda:
        'environment.yml'
    shell:
        """
        export KERAS_BACKEND=tensorflow
        python train_models_single_zdim.py   {input.train} {input.test} --mad_train {input.mad_train} --mad_test {input.mad_test} --basename {wildcards.sample} --paramsfile {params.paramsF} --zdim {wildcards.zdim} --outdir {params.out_dir} --shuffle
        """


rule reconstruct_results:
    input:
         expand("model_results/ensemble_z_results/{zdim}_components/{sample}_{zdim}_components_{ext}", sample=SAMPLE, zdim=ZDIMS, ext=EXTS),
         expand("model_results/ensemble_z_results/{zdim}_components_shuffled/{sample}_{zdim}_components_shuffled_{ext}", sample=SAMPLE, zdim=ZDIMS, ext=EXTS),
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

rule reconstruct_results_mad:
    input:
         expand("model_results_mad/ensemble_z_results/{zdim}_components/{sample}_{zdim}_components_{ext}", sample=SAMPLE, zdim=ZDIMS, ext=EXTS),
         expand("model_results_mad/ensemble_z_results/{zdim}_components_shuffled/{sample}_{zdim}_components_shuffled_{ext}", sample=SAMPLE, zdim=ZDIMS, ext=EXTS)
    output:
        recon_cost = "model_results_mad/reconstruction_cost_{sample}.tsv",
        recon_cost_fig = "model_results_mad/figures/reconstruction_cost_{sample}.tsv",
        vae_recon_cost_fig = "model_results_mad/figures/vae_training_reconstruction_{sample}.tsv"
    params:
        results_dir = "model_results_mad",
        figures_dir = "model_results_mad/figures"
        #basename = {wildcards.sample}, # would this work?
    conda:
        "environment.yml"
    shell:
        """
        "visualize_reconstruction.r"
        """
