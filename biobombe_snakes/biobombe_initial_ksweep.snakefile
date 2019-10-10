"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  biobombe_initial_ksweep.snakefile --use-conda -n
"""

import os
import pandas as pd
from biobombe_snakemake_utils import read_params


adage_params = read_params("config/initial_z_parameter_sweep_adage_MMETSP.tsv")
adage_params['sweep_values'] = adage_params['sweep_values'].str.split(',')
adage_paramsD = adage_params.to_dict()

tybalt_params = read_params("config/initial_z_parameter_sweep_tybalt_MMETSP.tsv")
tybalt_params['sweep_values'] = tybalt_params['sweep_values'].str.split(',')
tybalt_paramsD = tybalt_params.to_dict()



rule all:
    input: 
        "figures/viz_results/z_parameter_final_loss_adage.png", "figures/viz_results/z_parameter_final_loss_tybalt.png"


rule run_adage:
    input: 
        expand("data/{sample}.scaled.tsv", sample = ['haptophyta_orthogroup'])
    output: 
        "results/adage/{sample}_lr{learning_rate}_bs{batch_size}_e{epochs}_sp{sparsity}_ns{noise}_numc{num_components}.tsv",
    conda: 'environment.yml'
    shell:
       """
       python adage.py  --input_data {input}
                        --learning_rate {wildcards.learning_rate}
                        --batch_size {wildcards.batch_size}
                        --epochs {wildcards.epochs}
                        --sparsity {wildcards.sparsity}
                        --noise {wildcards.noise}
                        --output_filename {output}
                        --num_components {wildcards.num_components}
                        --subset_mad_genes
        """

rule summarize_paramsweep_adage:
    input:
        all_adage_runs = expand("results/adage/{sample}_lr{learning_rate}_bs{batch_size}_e{epochs}_sp{sparsity}_ns{noise}_numc{num_components}.tsv", sample= "haptophyta_orthogroup", learning_rate =        adage_paramsD['sweep_values']['learning_rate'], batch_size = adage_paramsD['sweep_values']['batch_size'], epochs = adage_paramsD['sweep_values']['epochs'], sparsity = adage_paramsD['sweep_values']['sparsity'],   noise = adage_paramsD['sweep_values']['noise'], num_components = adage_paramsD['sweep_values']['num_components']),
    output: "results/adage/paramsweep_summary.txt"
    params: 
        results_dir = directory("results/adage")
    shell:
        """ 
        python biobombe_scripts/summarize_paramsweep.py -r {params.results_dir} -f {output}
        """

rule run_tybalt:
    input: 
        expand("data/{sample}.scaled.tsv", sample = ['haptophyta_orthogroup'])
    output:
        "results/tybalt/{sample}_lr{learning_rate}_bs{batch_size}_e{epochs}_k{kappa}_numc{num_components}.tsv"
    conda: 'environment.yml'
    shell:
       """
       python vae.py    --input_data {input}
                        --learning_rate {wildcards.learning_rate}
                        --batch_size {wildcards.batch_size}
                        --epochs {wildcards.epochs}
                        --kappa {wildcards.kappa}
                        --output_filename {output}
                        --num_components {wildcards.num_components}
                        --subset_mad_genes
        """

rule summarize_paramsweep_tybalt:
    input: 
        all_tybalt_runs = expand("results/tybalt/{sample}_lr{learning_rate}_bs{batch_size}_e{epochs}_k{kappa}_numc{num_components}.tsv", sample= "haptophyta_orthogroup", learning_rate = tybalt_paramsD['sweep_values']['learning_rate'], batch_size = tybalt_paramsD['sweep_values']['batch_size'], epochs = tybalt_paramsD['sweep_values']['epochs'], kappa = tybalt_paramsD['sweep_values']['kappa'], num_components = tybalt_paramsD['sweep_values']['num_components']),
    params:
        results_dir = directory("results/tybalt")
    output: "results/tybalt/paramsweep_summary.txt"
    shell:
        """ 
        python biobombe_scripts/summarize_paramsweep.py -r {params.results_dir} -f {output}
        """


rule visualize_paramsweep:
    input: 
        adage = "results/adage/paramsweep_summary.txt", 
        tybalt = "results/tybalt/paramsweep_summary.txt",
    params:
        dataset_name = "haptophyta",
    output: 
        "figures/viz_results/z_parameter_final_loss_adage.png",
        "figures/viz_results/z_parameter_final_loss_tybalt.png"
    script:
        "visualize-parameter-sweep.r"


