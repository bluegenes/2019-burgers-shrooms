"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s biobombe.snakefile --use-conda
"""

import os
import pandas as pd


def read_params(pfile):
    if '.tsv' in pfile:
        separator = '\t'
        try:
            paramsDF = pd.read_csv(pfile, dtype=str, sep=separator, index_col=0)
            paramsDF['sweep_values'] = paramsDF['sweep_values'].str.split(',')
        except Exception as e:
            sys.stderr.write(f"\n\tError: {pfile} file is not properly formatted. Please fix.\n\n")
            print(e)
    elif '.xls' in countfile:
        try:
            paramsDF = pd.read_excel(pfile, dtype=str, index_col=0)
        except Exception as e:
            sys.stderr.write(f"\n\tError: {pfile} file is not properly formatted. Please fix.\n\n")
            print(e)
    return paramsDF


adage_params = read_params("config/z_parameter_sweep_adage_MMETSP.tsv")
adage_paramsD = adage_params.to_dict()

tybalt_params = read_params("config/z_parameter_sweep_tybalt_MMETSP.tsv")
tybalt_paramsD = tybalt_params.to_dict()



rule all:
    input: 
        "results/adage/paramsweep_summary.txt", "results/tybalt/paramsweep_summary.txt", "figures/viz_results/z_parameter_final_loss_adage.png"


rule run_adage:
    input: 
        expand("data/{sample}.scaled.tsv", sample = ['haptophyta_orthogroup'])
    output: 
        "results/adage/{sample}_lr{learning_rate}_bs{batch_size}_e{epochs}_sp{sparsity}_ns{noise}_numc{num_components}.tsv",
    params: 
        learning_rate = adage_paramsD['sweep_values']['learning_rate'],
        batch_size = adage_paramsD['sweep_values']['batch_size'],
        epochs = adage_paramsD['sweep_values']['epochs'],
        sparsity = adage_paramsD['sweep_values']['sparsity'],
        noise = adage_paramsD['sweep_values']['noise'],
        num_components = adage_paramsD['sweep_values']['num_components']
    conda: 'environment.yaml'
    shell:
       """
       python adage.py  --input_data {input}
                        --learning_rate {params.learning_rate}
                        --batch_size {params.batch_size}
                        --epochs {params.epochs}
                        --sparsity {params.sparsity}
                        --noise {params.noise}
                        --output_filename {output}
                        --num_components {params.num_components}
                        --subset_mad_genes
        """

rule summaraize_paramsweep_adage:
    input:
        all_adage_runs = expand("results/adage/{sample}_lr{learning_rate}_bs{batch_size}_e{epochs}_sp{sparsity}_ns{noise}_numc{num_components}.tsv", sample= "haptophyta_orthogroup", learning_rate =        adage_paramsD['sweep_values']['learning_rate'], batch_size = adage_paramsD['sweep_values']['batch_size'], epochs = adage_paramsD['sweep_values']['epochs'], sparsity = adage_paramsD['sweep_values']['sparsity'],   noise = adage_paramsD['sweep_values']['noise'], num_components = adage_paramsD['sweep_values']['num_components']),
    output: "results/adage/paramsweep_summary.txt"
    params: 
        results_dir = directory("results/adage")
    shell:
        """ 
        python biobombe_scripts/summarize_paramsweep.py -r {input} -f {output}
        """

rule run_tybalt:
    input: 
        expand("data/{sample}.scaled.tsv", sample = ['haptophyta_orthogroup'])
    output:
        "results/tybalt/{sample}_lr{learning_rate}_bs{batch_size}_e{epochs}_k{kappa}_numc{num_components}.tsv"
    params: 
        learning_rate = tybalt_paramsD['sweep_values']['learning_rate'],
        batch_size = tybalt_paramsD['sweep_values']['batch_size'],
        epochs = tybalt_paramsD['sweep_values']['epochs'],
        kappa = tybalt_paramsD['sweep_values']['kappa'],
        num_components = tybalt_paramsD['sweep_values']['num_components']
    conda: 'environment.yaml'
    shell:
       """
       python vae.py    --input_data {input}
                        --learning_rate {params.learning_rate}
                        --batch_size {params.batch_size}
                        --epochs {params.epochs}
                        --kappa {params.kappa}
                        --output_filename {output}
                        --num_components {params.num_components}
                        --subset_mad_genes
        """

rule summaraize_paramsweep_tybalt:
    input: 
        all_tybalt_runs = expand("results/tybalt/{sample}_lr{learning_rate}_bs{batch_size}_e{epochs}_k{kappa}_numc{num_components}.tsv", sample= "haptophyta_orthogroup", learning_rate = tybalt_paramsD['sweep_values']['learning_rate'], batch_size = tybalt_paramsD['sweep_values']['batch_size'], epochs = tybalt_paramsD['sweep_values']['epochs'], kappa = tybalt_paramsD['sweep_values']['kappa'], num_components = tybalt_paramsD['sweep_values']['num_components']),
    params:
        results_dir = directory("results/tybalt")
    output: "results/tybalt/paramsweep_summary.txt"
    shell:
        """ 
        python biobombe_scripts/summarize_paramsweep.py -r {input} -f {output}
        """


rule visualize_tybalt_paramsweep:
    input: 
        adage = "results/adage/paramsweep_summary.txt", 
        tybalt = "results/tybalt/paramsweep_summary.txt",
    params:
        dataset_name = "haptophyta",
    output: 
        "figures/viz_results/z_parameter_final_loss_adage.png"
    script:
        "visualize-parameter-sweep.r"


