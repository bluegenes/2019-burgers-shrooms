"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  biobombe_initial_ksweep.snakefile --use-conda -n
"""

import os
import pandas as pd
from scripts.biobombe_snakemake_utils import read_params
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()

adage_params = read_params("config/initial_z_parameter_sweep_adage_MMETSP.tsv")
#adage_params = read_params("config/initial_z_parameter_sweep_adage_MMETSP_test.tsv")
adage_params['sweep_values'] = adage_params['sweep_values'].str.split(',')
adage_paramsD = adage_params.to_dict()

tybalt_params = read_params("config/initial_z_parameter_sweep_tybalt_MMETSP.tsv")
#tybalt_params = read_params("config/initial_z_parameter_sweep_tybalt_MMETSP_test.tsv")
tybalt_params['sweep_values'] = tybalt_params['sweep_values'].str.split(',')
tybalt_paramsD = tybalt_params.to_dict()


rule all:
    input: 
        "figures/viz_results/z_parameter_final_loss_adage.png", "figures/viz_results/z_parameter_final_loss_tybalt.png",
        expand("results/tybalt/{sample}_lr{learning_rate}_bs{batch_size}_e{epochs}_k{kappa}_numc{num_components}.tsv", sample= "haptophyta_orthogroup", learning_rate =                   tybalt_paramsD['sweep_values']['learning_rate'], batch_size = tybalt_paramsD['sweep_values']['batch_size'], epochs = tybalt_paramsD['sweep_values']['epochs'], kappa =                   tybalt_paramsD['sweep_values']['kappa'], num_components = tybalt_paramsD['sweep_values']['num_components']), 
        expand("results/adage/{sample}_lr{learning_rate}_bs{batch_size}_e{epochs}_sp{sparsity}_ns{noise}_numc{num_components}.tsv", sample= "haptophyta_orthogroup", learning_rate         =        adage_paramsD['sweep_values']['learning_rate'], batch_size = adage_paramsD['sweep_values']['batch_size'], epochs = adage_paramsD['sweep_values']['epochs'], sparsity =                             adage_paramsD['sweep_values']['sparsity'],   noise = adage_paramsD['sweep_values']['noise'], num_components = adage_paramsD['sweep_values']['num_components']),
        
        # use this for specific param compbas version
        #"figures/viz_results/{sample}_z_parameter_final_loss_dae.png", "figures/viz_results/{sample}_z_parameter_final_loss_vae.png"

rule download_data:
    input: HTTP.remote("https://osf.io/ek9nu/download")
    output: "data/haptophyta_orthogroup.quant.tsv" 
    log: "logs/download_haptophyta_orthogroup.log"
    shell: "mv {input} {output} 2> {log}"

rule preprocess_data:
    input:
        "data/{sample}.quant.tsv"
    output:
        processed = "data/{sample}.processed.tsv.gz",
        train = "data/{sample}.processed.train.tsv.gz",
        test = "data/{sample}.processed.test.tsv.gz",
        mad = "data/{sample}.processed.mad.tsv.gz"
    params:
        outdir = "data"
    conda:
        "environment.yml"
    shell:
        """
        python scripts/process_expression_data.py {input} --mad --output_folder {params.outdir}
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
        python scripts/process_expression_data.py {input} --mad --output_folder {params.outdir} --scale --scale_method "min_max"
        """


rule run_adage:
    input: 
        expand("data/{sample}.processed.tsv.gz", sample = ['haptophyta_orthogroup'])
    output: 
        "results/adage/{sample}_lr{learning_rate}_bs{batch_size}_e{epochs}_sp{sparsity}_ns{noise}_numc{num_components}.tsv",
    conda: 'environment.yml'
    shell:
       """
        export KERAS_BACKEND=tensorflow
       python scripts/adage.py  --input_data {input} --learning_rate {wildcards.learning_rate} --batch_size {wildcards.batch_size} --epochs {wildcards.epochs} --sparsity {wildcards.sparsity} --noise {wildcards.noise} --output_filename {output} --num_components {wildcards.num_components} --subset_mad_genes 8000 --scale
        """

rule summarize_paramsweep_adage:
    input:
        all_adage_runs = expand("results/adage/{sample}_lr{learning_rate}_bs{batch_size}_e{epochs}_sp{sparsity}_ns{noise}_numc{num_components}.tsv", sample= "haptophyta_orthogroup", learning_rate =        adage_paramsD['sweep_values']['learning_rate'], batch_size = adage_paramsD['sweep_values']['batch_size'], epochs = adage_paramsD['sweep_values']['epochs'], sparsity = adage_paramsD['sweep_values']['sparsity'],   noise = adage_paramsD['sweep_values']['noise'], num_components = adage_paramsD['sweep_values']['num_components']),
    output: "results/adage/paramsweep_summary.txt"
    params: 
        results_dir = directory("results/adage")
    shell:
        """ 
        python scripts/summarize_paramsweep.py -r {params.results_dir} -f {output}
        """

rule run_tybalt:
    input: 
        expand("data/{sample}.processed.tsv.gz", sample = ['haptophyta_orthogroup'])
    output:
        "results/tybalt/{sample}_lr{learning_rate}_bs{batch_size}_e{epochs}_k{kappa}_numc{num_components}.tsv"
    conda: 'environment.yml'
    shell:
       """
        export KERAS_BACKEND=tensorflow
       python scripts/vae.py    --input_data {input} --learning_rate {wildcards.learning_rate} --batch_size {wildcards.batch_size} --epochs {wildcards.epochs} --kappa {wildcards.kappa} --output_filename {output} --num_components {wildcards.num_components} --subset_mad_genes 8000 --scale
        """

rule summarize_paramsweep_tybalt:
    input: 
        all_tybalt_runs = expand("results/tybalt/{sample}_lr{learning_rate}_bs{batch_size}_e{epochs}_k{kappa}_numc{num_components}.tsv", sample= "haptophyta_orthogroup", learning_rate = tybalt_paramsD['sweep_values']['learning_rate'], batch_size = tybalt_paramsD['sweep_values']['batch_size'], epochs = tybalt_paramsD['sweep_values']['epochs'], kappa = tybalt_paramsD['sweep_values']['kappa'], num_components = tybalt_paramsD['sweep_values']['num_components']),
    params:
        results_dir = directory("results/tybalt")
    output: "results/tybalt/paramsweep_summary.txt"
    shell:
        """ 
        python scripts/summarize_paramsweep.py -r {params.results_dir} -f {output}
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
    conda:
        "environment.yml"
    script:
        "scripts/visualize-parameter-sweep.R"


### code below can be used to run individual dae/vae with specified parameter combos


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


def generate_vae_combos(w):
    vae_combos = []
    for component, zparams in zsweep_paramsD.items():
        learning_rate = zparams["vae_lr"]
        batch_size = zparams["vae_batch_size"]
        epochs = zparams["vae_epochs"]
        kappa = zparams["vae_kappa"]
        vae_combos.append(f"results/vae/{w.sample}_comp{component}_lr{learning_rate}_bs{batch_size}_e{epochs}_k{kappa}.tsv")
    return vae_combos

rule run_dae:
    input: 
        "data/{sample}.scaled.tsv"
    output: 
        "results/dae/{sample}_comp{component}_lr{learning_rate}_bs{batch_size}_e{epochs}_sp{sparsity}_ns{noise}.tsv"
    conda: 'environment.yml'
    shell:
       """
       python scripts/dae.py    --input_data {input}
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
        python scripts/summarize_paramsweep.py -r {params.results_dir} -f {output}
        """

rule run_vae:
    input:
        "data/{sample}.scaled.tsv"
    output:
        "results/vae/{sample}_comp{component}_lr{learning_rate}_bs{batch_size}_e{epochs}_k{kappa}.tsv"
    conda: 'environment.yml'
    shell:
       """
       python scripts/vae.py    --input_data {input}
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
        python scripts/summarize_paramsweep.py -r {params.results_dir} -f {output}
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
        "scripts/visualize-parameter-sweep.r"

