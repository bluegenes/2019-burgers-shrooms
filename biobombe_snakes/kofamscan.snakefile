"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  diamond_blast.snakefile --use-conda # use -n for dry run
Function: Download and run kofamscan 
"""


import os
import sys
import itertools

from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
FTP = FTPRemoteProvider()
HTTP = HTTPRemoteProvider()

#BASE_DIR = "/Users/tessa/Dropbox/dib-lab/2019-burgers-shrooms/"
BASE_DIR = "/pylon5/mc5phkp/ntpierce/"

MMETSP_INFO_DIR = os.path.join(BASE_DIR, "2019-burgers-shrooms", "mmetsp_info")
#ORTHO_DIR = os.path.join(BASE_DIR, "Orthogroups_Haptophyta/Results_june_2019/Orthogroup_Sequences")
ORTHO_DIR = os.path.join(BASE_DIR, "haptophyta_orthofinder_results/Results_june_2019/Orthogroup_Sequences")
KOFAM_DIR = os.path.join(BASE_DIR,"haptophyta_kofam")


#ortholog_quantfile = os.path.join(MMETSP_INFO_DIR, "haptophyta_orthogroup.quant.tsv")
ortholog_quantfile = os.path.join(MMETSP_INFO_DIR, "haptophyta_orthogroup.quant_test.tsv")

with open(ortholog_quantfile, 'r') as f:
    ORTHOLOG_NAMES = [x.strip().split('\t')[0] for x in f.readlines()[1:]]


rule all:
    input: 
        expand(os.path.join(KOFAM_DIR, "kofamscan_results", "{ortholog}.txt"), ortholog = ORTHOLOG_NAMES)

##rule download_quantfile:
#    input: HTTP.remote("https://osf.io/ek9nu/download")
#    output: os.path.join(KOFAM_DIR, "data/haptophyta_orthogroup.quant.tsv")
#    log: "logs/download_haptophyta_orthogroup.log"
#    shell: "mv {input} {output} 2> {log}"


# download ko list
rule download_ko_list:
    input: FTP.remote("ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz", static=True, keep_local=True, immediate_close=True)
    output: 
        gz_file = os.path.join(KOFAM_DIR,"ko_list.gz"),
        unzipped = os.path.join(KOFAM_DIR, "ko_list"),
    log: os.path.join(KOFAM_DIR,"logs/download_ko_list.log")
    shell: 
        """
        mv {input} {output.gz_file}
        gunzip -c {output.gz_file} > {output.unzipped} 
        """

# download hmm profiles
rule download_ko_profiles:
    input: FTP.remote("ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz", static=True, keep_local=True, immediate_close=True)
    output: 
        gz_file = os.path.join(KOFAM_DIR, "profiles.tar.gz"),
        unzipped = directory(os.path.join(KOFAM_DIR,"profiles"))
    log: os.path.join(KOFAM_DIR,"logs/download_profiles.log")
    params:
        kofam_dir = KOFAM_DIR,
    shell: 
        """
        mv {input} {output.gz_file} 2> {log}
        tar xzf {output.gz_file} --directory {params.kofam_dir}
        """

rule download_kofam_readme:
    input: FTP.remote("ftp://ftp.genome.jp/pub/tools/kofamscan/README.md", static=True, keep_local=True, immediate_close=True)
    output: os.path.join(KOFAM_DIR,"README.md")
    log: os.path.join(KOFAM_DIR, "logs/download_kofam_readme.log")
    shell: "mv {input} {output} 2> {log}"

# download program (replace with conda install in the future!?)
rule download_kofamscan_program:
    input: FTP.remote("ftp://ftp.genome.jp/pub/tools/kofamscan/kofamscan-1.1.0.tar.gz", static=True, keep_local=True, immediate_close=True)
    output: 
        unzipped = directory(os.path.join(KOFAM_DIR,"kofamscan-1.1.0")),
        gz_file = os.path.join(KOFAM_DIR, "kofamscan-1.1.0.tar.gz"),
        executable = os.path.join(KOFAM_DIR,"kofamscan-1.1.0/exec_annotation")
    log: os.path.join(KOFAM_DIR,"logs/download_kofamscan_program.log")
    params:
        kofam_dir = KOFAM_DIR,
    shell: 
        """
        mv {input} {output.gz_file} 
        tar xzf {output.gz_file} --directory {params.kofam_dir}
        """

# config.yml needs to be in same dir as the kofamscan executable
# https://github.com/takaram/kofam_scan
rule print_config:
    input: 
        profile_path = os.path.join(KOFAM_DIR,"profiles"),
        ko_list = os.path.join(KOFAM_DIR, "ko_list") 
    output: os.path.join(KOFAM_DIR,"kofamscan-1.1.0", "config.yml")
    log: os.path.join(KOFAM_DIR,"logs/print_config.txt")
    threads: 8
    shell:
        """
        echo "# Path to your KO-HMM database" >> {output}
        echo "profile: {input.profile_path}" >> {output}
        echo "# Path to the KO list file" >> {output}
        echo "ko_list: {input.ko_list}" >> {output}
        echo "# Number of hmmsearch processes to be run parallelly" >> {output}
        echo "cpu: {threads}" >> {output}
        """

rule run_kofamscan:
    input:
        fasta = os.path.join(ORTHO_DIR, "{ortholog}.fa"),
        executable = os.path.join(KOFAM_DIR, "kofamscan-1.1.0/exec_annotation"),
        config = os.path.join(KOFAM_DIR, "kofamscan-1.1.0", "config.yml"), # don't actually need this anymore, bc just calling params in executable command
        profile_dir = os.path.join(KOFAM_DIR,"profiles"),
        ko_list = os.path.join(KOFAM_DIR, "ko_list") 
    output:
        #os.path.join(ORTHO_DIR, "kofamscan_results/{ortholog}.txt")
        os.path.join(KOFAM_DIR, "kofamscan_results/{ortholog}.txt")
    log:
        os.path.join(KOFAM_DIR, "logs/kofamscan_results_{ortholog}.log")
    threads: 28
    conda:
        "kofamscan-env.yml"
    shell:
        """
        {input.executable} --ko-list {input.ko_list} --profile {input.profile_dir} --cpu {threads} -f mapper -o {output} {input.fasta}
        """

