
"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s  biobombe_build_hetnets.snakefile --use-conda
"""

import os
import pandas as pd
from scripts.biobombe_snakemake_utils import read_params
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()


# NOTE - This fill is not added to the repository because it contains
# gene sets with restrictive licenses
#url_prefix = 'https://www.broadinstitute.org/gsea/resources/msigdb/'
#url = '{}{}/msigdb.v{}.{}.gmt'.format(url_prefix, version, version, gene_id)
#! wget --timestamping --no-verbose --directory-prefix 'data' $url

# MSigDB version
msigdb_version = '6.1'
msigdb_gene_id = 'entrez'
        
# Many of the genesets have sub gene sets - process these as well
msigdb_dict = {
    'msigdb': 'full msigdb', # uses same link - can dl full using same rule as subsets
    'c1.all': 'positional gene sets',
    'c2.cgp': 'chemical and genetic perturbations',
    'c2.cp.reactome': 'Reactome gene sets',
    'c3.mir': 'microRNA targets',
    'c3.tft': 'transcription factor targets',
    'c4.cgn': 'cancer gene neighborhoods',
    'c4.cm': 'cancer modules',
    'c5.bp': 'GO biological processes',
    'c5.cc': 'GO cellular components',
    'c5.mf': 'GO molecular functions',
    'c6.all': 'oncogenic signatures',
    'c7.all': 'immunologic signatrues'
}

gene_set_names = msigdb_dict.keys()

rule all:
    input: "msigdb/msigdb.v6.1.entrez.gmt", 
            expand("msigdb/{gene_set}.v{vers}.{gene_id}.gmt", vers = msigdb_version, gene_id=msigdb_gene_id, gene_set=gene_set_names),
            expand("msigdb/{gene_set}_v{vers}.{gene_id}_binary_matrix.tsv.bz2",vers = msigdb_version, gene_id= msigdb_gene_id, gene_set=gene_set_names)


# Download full MSigDB matrix
# NOTE - This fill(full?) is not added to the repository because it contains gene sets with restrictive licenses
#rule httpget_mSigDB:
    #input: HTTP.remote(f"https://www.broadinstitute.org/gsea/resources/msigdb/{vers}/msigdb.v{vers}.{gene_id}.gmt", static=True, keep_local=True, allow_redirects=True)
    #output: "msigdb/msigdb.v{vers}.{gene_id}.gmt" 
    #log: "logs/httpget_msigdb_v{vers}.{gene_id}.log"
    #shell: "mv {input} {output} 2> {log}"

#for gene_set in msigdb_dict:
#    url = '{}{}/{}.v{}.{}.gmt'.format(url_prefix, version, gene_set, version, gene_id)
#    ! wget --timestamping --no-verbose --directory-prefix 'data' $url

# httpget full and geneset-specific msigdb databases
rule httpget_mSigDB:
    input: lambda wildcards: HTTP.remote(f"https://www.broadinstitute.org/gsea/resources/msigdb/{wildcards.vers}/{wildcards.gene_set}.v{wildcards.vers}.{wildcards.gene_id}.gmt", static=True, keep_local=True, allow_redirects=True)
    output: "msigdb/{gene_set}.v{vers}.{gene_id}.gmt"
    log: "logs/httpget_{gene_set}.v{vers}.{gene_id}.log" 
    shell: "mv {input} {output} 2> {log}"


# currently this only does the full one - do we want to do this for all gene sets?
#rule build_mSigDB_binary_matrix:
#    input: "msigdb/msigdb.v{vers}.{gene_id}.gmt"
#    output: "msigdb/fullmsigdb_v{vers}.{gene_id}_binary_matrix.tsv.bz2"
#    log: "logs/make_binary_matrix_msigdb_v{vers}.{gene_id}.log"
#    shell:
#        """
#        python generate_msigdb_matrix.py --input_msigdb_file {input} --out_matrix {output} --checkblacklist
#        """
rule build_mSigDB_genesets_binary_matrix:
    input: "msigdb/{gene_set}.v{vers}.{gene_id}.gmt" 
    output: "msigdb/{gene_set}_v{vers}.{gene_id}_binary_matrix.tsv.bz2"
    log: "logs/make_binary_matrix_msigdb_{gene_set}_v{vers}.{gene_id}.log"
    shell:
        """
        python scripts/generate_msigdb_matrix.py --input_msigdb_file {input} --out_matrix {output} --checkblacklist
        """

#### still need to add::
# integrate-compression-hetnet.ipynb
# permuted-hetnets.ipynb
# process_xCell.ipynb
