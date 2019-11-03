"""
2018 Gregory Way, modified 2019 N. Tessa Pierce
#2.ensemble-z-analysis/scripts/train_models_given_z.py

# train_models_single_zdim.py

This script will train various compression models given a specific z dimension.
Each model will train several times with different initializations.

The script pulls hyperparameters from a parameter file that was determined
after initial hyperparameter sweeps testing latent dimensionality.

Usage:

    python train_models.py

    With required command line arguments:

        --input_train       Input train dataset
        --input_test        Input test dataset
        --basename          Basename for file naming
        --zdim    The z dimensionality we're testing
        --paramsfile      A tsv file (param by z dimension) indicating the
                            specific parameter combination for the z dimension
        --out_dir           The directory to store the output files

    And optional command line arguments

        --num_seeds         The number of specific models to generate
                                default: 5
        --shuffle           If provided, shuffle the input expression matrix

Output:

The script will save associated weights and z matrices for each permutation as
well as reconstruction costs for all algorithms, sample specific correlations,
and training histories for Tybalt and ADAGE models. The population of models
(ensemble) are saved, across dimensions z, for downstream evaluation.
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from tybalt.data_models import DataModel
from biobombe_snakemake_utils import *


def get_recon_correlation(df, recon_mat_dict, algorithm, cor_type,
                          genes=False):
    """
    Get gene or sample correlations between input and reconstructed input

    Arguments:
    df - the input dataframe
    recon_mat_dict - dictionary of different algorithms reconstructions
    algorithm - string representing the compression algorithm
    cor_type - string representing Pearson or Spearman correlations
    genes - boolean if to calculate correaltion over genes (sample by default)
    """
    recon_mat = recon_mat_dict[algorithm]
    if genes:
        df = df.T
        recon_mat = recon_mat.T
    if cor_type == 'pearson':
        r = [pearsonr(recon_mat.iloc[x, :],
                      df.iloc[x, :])[0] for x in range(df.shape[0])]
    elif cor_type == 'spearman':
        r = [spearmanr(recon_mat.iloc[x, :],
                       df.iloc[x, :])[0] for x in range(df.shape[0])]
    return r


def compile_corr_df(pearson_list, spearman_list, algorithm_list, column_names,
                    seed, data_type):
    """
    Compile together correlations across algorithms

    Arguments:
    pearson_list - a list of pearson correlations across algorithms
    spearman_list - a list of spearman correlations across algorithms
    algorithm_list - list of algorithm names
    column_names - list of names supplied to the compiled dataframe
    seed - the current random seed
    data_type - training or testing set
    """
    pearson_df = pd.DataFrame(pearson_list,
                              index=algorithm_list,
                              columns=column_names)
    pearson_df.index.name = 'algorithm'
    spearman_df = pd.DataFrame(spearman_list,
                               index=algorithm_list,
                               columns=column_names)
    spearman_df.index.name = 'algorithm'
    pearson_df = pearson_df.reset_index().melt(id_vars=['algorithm'],
                                               var_name='id',
                                               value_name='correlation')
    spearman_df = spearman_df.reset_index().melt(id_vars=['algorithm'],
                                                 var_name='id',
                                                 value_name='correlation')
    corr_df = pd.concat([pearson_df.assign(cor_type='pearson'),
                         spearman_df.assign(cor_type='spearman')])
    corr_df = corr_df.assign(seed=seed)
    corr_df = corr_df.assign(data=data_type)
    return corr_df


def train_models(basename, input_train, input_test, zdim, paramsD, out_dir, num_seeds, shuffle, madfile, num_mad_genes, algorithms = ['pca', 'ica', 'nmf', 'dae', 'vae']):

    # Set output directory and file names
    train_dir = os.path.join(out_dir, 'ensemble_z_results', f'{zdim}_components')
    if shuffle:
        train_dir = f'{train_dir}_shuffled'

    if not os.path.exists(train_dir):
        os.makedirs(train_dir)

    if shuffle:
        file_pre = f'{basename}_{zdim}_components_shuffled_'
    else:
        file_pre = f'{basename}_{zdim}_components_'

    recon_file = os.path.join(train_dir, f'{file_pre}reconstruction.tsv')
    co_file = os.path.join(train_dir, f'{file_pre}sample_corr.tsv.gz')
    co_g_file = os.path.join(train_dir, f'{file_pre}gene_corr.tsv.gz')
    tybalt_hist_file = os.path.join(train_dir,
                                    f'{file_pre}tybalt_training_hist.tsv')
    adage_hist_file = os.path.join(train_dir,
                                   f'{file_pre}adage_training_hist.tsv')

    # Load Preprocessed Data --> could provide option to input raw data and just call process data from here,
    # but don't want to process multiple times, esp since we're running this independently on each zdim
    rnaseq_train_df = read_counts_or_params(input_train)
    rnaseq_test_df = read_counts_or_params(input_test)

    # Determine most variably expressed genes and subset
    if madfile is not None:
        mad_genes_df = read_counts_or_params(madfile)

        #data_base = os.path.join('..', '0.expression-download', 'data')
        #mad_file = os.path.join(data_base, '{}_mad_genes.tsv'.format(dataset))

        #mad_genes_df = pd.read_table(mad_file)
        mad_genes = mad_genes_df.iloc[0:num_mad_genes, ].index.astype(str)
        rnaseq_train_df = rnaseq_train_df.reindex(mad_genes, axis='columns')
        rnaseq_test_df = rnaseq_test_df.reindex(mad_genes, axis='columns')


# Initialize DataModel class

    dm = DataModel(df=rnaseq_train_df, test_df=rnaseq_test_df)
    dm.transform(how='zeroone') # data normalization happens here, don't need to feed in normalized data
# Set seed and list of algorithms for compression
    np.random.seed(1234)
    random_seeds = np.random.randint(0, high=1000000, size=num_seeds)

    #algorithms = ['pca', 'ica', 'nmf', 'dae', 'vae']

# Save population of models in specific folder
    comp_out_dir = os.path.join(out_dir, 'ensemble_z_matrices',
                                f'{basename}_components_{zdim}')
    if not os.path.exists(comp_out_dir):
        os.makedirs(comp_out_dir)

    reconstruction_results = []
    test_reconstruction_results = []
    sample_correlation_results = []
    tybalt_training_histories = []
    adage_training_histories = []
    for seed in random_seeds:

        np.random.seed(seed)

        seed_file = os.path.join(comp_out_dir, f'model_{seed}')

        if shuffle:
            seed_file = f'{seed_file}_shuffled'

            # randomly permute genes of each sample in the rnaseq matrix
            shuf_df = rnaseq_train_df.apply(lambda x:
                                            np.random.permutation(x.tolist()),
                                            axis=1)

            # Setup new pandas dataframe
            shuf_df = pd.DataFrame(shuf_df, columns=['gene_list'])
            shuf_df = pd.DataFrame(shuf_df.gene_list.values.tolist(),
                                   columns=rnaseq_train_df.columns,
                                   index=rnaseq_train_df.index)

            # Initiailze a new DataModel, with different shuffling each permutation
            dm = DataModel(df=shuf_df, test_df=rnaseq_test_df)
            dm.transform(how='zeroone')

        # Fit models
        if "pca" in algorithms:
            sys.stdout.write("\nrunning pca\n")
            dm.pca(n_components=zdim, transform_test_df=True)
        if "ics" in algorithms:
            sys.stdout.write("\nrunning ics\n")
            dm.ica(n_components=zdim, transform_test_df=True)
        if "nmf" in algorithms:
            sys.stdout.write("\nrunning nmf\n")
            dm.nmf(n_components=zdim, transform_test_df=True)

        if "vae" in algorithms:
            # run tybalt vae
            sys.stdout.write("\nrunning tybalt vae\n")
            dm.nn(n_components=zdim,
                  model='tybalt',
                  loss='binary_crossentropy',
                  epochs=int(paramsD['vae_epochs']),
                  batch_size=int(paramsD['vae_batch_size']),
                  learning_rate=float(paramsD['vae_lr']),
                  separate_loss=True,
                  verbose=True,
                  transform_test_df=True)

        if "dae" in algorithms:
            # run adage dae
            sys.stdout.write("\nrunning adage dae\n")
            dm.nn(n_components=zdim,
                  model='adage',
                  loss='binary_crossentropy',
                  epochs=int(paramsD['dae_epochs']),
                  batch_size=int(paramsD['dae_batch_size']),
                  learning_rate=float(paramsD['dae_lr']),
                  noise=float(paramsD['dae_noise']),
                  sparsity=float(paramsD['dae_sparsity']),
                  verbose=True,
                  transform_test_df=True)

        # Obtain z matrix (sample scores per latent space feature) for all models
        full_z_file = f'{seed_file}_z_matrix.tsv.gz'
        dm.combine_models().to_csv(full_z_file, sep='\t', compression='gzip')

        full_test_z_file = f'{seed_file}_z_test_matrix.tsv.gz'

        sys.stdout.write("combining trained models")
        dm.combine_models(test_set=True).to_csv(full_test_z_file, sep='\t',
                                                compression='gzip')

        # Obtain weight matrices (gene by latent space feature) for all models
        full_weight_file = f'{seed_file}_weight_matrix.tsv.gz'
        dm.combine_weight_matrix().to_csv(full_weight_file, sep='\t',
                                          compression='gzip')

        # Store reconstruction costs and reconstructed input at training end

        sys.stdout.write("compiling reconstruction costs")
        full_reconstruction, reconstructed_matrices = dm.compile_reconstruction()

        # Store reconstruction evaluation and data for test set
        full_test_recon, test_recon_mat = dm.compile_reconstruction(test_set=True)

        # Get correlations across samples and genes between input and output data
        pearson_corr = []
        spearman_corr = []
        pearson_corr_test = []
        spearman_corr_test = []

        sys.stdout.write("\ncalculating correlations across samples and genes\n")

        for algorithm in algorithms:
            # Training Sample Correlations
            sys.stdout.write(f"training: calculating pearson correlation for {algorithm}") # f string requires py >=3.6
            pearson_corr.append(
                get_recon_correlation(df=dm.df,
                                      recon_mat_dict=reconstructed_matrices,
                                      algorithm=algorithm,
                                      cor_type='pearson',
                                      genes=False)
                                )

            sys.stdout.write(f"training: calculating spearman correlation for {algorithm}") # f string requires py >=3.6
            spearman_corr.append(
                get_recon_correlation(df=dm.df,
                                      recon_mat_dict=reconstructed_matrices,
                                      algorithm=algorithm,
                                      cor_type='spearman',
                                      genes=False)
                                )

            # Testing Sample Correlations
            sys.stdout.write(f"testing: calculating pearson correlation for {algorithm}") # f string requires py >=3.6
            pearson_corr_test.append(
                get_recon_correlation(df=dm.test_df,
                                      recon_mat_dict=test_recon_mat,
                                      algorithm=algorithm,
                                      cor_type='pearson',
                                      genes=False)
                                    )
            sys.stdout.write(f"testing: calculating spearman correlation for {algorithm}") # f string requires py >=3.6
            spearman_corr_test.append(
                get_recon_correlation(df=dm.test_df,
                                      recon_mat_dict=test_recon_mat,
                                      algorithm=algorithm,
                                      cor_type='spearman',
                                      genes=False)
                                     )

        # Training - Sample correlations between input and reconstruction
        sys.stdout.write(f"training: calculating sample correlation for {algorithms}") # f string requires py >=3.6
        sample_correlation_results.append(
            compile_corr_df(pearson_list=pearson_corr,
                            spearman_list=spearman_corr,
                            algorithm_list=algorithms,
                            column_names=dm.df.index,
                            seed=seed,
                            data_type='training')
                            )

        # Testing - Sample correlations between input and reconstruction
        sys.stdout.write(f"testing: calculating sample correlation for {algorithms}")
        sample_correlation_results.append(
            compile_corr_df(pearson_list=pearson_corr_test,
                            spearman_list=spearman_corr_test,
                            algorithm_list=algorithms,
                            column_names=dm.test_df.index,
                            seed=seed,
                            data_type='testing')
                            )

        # Store training histories and intermediate results for neural networks
        reconstruction_results.append(
            full_reconstruction.assign(seed=seed, shuffled=shuffle)
            )
        test_reconstruction_results.append(
            full_test_recon.assign(seed=seed, shuffled=shuffle)
            )
        tybalt_training_histories.append(
            dm.tybalt_fit.history_df.assign(seed=seed, shuffle=shuffle)
            )
        adage_training_histories.append(
            dm.adage_fit.history_df.assign(seed=seed, shuffle=shuffle)
            )

# Save reconstruction and neural network training results
    pd.concat([
        pd.concat(reconstruction_results).assign(data_type='training'),
        pd.concat(test_reconstruction_results).assign(data_type='testing')
        ]).reset_index(drop=True).to_csv(recon_file, sep='\t', index=False)
    pd.concat(sample_correlation_results).to_csv(co_file,
                                                 sep='\t',
                                                 index=False,
                                                 float_format='%.3f',
                                                 compression='gzip')
    (pd.concat(tybalt_training_histories)
       .reset_index()
       .rename({'index': 'epoch'}, axis='columns')
       .to_csv(tybalt_hist_file, sep='\t', index=False))
    (pd.concat(adage_training_histories)
       .reset_index()
       .rename({'index': 'epoch'}, axis='columns')
       .to_csv(adage_hist_file, sep='\t', index=False))


if __name__ == '__main__':

    p = argparse.ArgumentParser()
    p.add_argument('input_train', help="input train data tsv")
    p.add_argument('input_test', help="input test data tsv")
    p.add_argument('--input_mad', help='filename containing genes sorted by median absolute deviation', default=None)
    p.add_argument('--basename', help= "basename to use in file names")
    p.add_argument('-z', '--zdim', help='dimensionality of z. prev called num_components')
    p.add_argument('-p', '--paramsfile',
                        help='text file optimal hyperparameter assignment for z')
    p.add_argument('-o', '--outdir', help='where to save the output files')
    p.add_argument('-s', '--num_seeds', default=5,
                        help='number of different seeds to run on current data')
    p.add_argument('-r', '--shuffle', action='store_true',
                        help='randomize gene expression data for negative control')
    p.add_argument('-m', '--num_mad_genes', default=8000,
                        help='subset num genes based on median absolute deviation')
    p.add_argument('--algorithms', nargs='+', default= ['pca', 'ica', 'nmf', 'dae', 'vae'],
                        help='compression algorithms to run')
    args = p.parse_args()

    zsweep_params = read_params(args.paramsfile)

    zsweep_paramsD = zsweep_params.to_dict()

    # check that the chosen zdims exist in the paramsD
    zdim = str(args.zdim)
    component_error = f"{zdim} is not found in {zsweep_paramsD} - either add it to the file or choose a different number of components"
    #assert str(num_components) in param_df.columns, component_error
    assert zdim in zsweep_params.keys(), component_error

    # start by just allowing a single zdim? Call train models once for each zdim.
    paramsD = zsweep_paramsD[zdim] #subset to just this z dimension

    sys.exit(train_models(args.basename, args.input_train, args.input_test, int(zdim), paramsD, args.outdir, int(args.num_seeds), args.shuffle, args.input_mad, int(args.num_mad_genes), algorithms))


