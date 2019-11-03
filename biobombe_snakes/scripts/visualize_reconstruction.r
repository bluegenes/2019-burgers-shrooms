# modified for snakemake by ntpierce, 2019

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))

# Load helper functions
source(file.path("scripts", "util.R"))

# Define the dataset to compile results for
#dataset <- 'TARGET'
basename <- snakemake@params[["basename"]]
base_dir <- file.path("figures", dataset) # make a snakemake param, out_dir

# compile_reconstruction data crawls through given folder structure 
# to obtain dataset specific reconstruction results

results_dir_base <- snakemake@params[["results_dir"]]
target_recon_cost_df <- compile_reconstruction_data(basename, results_dir_base) 
# Write results to file
recon_cost_file <- file.path(snakemake@output[["recon_cost"]])
readr::write_tsv(target_recon_cost_df, path = recon_file)

# plot the reconstruction loss
target_recon_gg <- plot_reconstruction_loss(target_recon_cost_df)

#target_path <- file.path(base_dir, paste0("reconstruction_cost_", dataset))
target_path <- file.path(snakemake@output[["recon_cost_fig"]]

save_png_pdf(p = target_recon_gg,
             path_prefix = target_path,
             height = 70,
             width = 170)

target_recon_gg

# Compile VAE specific reconstruction loss
target_vae_recon_cost_df <- compile_reconstruction_data(basename, results_dir_base, data_focus = "vae")
target_vae_loss_gg <- plot_vae_training(target_vae_recon_cost_df)

#target_path <- file.path(base_dir, paste0("vae_training_reconstruction_", dataset))
target_path <- file.path(snakemake@output[["vae_recon_cost_fig"]]

save_png_pdf(p = target_vae_loss_gg,
             path_prefix = target_path,
             height = 130,
             width = 100)

target_vae_loss_gg

######
### here down is same, except the other datasets? double check!


# Subset to iterations that may have converged
target_recon_cost_df <- target_recon_cost_df %>% dplyr::filter(reconstruction_cost < 4000) #5000?

target_recon_filter_gg <- plot_reconstruction_loss(target_recon_cost_df)

filter_path <- file.path(base_dir, paste0("reconstruction_cost_subset_converge_", basename))

save_png_pdf(p = target_recon_filter_gg,
             path_prefix = filter_path,
             height = 70,
             width = 170)

target_recon_filter_gg

# Subset to testing non-shuffled data
target_recon_cost_df <- target_recon_cost_df %>%  #SHOULD THIS BE PRE OR POST FILTERING (ABOVE)? I think post, based on "converge" in filename below
    dplyr::filter(data_type == 'testing', shuffled == 'False')

target_recon_filter_test_gg <- plot_reconstruction_loss(target_recon_cost_df)

non_shuffle_path <- file.path(base_dir, paste0("reconstruction_cost_subset_converge_testing_", basename))

save_png_pdf(p = target_recon_filter_test_gg,
             path_prefix = non_shuffle_path,
             height = 70,
             width = 170)

target_recon_filter_test_gg

# Remove shuffled data and replot
target_vae_recon_cost_df <- target_vae_recon_cost_df %>% dplyr::filter(shuffle == "False")

target_vae_loss_filter_test_gg <- plot_vae_training(tcga_vae_recon_cost_df)

target_path <- file.path(base_dir, paste0("vae_training_reconstruction_subset_converge_", basename))

save_png_pdf(p = tcga_vae_loss_filter_test_gg,
             path_prefix = tcga_path,
             height = 130,
             width = 100)

tcga_vae_loss_filter_test_gg


#### below plots all three datasets together --> don't need. 
# do we already have this plot for each individual? Or do I need to rework it for that?

### ORRR can I rework it to aggregate multiple dataset plots (diff script, diff snakemake rule?)


legend <- get_legend(target_recon_gg) 

main_plot <- (
    cowplot::plot_grid(
        #gtex_recon_filter_test_gg + ggtitle('GTEX') + xlab('') +
        #    theme(plot.margin = margin(t = 0.5, r = 0.2, b = 0, l = 0.4),
        #          legend.position = "none",
        #          panel.grid.major = element_line(size = 0.25),
        #          panel.grid.minor = element_line(size = 0.175)),
        #tcga_recon_filter_test_gg + ggtitle('TCGA') + xlab('') +
        #    theme(plot.margin = margin(t = 0, r = 0.2, b = 0, l = 0.4),
        #          legend.position = "none",
        #          panel.grid.major = element_line(size = 0.25),
        #          panel.grid.minor = element_line(size = 0.175)),
        target_recon_gg + ggtitle('TARGET') +
            theme(plot.margin = margin(t = 0, r = 0.2, b = 0.3, l = 0.4),
                  legend.position = "none",
                  panel.grid.major = element_line(size = 0.25),
                  panel.grid.minor = element_line(size = 0.175)),
        labels = c("a" )#, "b", "c"),
        ncol = 1,
        nrow = 1 #3
    )
)

main_plot = cowplot::plot_grid(main_plot, legend, rel_widths = c(1, 0.15), ncol = 2)
main_plot

main_path <- file.path("figures", "reconstruction_summary")

save_png_pdf(p = main_plot,
             path_prefix = main_path,
             height = 130,
             width = 170)
