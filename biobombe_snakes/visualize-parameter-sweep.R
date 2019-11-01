
library(ggplot2)

set.seed(123)

`%>%` <- dplyr::`%>%`

source(file.path('viz_utils.R'))

adage_file <- file(snakemake@input[["adage"]])
tybalt_file <- file(snakemake@input[["tybalt"]])
dataset_name <- snakemake@params[["dataset_name"]]

# Set output directories
output_fig <- file.path("figures", "viz_results")

# Load Data
tybalt_data <- readParamSweep(tybalt_file, algorithm = 'Tybalt')

tybalt_param_z_png <- file.path(output_fig, "z_parameter_final_loss_tybalt.png")
tybalt_param_z_pdf <- file.path(output_fig, "z_parameter_final_loss_tybalt.pdf")

p <- plotFinalLoss(tybalt_data$select_df, algorithm = 'Tybalt', dataset = dataset_name)
p

ggsave(tybalt_param_z_png, plot = p, height = 3, width = 5.5)
ggsave(tybalt_param_z_pdf, plot = p, height = 3, width = 5.5)

tybalt_one_model <- tybalt_data$full_df %>% dplyr::filter(learning_rate == "0.0005", batch_size == "50", epochs == 100)

tybalt_one_model$num_components <-
  dplyr::recode_factor(tybalt_one_model$num_components, 
                       `5` = "Latent Dim: 5",
                       `25` = "Latent Dim: 25", 
                       `50` = "Latent Dim: 50",
                       `75` = "Latent Dim: 75",
                       `100` = "Latent Dim: 100",
                       `125` = "Latent Dim: 125")

tybalt_one_model_png <- file.path(output_fig, "z_parameter_tybalt_training.png")
tybalt_one_model_pdf <- file.path(output_fig, "z_parameter_tybalt_training.pdf")

p <- plotOneModel(tybalt_data$one_model_df, algorithm = 'Tybalt', dataset = dataset_name)
p

ggsave(tybalt_one_model_png, plot = p, height = 2.5, width = 5)
ggsave(tybalt_one_model_pdf, plot = p, height = 2.5, width = 5)

tybalt_best_params <- tybalt_data$best_params
tybalt_best_params

best_param_file <- file.path("results" , "z_latent_dim_best_tybalt_params.tsv")
readr::write_tsv(tybalt_best_params, best_param_file)

tybalt_good_training_df <- tybalt_data$melt_df %>%
  dplyr::filter(batch_size == 50, epochs == 100, kappa == "0.0") %>%
  dplyr::filter(
    (learning_rate == "0.002" & num_components == 5) |
      (learning_rate == "0.0015" & num_components == 25) |
      (learning_rate == "0.0015" & num_components == 50) |
      (learning_rate == "0.0015" & num_components == 75) |
      (learning_rate == "0.001" & num_components == 100) |
      (learning_rate == "0.0005" & num_components == 125)
  )

# Reorder the latent space dimensionality for plotting
num_com <- tybalt_good_training_df$num_components
num_com <- factor(num_com, levels = sort(as.numeric(paste(unique(num_com)))))
tybalt_good_training_df$num_components <- num_com

tybalt_best_model_png <- file.path(output_fig, "z_parameter_best_model_tybalt.png")
tybalt_best_model_pdf <- file.path(output_fig, "z_parameter_best_model_tybalt.pdf")

p <- plotBestModel(tybalt_good_training_df, algorithm = 'Tybalt', dataset = dataset_name, output_fig_dir=output_fig)
p

ggsave(tybalt_best_model_png, plot = p, height = 2.5, width = 4)
ggsave(tybalt_best_model_pdf, plot = p, height = 2.5, width = 4)

# Load Data
adage_data <- readParamSweep(adage_file, algorithm = "ADAGE")

# Specify that the model has tied weights
adage_param_z_png <- file.path(output_fig, "z_parameter_final_loss_adage.png")
adage_param_z_pdf <- file.path(output_fig, "z_parameter_final_loss_adage.pdf")

p <- plotFinalLoss(adage_data$select_df, algorithm = 'ADAGE', dataset = dataset_name)
p

ggsave(adage_param_z_png, plot = p, height = 2.5, width = 5.5)
ggsave(adage_param_z_pdf, plot = p, height = 2.5, width = 5.5)

# Several hyperparameter combinations did not converge
# This was particularly a result of the low learning rates - filter and replot
adage_data$select_df <- adage_data$select_df %>%
    dplyr::filter(end_loss < 0.01, learning_rate != 'Learn: 1e-05')

adage_sub_png <- file.path(output_fig, "z_parameter_final_loss_remove_converge_adage.png")
adage_sub_pdf <- file.path(output_fig, "z_parameter_final_loss_remove_converge_adage.pdf")

p <- plotFinalLoss(adage_data$select_df, algorithm = 'ADAGE', dataset = dataset_name)
p

ggsave(adage_sub_png, plot = p, height = 2.5, width = 5.5)
ggsave(adage_sub_pdf, plot = p, height = 2.5, width = 5.5)

adage_best_params <- adage_data$best_params
adage_best_params

best_param_file <- file.path("results" , "z_latent_dim_best_adage_params.tsv")
readr::write_tsv(adage_best_params, best_param_file)

adage_tied_good_training_df <- adage_data$melt_df %>%
  dplyr::filter(sparsity == "0.0",
                epochs == 100,
                batch_size == 50,
                noise == "0.0") %>%
  dplyr::filter(
    (num_components == 5 & learning_rate == "0.0015") |
      (num_components == 25 & learning_rate == "0.0015") |
      (num_components == 50 & learning_rate == "0.0005") |
      (num_components == 75 & learning_rate == "0.0005") |
      (num_components == 100 & learning_rate == "0.0005") |
      (num_components == 125 & learning_rate == "0.0005"))

num_com <- adage_tied_good_training_df$num_components
num_com <- factor(num_com, levels = sort(as.numeric(paste(unique(num_com)))))
adage_tied_good_training_df$num_components <- num_com

best_model_png <- file.path(output_fig, "z_parameter_best_model_adage.png")
best_model_pdf <- file.path(output_fig, "z_parameter_best_model_adage.pdf")

p <- plotBestModel(adage_tied_good_training_df, algorithm = 'ADAGE', dataset = dataset_name)
p

ggsave(best_model_png, plot = p, height = 2.5, width = 4)
ggsave(best_model_pdf, plot = p, height = 2.5, width = 4)

# Load Data
#gtex_tybalt <- readParamSweep(gtex_tybalt_file, algorithm = 'tybalt')

#tybalt_param_z_png <- file.path(gtex_fig, "z_parameter_final_loss_tybalt_GTEX.png")
#tybalt_param_z_pdf <- file.path(gtex_fig, "z_parameter_final_loss_tybalt_GTEX.pdf")

#p <- plotFinalLoss(gtex_tybalt$select_df, algorithm = 'Tybalt', dataset = 'GTEx')
#p

#ggsave(tybalt_param_z_png, plot = p, height = 3, width = 5.5)
#ggsave(tybalt_param_z_pdf, plot = p, height = 3, width = 5.5)

#gtex_tybalt_best_params <- gtex_tybalt$best_params
#gtex_tybalt_best_params

#best_param_file <- file.path("results" , "z_latent_dim_best_tybalt_params_GTEX.tsv")
#readr::write_tsv(gtex_tybalt_best_params, best_param_file)

#gtex_tybalt_good_training_df <- gtex_tybalt$melt_df %>%
#  dplyr::filter(epochs == 100, kappa == 0.5) %>%
#  dplyr::filter(
#    (learning_rate == "0.0025" & batch_size == 100 & num_components == 5) |
#      (learning_rate == "0.0025" & batch_size == 100 & num_components == 25) |
#      (learning_rate == "0.002" & batch_size == 100 & num_components == 50) |
#      (learning_rate == "0.002" & batch_size == 50 & num_components == 75) |
#      (learning_rate == "0.0015" & batch_size == 50 & num_components == 100) |
#      (learning_rate == "0.0015" & batch_size == 50 & num_components == 125)
#  )

# Reorder the latent space dimensionality for plotting
#num_com <- gtex_tybalt_good_training_df$num_components
#num_com <- factor(num_com, levels = sort(as.numeric(paste(unique(num_com)))))
#gtex_tybalt_good_training_df$num_components <- num_com

#best_model_png <- file.path(gtex_fig, "z_parameter_best_model_tybalt_GTEX.png")
#best_model_pdf <- file.path(gtex_fig, "z_parameter_best_model_tybalt_GTEX.pdf")

#p <- plotBestModel(gtex_tybalt_good_training_df, algorithm = 'Tybalt', dataset = 'GTEx')
#p

#ggsave(best_model_png, plot = p, height = 2.5, width = 4)
#ggsave(best_model_pdf, plot = p, height = 2.5, width = 4)

# Load Data
#gtex_adage <- readParamSweep(gtex_adage_file, algorithm = 'adage')

#adage_param_z_png <- file.path(gtex_fig, "z_parameter_final_loss_adage_GTEX.png")
#adage_param_z_pdf <- file.path(gtex_fig, "z_parameter_final_loss_adage_GTEX.pdf")

#p <- plotFinalLoss(gtex_adage$select_df, algorithm = 'ADAGE', dataset = 'GTEx')
#p

#ggsave(adage_param_z_png, plot = p, height = 3, width = 5.5)
#ggsave(adage_param_z_pdf, plot = p, height = 3, width = 5.5)

# Several hyperparameter combinations did not converge
# This was particularly a result of the low learning rates - filter and replot
#gtex_adage$select_df <- gtex_adage$select_df %>%
#    dplyr::filter(end_loss < 0.01, learning_rate != 'Learn: 1e-05')

#adage_sub_png <- file.path(gtex_fig, "z_parameter_final_loss_remove_converge_adage_GTEX.png")
#adage_sub_pdf <- file.path(gtex_fig, "z_parameter_final_loss_remove_converge_adage_GTEX.pdf")

#p <- plotFinalLoss(gtex_adage$select_df, algorithm = 'ADAGE', dataset = 'GTEx')
#p

#ggsave(adage_sub_png, plot = p, height = 2.5, width = 5.5)
#ggsave(adage_sub_pdf, plot = p, height = 2.5, width = 5.5)

#gtex_adage_best_params <- gtex_adage$best_params
#gtex_adage_best_params

#best_param_file <- file.path("results" , "z_latent_dim_best_adage_params_GTEX.tsv")
#readr::write_tsv(gtex_adage_best_params, best_param_file)

#gtex_adage_tied_good_training_df <- gtex_adage$melt_df %>%
#  dplyr::filter(sparsity == "0.0",
#                epochs == 100,
#                batch_size == 50) %>%
#  dplyr::filter(
#    (num_components == 5 & learning_rate == "0.001" & noise == "0.1") |
#      (num_components == 25 & learning_rate == "0.001" & noise == "0.0") |
#      (num_components == 50 & learning_rate == "0.0005" & noise == "0.0") |
#      (num_components == 75 & learning_rate == "0.0005" & noise == "0.0") |
#      (num_components == 100 & learning_rate == "0.0005" & noise == "0.0") |
#      (num_components == 125 & learning_rate == "0.0005" & noise == "0.0"))

# Reorder the latent space dimensionality for plotting
#num_com <- gtex_adage_tied_good_training_df$num_components
#num_com <- factor(num_com, levels = sort(as.numeric(paste(unique(num_com)))))
#gtex_adage_tied_good_training_df$num_components <- num_com

#best_model_png <- file.path(gtex_fig, "z_parameter_best_model_adage_GTEX.png")
#best_model_pdf <- file.path(gtex_fig, "z_parameter_best_model_adage_GTEX.pdf")

#p <- plotBestModel(gtex_adage_tied_good_training_df, dataset = 'GTEx', algorithm = 'ADAGE')
#p

#ggsave(best_model_png, plot = p, height = 2.5, width = 4)
#ggsave(best_model_pdf, plot = p, height = 2.5, width = 4)

# Load Data
#target_tybalt <- readParamSweep(target_tybalt_file, algorithm = 'tybalt')

#tybalt_param_z_png <- file.path(target_fig, "z_parameter_final_loss_tybalt_TARGET.png")
#tybalt_param_z_pdf <- file.path(target_fig, "z_parameter_final_loss_tybalt_TARGET.pdf")

#p <- plotFinalLoss(target_tybalt$select_df, algorithm = 'Tybalt', dataset = 'TARGET')
#p

#ggsave(tybalt_param_z_png, plot = p, height = 3, width = 5.5)
#ggsave(tybalt_param_z_pdf, plot = p, height = 3, width = 5.5)

#target_tybalt_best_params <- target_tybalt$best_params
#target_tybalt_best_params

#best_param_file <- file.path("results" , "z_latent_dim_best_tybalt_params_TARGET.tsv")
#readr::write_tsv(tcga_adage_best_params, best_param_file)

#target_tybalt_good_training_df <- target_tybalt$melt_df %>%
#  dplyr::filter(batch_size == 25,
#                epochs == 100,
#                kappa == '0.5') %>%
#  dplyr::filter(
#    (num_components == 5 & learning_rate == "0.0015") |
#      (num_components == 25 & learning_rate == "0.0015") |
#      (num_components == 50 & learning_rate == "0.0015") |
#      (num_components == 75 & learning_rate == "0.0015") |
#      (num_components == 100 & learning_rate == "0.0015") |
#      (num_components == 125 & learning_rate == "0.0005"))

# Reorder the latent space dimensionality for plotting
#num_com <- target_tybalt_good_training_df$num_components
#num_com <- factor(num_com, levels = sort(as.numeric(paste(unique(num_com)))))
#target_tybalt_good_training_df$num_components <- num_com

#best_model_png <- file.path(target_fig, "z_parameter_best_model_tybalt_TARGET.png")
#best_model_pdf <- file.path(target_fig, "z_parameter_best_model_tybalt_TARGET.pdf")

#p <- plotBestModel(target_tybalt_good_training_df, algorithm = 'Tybalt', dataset = 'TARGET')
#p

#ggsave(best_model_png, plot = p, height = 2.5, width = 4)
#ggsave(best_model_pdf, plot = p, height = 2.5, width = 4)

# Load Data
#target_adage <- readParamSweep(target_adage_file, algorithm = 'adage')

#adage_param_z_png <- file.path(target_fig, "z_parameter_final_loss_adage_TARGET.png")
#adage_param_z_pdf <- file.path(target_fig, "z_parameter_final_loss_adage_TARGET.pdf")

#p <- plotFinalLoss(target_adage$select_df, algorithm = 'ADAGE', dataset = 'TARGET')
#p

#ggsave(adage_param_z_png, plot = p, height = 3, width = 5.5)
#ggsave(adage_param_z_pdf, plot = p, height = 3, width = 5.5)

#target_adage_best_params <- target_adage$best_params
#target_adage_best_params

#best_param_file <- file.path("results" , "z_latent_dim_best_adage_params_TARGET.tsv")
#readr::write_tsv(tcga_adage_best_params, best_param_file)

#target_adage_tied_good_training_df <- target_adage$melt_df %>%
#  dplyr::filter(sparsity == "0.0",
#                epochs == 100,
#                batch_size == 50,
#                noise == "0.1",
#                learning_rate == "0.0005")

# Reorder the latent space dimensionality for plotting
#num_com <- target_adage_tied_good_training_df$num_components
#num_com <- factor(num_com, levels = sort(as.numeric(paste(unique(num_com)))))
#target_adage_tied_good_training_df$num_components <- num_com

#best_model_png <- file.path(target_fig, "z_parameter_best_model_adage_TARGET.png")
#best_model_pdf <- file.path(target_fig, "z_parameter_best_model_adage_TARGET.pdf")

#p <- plotBestModel(target_adage_tied_good_training_df, algorithm = 'ADAGE', dataset = 'TARGET')
#p

#ggsave(best_model_png, plot = p, height = 2.5, width = 4)
#ggsave(best_model_pdf, plot = p, height = 2.5, width = 4)
