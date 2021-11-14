#### covariate_FDA.R ####

# This script performs functional data analysis involving group comparisons 
# while potentially controlling for one or more non-functional covariates. 
# This is done through Bayesian P-spline model using the bayesFDA package 
# found at github.com/wzhorton/bayesFDA. 

#-------------------------------------------------------------------------------
# Data Loading #
#-------------------------------------------------------------------------------

# There are 3 data components for this model and each should be given in a
# separate csv file:
#   - Observed curves/waveforms
#   - Individual level covariates
#   - Trial level covariates
# In this section you will specify paths and column names for data to be loaded

#-----------------#
# Observed Curves #
#-----------------#
# This script expects the waveform file to be organized into columns. That is,
# a single observed curve is contained in one column. Additionally, the top row
# should contain the ID number corresponding to the individual. Further, missing
# entries are not allowed in any part - ID numbers need to be repeated if 
# applicable, missing trials/replicates must be deleted, and curves must be
# time normalized to contain the same number of points. Combine all groups into
# one file. Group demarcation is handled in the covariate script.

waveform_path <- "path/to/curves.csv"

#-----------------------------#
# Individual Level Covariates #
#-----------------------------#
# This file is expected to have at least 2 columns: an ID number column and a 
# group label column. More columns may be used for additional covariates. The
# ID numbers should never repeat and must match the ID's found in the waveform
# file. Note that this file is always required. 

individual_covs_path <- "path/to/indvcovs.csv"
group_column_label <- "grp"
id_column_label_indv <- "ID"

#----------------------------------#
# Trial/Replicate Level Covariates #
#----------------------------------#
# If the data contain multiple observations per individual, this file may 
# be used to provide replicate level covariates. This file is not required if 
# no trial level covariates are present, or if the data structure does not 
# involve subject replication. Leave the entry below NULL if so. If including, 
# an ID column must be provided which repeats according to the number of 
# replications. A trial/replicate number may also be added, but will be removed 
# during analysis. Do not include group membership here.

trial_covs_path <- "path/to/trialcovs.csv" # NULL if not applicable
id_column_label_trial <- "ID"
trial_number_column_label <- "Trial" # NULL if not included (not required)


#-------------------------------------------------------------------------------
# Comparison Specification #
#-------------------------------------------------------------------------------

# Here you provide the group comparisons you want the model to output. This is
# done by providing the group names in pairs. The group names here must match
# those found in the group column in the individual covariate file. Each pair
# is formatted as c("A","B"), which would produce a comparison of "B minus A".
# Add as many comparisons as you would like.

comparison_list <- list(
  c("group 2", "group 1"),
  c("group 2", "group 3")
)


#-------------------------------------------------------------------------------
# Covariate Baseline Values #
#-------------------------------------------------------------------------------

# Here you provide a baseline value for each covariate. This baseline is 
# important for interpreting the group averages. Common choices include generally
# accepted averages, values of interest, zero, and sample means. The sample mean 
# will automatically be computed given an entry of "avg". The names must match
# those given in the covariate files.

individual_covs_baselines <- list(
  weight = 80.5,
  height = "avg"
)

trial_covs_baselines <- list( # NULL if not applicable
  time = "avg",
  hope = 0
)


#-------------------------------------------------------------------------------
# Output Direction #
#-------------------------------------------------------------------------------

# Here you specify output parameters. Make sure these are right as any mistakes
# will likely involve having to run the full analysis again.

#-------------------------------#
# Output Destination and Naming #
#-------------------------------#
# Give both the folder path and a common name to append to each output file.

output_folder_path <- "path/to/folder"
output_prefix <- "VGRF"

#------------------#
# Output Selection #
#------------------#
# Choices include model output for comparisons, covariates, group averages, and 
# significant regions. To select, set to TRUE.

output_comparisons <- TRUE
output_covariates <- TRUE
output_group_avgs <- TRUE
output_sig_regions <- TRUE


#-------------------------------------------------------------------------------
# MCMC Parameters #
#-------------------------------------------------------------------------------

# These govern behavior surrounding model fit. The default values are generally
# fine given the incoming data are smooth. Slight modifications to these will
# have little impact on overall analysis, but can help with computation speed
# at the cost of some potential accuracy.

num_spline_bases <- 20
mcmc_iterations <- 10000
mcmc_burnin <- round(mcmc_iterations/10)


#-------------------------------------------------------------------------------
# Package Installation Notes #
#-------------------------------------------------------------------------------

# The bayesFDA package is hosted through GitHub at github.com/wzhorton/bayesFDA.
# The easiest method to install is via the Devtools package:
#   - Make sure you have the package "devtools" installed.
#   - Run the command devtools::github_install("wzhorton/bayesFDA")
# If the install fails, you are likely on MacOS. If so, follow these steps:
#   - Go to github.com/wzhorton/bayesFDA/releases/latest
#   - Download the bayesFDA MacOS binary file
#   - Within R run install.packages("path/to/binaryfile",repos=NULL)


#-------------------------------------------------------------------------------
# Main Script (Modify at your own risk) #
#-------------------------------------------------------------------------------

library(bayesFDA)

error <- FALSE
error_msg <- ""

#---------------------------------#
# Read, Check, and Aggregate Data #
#---------------------------------#

# Waveforms
curves <- read.csv(waveform_path, header = FALSE)
id_curves <- as.numeric(curves[1,])
curves <- curves[-1,]

if(any(id_curves%%1!=0)){
  error <- TRUE
  error_msg <- "ERROR: The curve ID numbers appear to contain decimals."
} else if(!error && any(is.na(curves))) {
  error <- TRUE
  error_msg <- "ERROR: Missing values found in waveforms"
}

# Individual Covariates
indv_covs <- read.csv(individual_covs_path)
indv_cov_labels <- colnames(indv_covs)

if(!error && !(group_column_label %in% indv_cov_labels)){
  error <- TRUE
  error_msg <- "ERROR: Group column label is not found in individual covariate file"
} else if(!error && !(id_column_label_indv %in% indv_cov_labels)){
  error <- TRUE
  error_msg <- "ERROR: ID column label is not found in individual covariate file"
}

id_indv <- indv_covs[,id_column_label_indv]
grp_indv <- as.factor(indv_covs[,group_column_label])
indv_covs[,id_column_label_indv] <- indv_covs[,group_column_label] <- NULL

if(!error && !setequal(id_indv,id_curves)){
  error <- TRUE
  error_msg <- "ERROR: Individual covariate IDs do not match waveform IDs"
} else if(!error && !setequal(names(individual_covs_baselines),colnames(indv_covs))){
  error <- TRUE
  error_msg <- "ERROR: Individual covariate names do not match baseline list"
}

for(i in 1:length(individual_covs_baselines)){
  cov_base <- individual_covs_baselines[[i]]
  cov_name <- names(individual_covs_baselines)[i]
  if(cov_base == "avg"){
    cov_base = mean(indv_covs[,cov_name])
  }
  indv_covs[,cov_name] <- indv_covs[,cov_name] - cov_base
}

grp_design_mat <- setNames(data.frame(matrix(as.numeric(model.matrix(~-1+grp_indv)), 
                                             ncol = length(levels(grp_indv)))), levels(grp_indv))
indv_covs <- cbind(indv_covs, grp_design_mat)

# Trial Covariates
if(!is.null(trial_covs_path)){
  trial_covs <- read.csv(trial_covs_path)
  trial_cov_labels <- colnames(trial_covs)
  
  if(!error && !(id_column_label_trial %in% trial_cov_labels)){
    error <- TRUE
    error_msg <- "ERROR: ID column label is not found in trial covariate file"
  }
  
  id_trial <- trial_covs[,id_column_label_trial]
  trial_covs[,id_column_label_trial] <- NULL
  
  if(!is.null(trial_number_column_label)){
    if(!error && !(trial_number_column_label %in% trial_cov_labels)){
      error <- TRUE
      error_msg <- 'ERROR: Trial number column is not found in trial covariates'
    }
    trial_covs[,trial_number_column_label] <- NULL
  }
  
  if(!error && !setequal(id_trial,id_curves)){
    error <- TRUE
    error_msg <- "ERROR: Trial covariate IDs do not match waveform IDs"
  } else if(!error && !setequal(names(trial_covs_baselines),colnames(trial_covs))){
    error <- TRUE
    error_msg <- "ERROR: Trial covariate names do not match baseline list"
  }
  
  for(i in 1:length(trial_covs_baselines)){
    cov_base <- trial_covs_baselines[[i]]
    cov_name <- names(trial_covs_baselines)[i]
    if(cov_base == "avg"){
      cov_base = mean(indv_covs[,cov_name])
    }
    trial_covs[,cov_name] <- trial_covs[,cov_name] - cov_base
  }
}

# Aggregate
cov_list <- lapply(unique(id_curves), function(id){
  ind_row <- indv_covs[which(id_indv == id), , drop = FALSE]
  if(!is.null(trial_covs_path)){
    trial_mat <- trial_covs[which(id_trial == id), , drop = FALSE]
    ind_mat <- ind_row[rep(1,nrow(trial_mat)),]
    return(cbind(trial_mat, ind_mat))
  } else {
    ind_mat <- ind_row[rep(1,sum(id_curves == id)), , drop = FALSE]
    return(ind_mat)
  }
})
cov_mat <- do.call(rbind, cov_list)

#-------------------------------#
# Fit Model and Produce Results #
#-------------------------------#

if(error){
  stop(error_msg)
} else {
  fit <- cov_fda(as.matrix(curves), as.matrix(cov_mat), p = num_spline_bases, 
                 niter = mcmc_iterations, nburn = mcmc_burnin)
  
  grp_labs <- names(grp_design_mat)
  cov_labs <- setdiff(names(cov_mat),names(grp_design_mat))
  all_labs <- names(cov_mat)
  basis <- splines::bs(seq(0,1,len = nrow(curves)), df = num_spline_bases, 
                       intercept = FALSE)
  
  setwd(output_folder_path)
  if(output_group_avgs){
    for(g in grp_labs){
      bchain <- basis%*%fit$B[match(g, all_labs),,]
      out <- cbind(apply(bchain, 1, function(bq) quantile(bq, 0.025)),
                   rowMeans(bchain),
                   apply(bchain, 1, function(bq) quantile(bq, 0.975)))
      write.csv(out, paste0(output_prefix,"-",g,"-GroupMean.csv"))
      if(output_sig_regions){
        write.csv(extract_sig_area(out), paste0(output_prefix,"-",g,"-GroupMean-Areas.csv"))
      }
    }
  }
  if(output_covariates){
    for(cc in cov_labs){
      bchain <- basis%*%fit$B[match(cc, all_labs),,]
      out <- cbind(apply(bchain, 1, function(bq) quantile(bq, 0.025)),
                   rowMeans(bchain),
                   apply(bchain, 1, function(bq) quantile(bq, 0.975)))
      write.csv(out, paste0(output_prefix,"-",cc,"-Effect.csv"))
      if(output_sig_regions){
        write.csv(extract_sig_area(out), paste0(output_prefix,"-",cc,"-Effect-Areas.csv"))
      }
    }
  }
  if(output_comparisons){
    for(i in 1:length(comparison_list)){
      g1 <- comparison_list[[i]][1]
      g2 <- comparison_list[[i]][2]
      bchain <- basis%*%(fit$B[match(g2, all_labs),,]-fit$B[match(g1, all_labs),,])
      out <- cbind(apply(bchain, 1, function(bq) quantile(bq, 0.025)),
                   rowMeans(bchain),
                   apply(bchain, 1, function(bq) quantile(bq, 0.975)))
      write.csv(out, paste0(output_prefix,"-",g2,"-minus-",g1,"-DiffMean.csv"))
      if(output_sig_regions){
        write.csv(extract_sig_area(out), paste0(output_prefix,"-",g2,"-minus-",g1,"-DiffMean-Areas.csv"))
      }
    }
  }
}

