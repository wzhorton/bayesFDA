# covariate_FDA.R #

# This script demonstrates a type of functional data analysis which incorporates
# non-functional covariates. This form of FDA can entirely superceed the 
# group to group analysis through indicator functions, albeit some conveniences
# are lacking here.

#-------------------------------------------------------------------------------
# Setup Information #
#-------------------------------------------------------------------------------

# Provide paths to data files.
# Data are generally expected in 3 files, however current functionality is 
# restricted to 2. The curve file or waveform file is organized such that
# each column is an observed (time normalized) waveform. The covariate file
# is organized such that each column corresponds to a different covariate. 
# An important connection between the two files is that the rows of the covariate
# file must correspond to the columns of the curve file. For example, individual
# 1 must have their waveform given in column 1 of the curve file and their 
# covariate observations in row 1 of the covariate file. In the future this 
# will be accomplished by ID matching, but right now it is done strictly based
# on order (do not include ID as a covariate).

curve_file_path <- "path/to/waveforms.csv"
covariate_file_path <- "path/to/covs.csv"

# Header rows are expected. Any number are allowed on the curve file, as long as
# you specify how many here. The covariate file expects exactly 1 header row 
# containing the name of the covariates.

curve_header_rows <- 3

# Provide a path to a folder where output will be saved. 

output_path = "path/to/folder"

# The final user inputs are model fitting parameters. The defaults are generally
# safe assuming the data are smooth and there are at least 2*num_spline_basis
# points per curve.

num_spline_basis = 25
num_mcmc_iterations = 10000
num_mcmc_burnin = round(num_mcmc_iterations*0.1)

#-------------------------------------------------------------------------------
# Package Installation Comments #
#-------------------------------------------------------------------------------

# To install the package, go to https://github.com/wzhorton/bayesFDA/releases/latest
# and download the right file for your operating system. Once downloaded, open
# R and run the following command:
# install.packages("path/to/downloaded/file.tar.gz", repos = NULL)

#-------------------------------------------------------------------------------
# Main Script (Modify at your own risk) #
#-------------------------------------------------------------------------------

library(bayesFDA)

#----------
# Read data
#----------
curves <- read.csv(curve_file_path, skip = curve_header_rows, header = FALSE)
covs <- read.csv(covariate_file_path)

#----------
# Fit model
#----------
model_fit <- cov_fda(as.matrix(curves), cbind(1,as.matrix(covs)), num_spline_basis, 
                     num_mcmc_iterations, num_mcmc_burnin)

#---------------
# Process output
#---------------
ncrv <- nrow(covs)
npts <- nrow(curves)

basis <- splines::bs(seq(0,1,len=npts), df = num_spline_basis, intercept = TRUE)




