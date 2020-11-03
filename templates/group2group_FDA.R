# group2group_FDA.R

# This script demonstrates group comparison functional data analysis using the
# bayesFDA package. As the simplest form of valuable FDA, this does NOT
# account for repeated measures, covariates, or other ANOVA components. However,
# it does handle an arbitrary number of groups and all desired  pairwise
# comparisons are handled efficiently via one model.

#-------------------------------------------------------------------------------
# Setup Information #
#-------------------------------------------------------------------------------

# Provide paths to data files.
# Data from each group need to be placed in different files. Each csv file is
# organized such that a column is an observed curve. Prior time normalization is
# currently required, thus each csv must have the same number of rows. Missing
# values or NA's are only allowed if the entire column is missing, other than
# header information like subject label or treatment. Add as many as desired,
# note that order will be used later.

data_paths <- c(
  groupA = "path/to/group1.csv",
  groupB = "path/to/group2.csv",
  groupR = "another/path/to/data3.csv"
)

# Provide the number of header rows. This must be common for all files.

num_header_rows = 3

# Provide desired comparisons as a matrix of couplets. Order matters as the
# second value will be subtracted from the first. Provide numbers corresponding
# to the ordering of the above paths. For example, to compare group R minus
# group A in the example above, you would enter c(3,1) as a couplet. Do no edit
# the other arguments.

comparison_couplets <- matrix(
  c(2,1),
  c(3,1)
, ncol = 2, byrow = TRUE)

# Provide a path and filenames for output to be saved.

