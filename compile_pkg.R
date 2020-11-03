#### compile_pkg.R ####
# This file is meant for recompilation during development.
# For regular installation use devtools functionality.
# Use build_vignettes = TRUE in devtools arguments to get html files

# Load needed libraries
library(roxygen2)
library(devtools)

# To include data, run this code then add constructive code to the file and run
# To update data, run the file again (not this code)
#library(usethis)
#use_data_raw("IndicatorList")


# Clean up session
rm(list = ls())
gc()


# pkgbuild::compile_dll() # run if new C functions don't appear
document() # devtools; creates namespace, documentation files, and compiles C


# setting build_vignettes = FALSE may increase speed of install
# if built with vignettes, use browseVignettes("pkgname") to view.
install(build_vignettes = FALSE) # devtools; installation complete
