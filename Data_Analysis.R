# 1. Installation and Preparation ----
# Some dependencies are large, so we increase the download time limit to be safe
# Already Done - Don't Have to do Every Time
# Press Cancel for Restart
# options(timeout = 600)
# install.packages("remotes")
# 
# Download cancereffectsizeR Package (Already Done - Don't Have to do Every Time)
# remotes::install_github("Townsend-Lab-Yale/cancereffectsizeR", dependencies = TRUE)
# 
# Obtain Reference Dataset (Already Done - Don't Have to do Every Time)
# options(timeout = 600)
# remotes::install_github("Townsend-Lab-Yale/ces.refset.hg19@*release")

# Just in case need to clear environment
# rm(list=ls())

# Set Working Directory
setwd("/Users/andrew/Desktop/Summer/Project/Code")

# Import Packages
library(cancereffectsizeR)
library(data.table)

# 2. Data Loading ----
# The file is in Excel
# Only need to install once
# install.packages("readxl")

library(readxl)

# Loading the data, skipping the first line because it is file description
# File contains the non-synonymous somatic mutations identified in thyroid cancers (TC) 
# subjected to targeted massively parallel sequencing
paper_tc_data <- read_excel("TC_Data.xlsx", skip = 1)

tc_maf <- preload_maf(maf = paper_tc_data, refset = "ces.refset.hg19", 
                      sample_col = "Sample ID",  start_col = "Position",
                      ref_col = "Reference", tumor_allele_col = "Alternate", 
                      keep_extra_columns = TRUE)


# Checking various metrics
library(dplyr)
# Checking observations
paper_tc_data %>% summarise(count = n_distinct(`Sample ID`))
# Checking if primary with metastases and no metastases are differentiated
paper_tc_data %>% summarise(count = n_distinct(`Sample type`))
# Check number of genes
unique(tc_maf$`Symbol`)

# Checking which samples are missing
unique(tc_maf$`Thyroid cancer subtype`)
tc_check = distinct(tc_maf, `Unique_Patient_Identifier`, .keep_all = TRUE)
table(tc_check$`Thyroid cancer subtype`)

