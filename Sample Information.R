library(cancereffectsizeR)
library(data.table)
library(readxl)
library(dplyr)
library(ces.refset.hg19)
library(stringr)
library(ggplot2)
library(ggrepel)
library(MutationalPatterns)
library(tidyverse)
library(patchwork)
library(openxlsx)

# Set Working Directory
setwd("/Users/andrew/Desktop/Summer/Project/Code")

## preparing the data
# Paper Data

paper_tc_data <- fread("Paper_TC_Data.txt", skip = 1)

# Refset determined using IGV
tc_maf <- preload_maf(maf = paper_tc_data, refset = "ces.refset.hg19", 
                      sample_col = "Sample ID",  start_col = "Position",
                      ref_col = "Reference", tumor_allele_col = "Alternate", 
                      keep_extra_columns =  TRUE)

# Loading TCGA Data ----
tcga_maf_file <- "TCGA-THCA.maf.gz"
if (!file.exists(tcga_maf_file)) {
  get_TCGA_project_MAF(project = "THCA", filename = tcga_maf_file)
}

tcga_clinical <- fread("TCGA_Clinical.txt")
setnames(tcga_clinical, "case_id", "Unique_Patient_Identifier")

names(tcga_clinical)[which(names(tcga_clinical) == "residual_disease")[2]] <- "residual_disease_2"

tcga_maf <- preload_maf(maf = tcga_maf_file, 
                        chain_file = "hg38ToHg19.over.chain", 
                        refset = "ces.refset.hg19")

# Loading Genie Data ----
genie <- fread("GENIE_Mutation_Data.txt")
genie_clinical <- fread("GENIE_Clinical.txt", skip = 4)
setnames(genie_clinical, "PATIENT_ID", "Unique_Patient_Identifier")
genie_clinical <- genie_clinical %>% filter(`CANCER_TYPE` == "Thyroid Cancer")
genie_sample <- unique(genie_clinical$SAMPLE_ID)
genie_maf <- preload_maf(maf = genie, refset = "ces.refset.hg19")
genie_maf <- genie_maf %>% filter(`Unique_Patient_Identifier` %in% c(genie_sample))

# Keep samples where column Problem is equal to NA:
tc_maf <- tc_maf %>% filter(is.na(problem))
tcga_maf <- tcga_maf %>% filter(is.na(problem))
genie_maf <- genie_maf %>% filter(is.na(problem))

# keeping only samples that do not occur at germline variant sites:
tc_maf <- tc_maf %>% filter (`germline_variant_site` == FALSE)
tcga_maf <- tcga_maf %>% filter (`germline_variant_site` == FALSE)
genie_maf <- genie_maf %>% filter (`germline_variant_site` == FALSE)

# keeping only samples that do not occur in repetitive regions 
tc_maf <- tc_maf %>% filter (`repetitive_region` == FALSE | cosmic_site_tier %in% 1:3)
tcga_maf <- tcga_maf %>% filter (`repetitive_region` == FALSE | cosmic_site_tier %in% 1:3)
genie_maf <- genie_maf %>% filter (`repetitive_region` == FALSE | cosmic_site_tier %in% 1:3)

# keeping snv:
tc_maf <- subset(tc_maf, variant_type == "snv")
tcga_maf <- subset(tcga_maf, variant_type == "snv")
genie_maf <- subset(genie_maf, variant_type == "snv")

# TXT File - Column 1: Sample ID | Column 2: Primary / Metastasis ----
sample_info <- tc_maf %>% select(`Unique_Patient_Identifier`, `Sample type`)

setnames(tcga_maf, "Unique_Patient_Identifier", "case_submitter_id", skip_absent=TRUE)
tcga <- full_join(tcga_maf, tcga_clinical, by = "case_submitter_id")
tcga <- tcga %>% select(`case_submitter_id`, `ajcc_pathologic_m`)
setnames(tcga, "case_submitter_id", "Unique_Patient_Identifier", skip_absent=TRUE)
setnames(tcga, "ajcc_pathologic_m", "Sample type", skip_absent=TRUE)

sample_info <- rbind(sample_info, tcga)

genie_full <- full_join(genie_maf, genie_clinical, by = "Unique_Patient_Identifier")
genie_full <- genie_full %>% select(`Unique_Patient_Identifier`, `SAMPLE_TYPE`)
setnames(genie_full, "SAMPLE_TYPE", "Sample type", skip_absent=TRUE)

sample_info <- rbind(sample_info, genie_full)

sample_info$`Sample type`[sample_info$`Sample type` == "M0"] <- "Primary"
sample_info$`Sample type`[sample_info$`Sample type` == "M1"] <- "Metastasis"
sample_info$`Sample type`[sample_info$`Sample type` == "MX" |
                            sample_info$`Sample type` == "'--" |
                            sample_info$`Sample type` == "Not Applicable or Heme" |
                            sample_info$`Sample type` == "Unspecified" |
                            sample_info$`Sample type` == "Not Collected"] <- NA

duplicates <- sample_info %>% group_by(`Unique_Patient_Identifier`) %>%
  filter(n_distinct(`Sample type`) > 1) %>% pull(`Unique_Patient_Identifier`) %>% unique()

sample_info <- sample_info %>% filter(!`Unique_Patient_Identifier` %in% duplicates)

sample_info <- sample_info %>% distinct(Unique_Patient_Identifier, .keep_all = TRUE)

write.xlsx(sample_info, file = "Sample Information.xlsx")
