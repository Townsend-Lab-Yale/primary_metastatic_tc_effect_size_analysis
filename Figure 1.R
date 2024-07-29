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

# Trinucleotide Mutation Profile for Primary and Metastatic Tumors ----
primary_cesa <- CESAnalysis(refset = "ces.refset.hg19")

primary_samples <- tc_maf %>% filter(`Sample type` == "Primary")
primary_cesa <- load_maf(cesa = primary_cesa, maf = primary_samples, 
                         coverage = "targeted", maf_name = "THCA", covered_regions = tgs_coverage, covered_regions_name = "top_genes",
                         covered_regions_padding = 10)

primary_tcga_clinical <- tcga_clinical %>% filter(`ajcc_pathologic_m` == "M0")
primary_tcga <- unique(primary_tcga_clinical$case_submitter_id)
primary_tcga_maf <- tcga_maf %>% filter(`Unique_Patient_Identifier` %in% c(primary_tcga))
primary_cesa <- load_maf(cesa = primary_cesa, maf = primary_tcga_maf)
primary_cesa <- load_sample_data(primary_cesa, primary_tcga_clinical)

primary_genie_clinical <- genie_clinical %>% filter(`SAMPLE_TYPE` == "Primary")
primary_genie <- unique(primary_genie_clinical$SAMPLE_ID)
primary_genie_maf <- genie_maf %>% filter(`Unique_Patient_Identifier` %in% c(primary_genie))

top_tgs_genes <- c("TP53", "PIK3CA", "TERT", "NF1", "NF2", "NRAS", "BRAF", "CDKN2A", "CDKN2B", 
                   "NKX2-1","RET", "KMT2C", "KMT2D", "BCOR", "TBX3", "PTEN", "EIF1AX", "RBM10", 
                   "ATM", "ARID1A")

tgs_coverage <- ces.refset.hg19$gr_genes[ces.refset.hg19$gr_genes$names %in% top_tgs_genes]

primary_cesa <- load_maf(primary_cesa, maf = primary_genie_maf, 
                         coverage = "targeted",
                         covered_regions = tgs_coverage, covered_regions_name = "topgenes",
                         covered_regions_padding = 10)
primary_cesa <- load_sample_data(primary_cesa, primary_genie_clinical)

signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "THCA",
                                                            treatment_naive = TRUE)

primary_cesa <- trinuc_mutation_rates(primary_cesa,
                                      signature_set = ces.refset.hg19$signatures$COSMIC_v3.2,
                                      signature_exclusions = signature_exclusions,
                                      sig_averaging_threshold = 0,
                                      assume_identical_mutational_processes = TRUE)

# Extract trinucleotide mutation rates from CESAnalysis object
primary_trinuc_rates <- primary_cesa$trinuc_rates
primary_trinuc_rates <- primary_trinuc_rates[, -1] 

# Prepare the data for plotting
primary_trinuc_rates_long <- primary_trinuc_rates %>%
  pivot_longer(cols = everything(), 
               names_to = "trinucleotide", 
               values_to = "rate")

# Plot the trinucleotide mutation rates using ggplot2
primary_trinuc_rate_plot <- ggplot(primary_trinuc_rates_long, 
                                   aes(x = trinucleotide, y = rate, 
                                       fill = substr(primary_trinuc_rates_long$trinucleotide, 1, 1))) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Trinucleotide Mutation Rates in Primary Tumors",
       x = NULL,
       y = "Trinucleotide Rates") + guides(fill = guide_legend(title = "First Base")) +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank())

# Metastasis ----
meta_cesa <- CESAnalysis(refset = "ces.refset.hg19")

meta_samples <- tc_maf %>% filter(`Sample type` == "Metastasis")
meta_cesa <- load_maf(cesa = meta_cesa, maf = meta_samples, 
                      coverage = "targeted", maf_name = "THCA",
                      covered_regions = tgs_coverage, covered_regions_name = "top_genes",
                      covered_regions_padding = 10)

meta_tcga_clinical <- tcga_clinical %>% filter(`ajcc_pathologic_m` == "M1")
meta_tcga <- unique(meta_tcga_clinical$case_submitter_id)
meta_tcga_maf <- tcga_maf %>% filter(`Unique_Patient_Identifier` %in% c(meta_tcga))
meta_cesa <- load_maf(cesa = meta_cesa, maf = meta_tcga_maf)
meta_cesa <- load_sample_data(meta_cesa, meta_tcga_clinical)

meta_genie_clinical <- genie_clinical %>% filter(`SAMPLE_TYPE` == "Metastasis")
meta_genie <- unique(meta_genie_clinical$SAMPLE_ID)
meta_genie_maf <- genie_maf %>% filter(`Unique_Patient_Identifier` %in% c(meta_genie))

top_tgs_genes <- c("TP53", "PIK3CA", "TERT", "NF1", "NF2", "NRAS", "BRAF", "CDKN2A", "CDKN2B", 
                   "NKX2-1","RET", "KMT2C", "KMT2D", "BCOR", "TBX3", "PTEN", "EIF1AX", "RBM10", 
                   "ATM", "ARID1A")

tgs_coverage <- ces.refset.hg19$gr_genes[ces.refset.hg19$gr_genes$names %in% top_tgs_genes]

meta_cesa <- load_maf(meta_cesa, maf = meta_genie_maf, 
                      coverage = "targeted",
                      covered_regions = tgs_coverage, covered_regions_name = "topgenes",
                      covered_regions_padding = 10)
meta_cesa <- load_sample_data(meta_cesa, meta_genie_clinical)

signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "THCA",
                                                            treatment_naive = TRUE)

meta_cesa <- trinuc_mutation_rates(meta_cesa,
                                   signature_set = ces.refset.hg19$signatures$COSMIC_v3.2,
                                   signature_exclusions = signature_exclusions,
                                   sig_averaging_threshold = 0,
                                   assume_identical_mutational_processes = TRUE)

# Extract trinucleotide mutation rates from CESAnalysis object
meta_trinuc_rates <- meta_cesa$trinuc_rates
meta_trinuc_rates <- meta_trinuc_rates[, -1] 

# Prepare the data for plotting
meta_trinuc_rates_long <- meta_trinuc_rates %>%
  pivot_longer(cols = everything(), 
               names_to = "trinucleotide", 
               values_to = "rate")

# Plot the trinucleotide mutation rates using ggplot2
meta_trinuc_rate_plot <- ggplot(meta_trinuc_rates_long, aes(x = trinucleotide, y = rate, 
                                                            fill = substr(meta_trinuc_rates_long$trinucleotide, 1, 1))) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Trinucleotide Mutation Rates in Metastatic Tumors",
       x = "Trinucleotide Variant",
       y = "Trinucleotide Rates") + guides(fill = guide_legend(title = "First Base")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Figure Merge----
combined_plot <- primary_trinuc_rate_plot / meta_trinuc_rate_plot + 
  plot_layout(guides = 'collect')

print(combined_plot)
