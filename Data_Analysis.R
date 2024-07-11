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
# remotes::install_github("Townsend-Lab-Yale/ces.refset.hg38@*release")

# Just in case need to clear environment
rm(list=ls())

# Set Working Directory
setwd("/Users/andrew/Desktop/Summer/Project/Code")

# Import Packages
library(cancereffectsizeR)
library(data.table)
library(dplyr)
library(stringr)

# 2. Data Loading and Cleaning ----
# The file is in Excel
# Only need to install once
# install.packages("readxl")
library(readxl)
library(readr)

# Loading the data, skipping the first line because it is file description
# File contains the non-synonymous somatic mutations identified in thyroid cancers (TC) 
# subjected to targeted massively parallel sequencing 
# Paper Data ----

paper_tc_data <- read_excel("TC_Data.xlsx", skip = 1)

# Refset determined using IGV
tc_maf <- preload_maf(maf = paper_tc_data, refset = "ces.refset.hg19", 
                      sample_col = "Sample ID",  start_col = "Position",
                      ref_col = "Reference", tumor_allele_col = "Alternate")

# Keep samples where column Problem is equal to NA:
tc_maf <- tc_maf %>% filter(is.na(problem))

# keeping only samples that do not occur at germline variant sites:
tc_maf <- tc_maf %>% filter (`germline_variant_site` == FALSE)

# keeping only samples that do not occur in repetitive regions 
tc_maf <- tc_maf %>% filter (`repetitive_region` == FALSE | cosmic_site_tier %in% 1:3)

# keeping snv:
tc_maf <- subset(tc_maf, variant_type == "snv")

# Loading TCGA Data ----
tcga_maf_file <- "TCGA-THCA.maf.gz"
if (!file.exists(tcga_maf_file)) {
  get_TCGA_project_MAF(project = "THCA", filename = tcga_maf_file)
}

tcga_clinical <- fread("clinical.tsv")

setnames(tcga_clinical, "case_id", "Unique_Patient_Identifier")

tcga_maf <- preload_maf(maf = tcga_maf_file, refset = "ces.refset.hg38")

tcga_maf <- tcga_maf %>% filter(is.na(problem))

tcga_maf <- tcga_maf %>% filter (`germline_variant_site` == FALSE)

tcga_maf <- tcga_maf %>% filter (`repetitive_region` == FALSE | cosmic_site_tier %in% 1:3)

tcga_maf <- subset(tcga_maf, variant_type == "snv")

# Loading Genie Data ----
genie_tc <- fread("data_mutations_extended.txt")

genie_clinical <- fread("data_clinical_sample.txt", skip = 4)

setnames(genie_clinical, "PATIENT_ID", "Unique_Patient_Identifier")

genie_clinical <- genie_clinical %>% filter(`CANCER_TYPE` == "Thyroid Cancer")

genie_sample <- unique(genie_clinical$SAMPLE_ID)

genie_maf <- preload_maf(maf = genie_tc, refset = "ces.refset.hg19")

genie_maf <- genie_maf %>% filter(`Unique_Patient_Identifier` %in% c(genie_sample))

genie_maf <- genie_maf %>% filter(is.na(problem))

genie_maf <- genie_maf %>% filter (`germline_variant_site` == FALSE)

genie_maf <- genie_maf %>% filter (`repetitive_region` == FALSE | cosmic_site_tier %in% 1:3)

genie_maf <- subset(genie_maf, variant_type == "snv")

# 3. CESAnalysis Creation and General Results ----
# Creating CESAnalysis
cesa <- CESAnalysis(refset = "ces.refset.hg19")

# Filter Variants
top_tgs_genes <- c("TP53", "PIK3CA", "TERT", "NF1", "NF2", "NRAS", "BRAF", "CDKN2A", "CDKN2B", 
               "NKX2-1","RET", "KMT2C", "KMT2D", "BCOR", "TBX3", "PTEN", "EIF1AX", "RBM10", 
               "ATM", "ARID1A")

tgs_coverage <- ces.refset.hg19$gr_genes[ces.refset.hg19$gr_genes$gene %in% top_tgs_genes]

# Loading MAF into CESAnalysis
cesa <- load_maf(cesa = cesa, maf = tc_maf, coverage = "genome", maf_name = "THCA")
# cesa <- load_maf(cesa = cesa, maf = tcga_maf, coverage = "genome", maf_name = "TCGA_THCA")
cesa <- load_maf(cesa, maf = genie_maf, maf_name = "Genie_THCA", coverage = "targeted",
                 covered_regions = tgs_coverage, covered_regions_name = "top_genes",
                 covered_regions_padding = 50)

# Loading Clinical Data
# cesa <- load_sample_data(cesa, tcga_clinical)
cesa <- load_sample_data(cesa, genie_clinical)

# We'll use all suggested exclusions (TCGA primary tumors are treatment-naive)
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "THCA",
                                                            treatment_naive = TRUE)

# Adding information about snv_counts, raw_attributions, biological_weights and trinuc_rates
# to the CESAnalysis
cesa <- trinuc_mutation_rates(cesa,
                              signature_set = ces.refset.hg19$signatures$COSMIC_v3.2,
                              signature_exclusions = signature_exclusions,
                              sig_averaging_threshold = 0,
                              assume_identical_mutational_processes = TRUE)

# Estimating regional rates of mutation in the absence of selection
cesa <- gene_mutation_rates(cesa, covariates = ces.refset.hg19$covariates$THCA)

# Including an optional run_name
cesa <- ces_variant(cesa = cesa, run_name = "recurrents")

# Plotting the most selected variants
plot_effects(effects = cesa$selection$recurrents)

# Effects of all recurrent variants across the most selected genes
plot_effects(cesa$selection$recurrents,
             group_by = "gene", label_individual_variants = FALSE)

# Epistasis Model on Gene Level
genes <- c("NRAS", "BRAF", "PIK3CA")

# Get consensus covered regions
combined_coverage <- intersect(cesa$coverage_ranges$exome$`exome+`, cesa$coverage_ranges$targeted$top_genes)

# Get variants in the genes of interest that have sequencing coverage in all samples
variants <- select_variants(cesa, genes = genes, gr = combined_coverage)

cesa <- ces_gene_epistasis(cesa = cesa, genes = genes, variants = variants, run_name = "gene_epistasis_example")

cesa <- ces_epistasis(cesa = cesa, variants = variants, run_name = "epistasis")

plot_effects(cesa$selection$gene_epistasis_example)

  # 3.1 Plotting the General Epistatic Model ----
  require(grid) # Need the grid package for this plot
  results <- cesa$epistasis$gene_epistasis_example
  results <- results[, .(
    v1 = variant_A, v2 = variant_B, ces_A0, ces_B0, ces_A_on_B,
    ces_B_on_A, p_A_change, p_B_change, p_epistasis
  )]
  
  # By change, we mean fold-change of selection on mutant background over wildtype background
  results[, change_in_v2 := ces_B_on_A / ces_B0]
  results[, change_in_v1 := ces_A_on_B / ces_A0]
  
  # Put in desired order for display
  results <- results[, pairname := paste(v1, v2, sep = ".")]
  
  results[, x := 1:.N]
  results[, v1_x := x - .2]
  results[, v2_x := x + .2]
  results[, alpha := .6]
  results[p_epistasis < .05, alpha := 1]
  
  results[, v1_signif := ""]
  results[p_A_change < .05, v1_signif := "*"]
  results[p_A_change < .01, v1_signif := "**"]
  results[p_A_change < .001, v1_signif := "***"]
  results[, v1_signif_y := change_in_v1 + (.13 * sign(change_in_v1 - 1))]
  
  results[, v2_signif := ""]
  results[p_B_change < .05, v2_signif := "*"]
  results[p_B_change < .01, v2_signif := "**"]
  results[p_B_change < .001, v2_signif := "***"]
  results[, v2_signif_y := change_in_v2 + (.13 * sign(change_in_v2 - 1))]
  
  x_labels <- unlist(S4Vectors::zipup(results$v1, results$v2))
  x_label_pos <- unlist(S4Vectors::zipup(results$v1_x, results$v2_x))
  
  # Have to get fancy to depict significance nicely in legend.
  draw_signif_key <- function(data, params, size) {
    grobTree(
      rectGrob(
        x = .25, y = .5, width = .5, height = 1,
        gp = gpar(col = NA, fill = alpha("plum4", data$alpha), lty = data$linetype)
      ),
      rectGrob(
        x = .75, y = .5, width = .5, height = 1,
        gp = gpar(col = NA, fill = alpha("sandybrown", data$alpha), lty = data$linetype)
      )
    )
  }
  
  ggplot(data = results) +
    # Put in a reference line depicting no change in selection
    geom_hline(yintercept = 1, color = "darkgrey") +
    geom_rect(aes(xmin = v1_x - .2, xmax = v1_x + .2, ymin = 1, ymax = change_in_v1, fill = "v1", alpha = alpha),
              show.legend = c(alpha = FALSE, fill = TRUE)
    ) +
    geom_rect(aes(xmin = 1, xmax = 1, ymin = 0, ymax = 0, alpha = alpha),
              show.legend = c(alpha = TRUE, fill = FALSE), key_glyph = draw_signif_key
    ) +
    geom_text(aes(x = v1_x, y = v1_signif_y, label = v1_signif), size = 7) +
    scale_alpha_identity(
      breaks = c(1, .6), labels = c("Significant", "Not significant"),
      guide = guide_legend(
        title = "Pairwise epistasis", override.aes = list(fill = "sandybrown", alpha = c(1, .6)),
        order = 1
      )
    ) +
    geom_rect(aes(xmin = v2_x - .2, xmax = v2_x + .2, ymin = 1, ymax = change_in_v2, fill = "v2", alpha = alpha),
              show.legend = c(alpha = FALSE, fill = TRUE)
    ) +
    geom_text(aes(x = v2_x, y = v2_signif_y, label = v2_signif), size = 7) +
    
    # Build legend
    scale_fill_manual(
      name = "Ratio of selection",
      breaks = c("v1", "v2"),
      labels = list(
        expression(frac("gene 1 on mutated gene 2", "gene 1 on wildtype gene 2")),
        expression(frac("gene 2 on mutated gene 1", "gene 2 on wildtype gene 1"))
      ),
      values = c("v1" = "plum4", "v2" = "sandybrown"),
      guide = guide_legend(label.theme = element_text(size = 6.5))
    ) +
    scale_x_continuous(breaks = x_label_pos, labels = x_labels) +
    scale_y_continuous(breaks = seq(from = 0, to = 3, by = .25)) +
    xlab("Gene pair") +
    ylab("Ratio of selection coefficients") +
    theme_classic() +
    theme(
      legend.position = "bottom", legend.title = element_text(size = 10),
      axis.ticks.length.x = unit(0, "cm")
    )
