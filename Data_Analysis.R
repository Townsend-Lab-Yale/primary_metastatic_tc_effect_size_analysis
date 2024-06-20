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

# Refset determined using IGV
tc_maf <- preload_maf(maf = paper_tc_data, refset = "ces.refset.hg19", 
                      sample_col = "Sample ID",  start_col = "Position",
                      ref_col = "Reference", tumor_allele_col = "Alternate", 
                      keep_extra_columns = TRUE)

# Filtering the data down to the top 20 frequent genes reported by the research
tc_maf <- tc_maf %>% filter(`Symbol` %in% c("TP53", "PIK3CA", "TERT", "NF1", "NF2", "NRAS", 
                                   "BRAF", "CDKN2A", "CDKN2B", "NKX2-1","RET", "KMT2C", 
                                   "KMT2D", "BCOR", "TBX3", "PTEN", "EIF1AX", "RBM10", 
                                   "ATM", "ARID1A"))

# 3. Checking Various Metrics ----
library(dplyr)
# # Checking observations
# paper_tc_data %>% summarise(count = n_distinct(`Sample ID`))
# # Checking if primary with metastases and no metastases are differentiated
# paper_tc_data %>% summarise(count = n_distinct(`Sample type`))
# # Check number of genes
# unique(tc_maf$`Symbol`)
# 
# # Checking which samples are missing
# unique(tc_maf$`Thyroid cancer subtype`)
# tc_check = distinct(tc_maf, `Unique_Patient_Identifier`, .keep_all = TRUE)
# table(tc_check$`Thyroid cancer subtype`)

# 4. CESAnalysis Creation ----
# Creating CESAnalysis
cesa <- CESAnalysis(refset = "ces.refset.hg19")

# Loading MAF into CESAnalysis
cesa <- load_maf(cesa = cesa, maf = tc_maf, coverage = "genome", maf_name = "TC")

# Improve accuracy of signature extraction by excluding signatures that can 
# be safely presumed to not be present in the samples
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "TC", 
                                                            treatment_naive = TRUE)

# Adding information about snv_counts, raw_attributions, biological_weights and trinuc_rates
# to the CESAnalysis
cesa <- trinuc_mutation_rates(cesa,
                              signature_set = ces.refset.hg19$signatures$COSMIC_v3,
                              signature_exclusions = signature_exclusions, 
                              assume_identical_mutational_processes = TRUE,
                              sig_averaging_threshold = 0)

# Estimating regional rates of mutation in the absence of selection
cesa <- gene_mutation_rates(cesa, covariates = ces.refset.hg19$covariates$thyroid)

# Including an optional run_name
cesa <- ces_variant(cesa = cesa, run_name = "recurrents")

# 5. General Results Across all TC Types ----
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
cesa$epistasis$gene_epistasis_example
  # 5.1 Plotting the General Metastatic Model ----
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