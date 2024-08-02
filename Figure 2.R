library(cancereffectsizeR)
library(data.table)
library(readxl)
library(dplyr)
library(ces.refset.hg19)
library(stringr)
library(ggrepel)

# Set Working Directory
setwd("/Users/andrew/Desktop/Summer/Project/Code")


Sample_Information <- read.delim("Sample_Information.txt")

setnames(Sample_Information, "Sample.type", "Sample_type")

# Keep only rows with "Primary" and "Metastasis":
filtered_samples <- Sample_Information %>% 
  filter(Sample_type %in% c("Primary", "Metastasis"))

# Remove duplicates based on "Unique_Patient_Identifier":
Sample_type <- filtered_samples %>% 
  distinct(Unique_Patient_Identifier, .keep_all = TRUE)

# Check for consistency in "Sample_type" for each "Unique_Patient_Identifier"
consistent_samples <- filtered_samples %>% 
  group_by(Unique_Patient_Identifier) %>% 
  filter(n_distinct(Sample_type) == 1) %>% 
  ungroup() %>% 
  distinct(Unique_Patient_Identifier, .keep_all = TRUE)

# Identify inconsistent samples
inconsistent_samples <- filtered_samples %>% 
  group_by(Unique_Patient_Identifier) %>% 
  filter(n_distinct(Sample_type) > 1) %>% 
  ungroup()

print("Consistent samples:")
print(consistent_samples)

print("Inconsistent samples:")
print(inconsistent_samples)

## preparing the data
# Paper Data

paper_tc_data <- fread("Paper_TC_Data.txt", skip = 1)

# Refset determined using IGV
tc_maf <- preload_maf(maf = paper_tc_data, refset = "ces.refset.hg19", 
                      sample_col = "Sample ID",  start_col = "Position",
                      ref_col = "Reference", tumor_allele_col = "Alternate")

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

# 3. CESAnalysis Creation and General Results ----
###creating CESAnalysis and loading data 
cesa <- CESAnalysis(refset = "ces.refset.hg19")

# Filter Variants
top_tgs_genes <- c("TP53", "PIK3CA", "TERT", "NF1", "NF2", "NRAS", "BRAF", "CDKN2A", "CDKN2B", 
                   "NKX2-1","RET", "KMT2C", "KMT2D", "BCOR", "TBX3", "PTEN", "EIF1AX", "RBM10", 
                   "ATM", "ARID1A")

tgs_coverage <- ces.refset.hg19$gr_genes[ces.refset.hg19$gr_genes$names %in% top_tgs_genes]

# Loading MAF into CESAnalysis
cesa <- load_maf(cesa = cesa, maf = tc_maf, coverage = "targeted", maf_name = "THCA",
                 covered_regions = tgs_coverage, covered_regions_name = "top_genes",
                 covered_regions_padding = 10)
cesa <- load_maf(cesa = cesa, maf = tcga_maf, coverage = "genome", maf_name = "TCGA_THCA")
cesa <- load_maf(cesa, maf = genie_maf, maf_name = "Genie_THCA", coverage = "targeted",
                 covered_regions = tgs_coverage, covered_regions_name = "top_genes",
                 covered_regions_padding = 10)

# Loading Clinical Data
cesa <- load_sample_data(cesa, tcga_clinical)
cesa <- load_sample_data(cesa, genie_clinical)
cesa <- load_sample_data(cesa, Sample_type)

#defining groups:
primary_samples <- consistent_samples[consistent_samples$Sample_type == "Primary", ]
primary_groups <- cesa$samples[Unique_Patient_Identifier %in% primary_samples$Unique_Patient_Identifier]

meta_samples <- consistent_samples[consistent_samples$Sample_type == "Metastasis", ]
meta_groups <- cesa$samples[Unique_Patient_Identifier %in% meta_samples$Unique_Patient_Identifier]

cesa <- gene_mutation_rates(cesa = cesa, covariates = "THCA", samples = primary_groups, save_all_dndscv_output = T)
cesa <- gene_mutation_rates(cesa = cesa, covariates = "THCA", samples = meta_groups, save_all_dndscv_output = T)

selected_genes <- top_tgs_genes

RefCDS = ces.refset.hg19$RefCDS
dndscv_gene_names <- cesa$gene_rates$gene
nsyn_sites = sapply(RefCDS[dndscv_gene_names], function(x) colSums(x[["L"]])[1])

# selecting mutation rate data for samples in primary_groups
samples_in_primary_groups <- length(unique(cesa$dNdScv_results$rate_grp_1$annotmuts$sampleID))

# selecting mutation rate data for samples in metastasis_groups
samples_in_metastasis_groups <- length(unique(cesa$dNdScv_results$rate_grp_2$annotmuts$sampleID))

library(tidyverse)

### creating a data frame with mutation rate data for primary_groups and metastasis_groups
mut_rate_df <- tibble(gene = cesa$dNdScv_results$rate_grp_1$genemuts$gene_name,
                      exp_primary_mu = cesa$dNdScv_results$rate_grp_1$genemuts$exp_syn_cv,
                      exp_metastasis_mu = cesa$dNdScv_results$rate_grp_2$genemuts$exp_syn_cv)

mut_rate_df$n_syn_sites = nsyn_sites[mut_rate_df$gene]

mut_rate_df %>% 
  mutate(primary_mu = (exp_primary_mu / n_syn_sites) / samples_in_primary_groups) %>%
  mutate(metastasis_mu = (exp_metastasis_mu / n_syn_sites) / samples_in_metastasis_groups) %>%
  mutate(cancer_greater = metastasis_mu > primary_mu) -> 
  mut_rate_df

# defining rate 1 and rate 2 as mutation rates for primary_groups and metastasis_groups
rate_1 <- mut_rate_df|>
  select(gene, primary_mu)
rate_2 <- mut_rate_df|>
  select(gene, metastasis_mu)

# change in mutation rate across stages
mut_rate_df <- mut_rate_df %>% 
  select(gene, primary_mu, metastasis_mu) %>% 
  mutate(p_1 = primary_mu / metastasis_mu) %>% 
  mutate(p_2 = 1 - p_1)

# saving "last" gene mutation rates into separate data frame, "last" rates meaning from last stage Metastasis_mu
set_cancer_rates <- mut_rate_df %>%
  select(gene, metastasis_mu) %>%
  data.table::setDT()

# clear the gene rates in the cesa object 
cesa <- clear_gene_rates(cesa = cesa)

# setting gene rates to highest rates from Metastasis_mu
setnames(set_cancer_rates, "metastasis_mu", "rate")

cesa <- set_gene_rates(cesa = cesa, rates = set_cancer_rates, missing_genes_take_nearest = T) 

# infer trinculeotide-context-specific relative rates of SNV mutation from a mutational signature analysis
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "THCA", 
                                                            treatment_naive = TRUE)

# estimating trinucleotide mutation rates
cesa <- trinuc_mutation_rates(cesa = cesa, signature_set = "COSMIC_v3.2", 
                                                signature_exclusions = signature_exclusions)

# defining compound variants
compound <- define_compound_variants(cesa = cesa, 
                                     variant_table = cesa$variants |>
                                       filter(intergenic == F, gene %in% selected_genes),
                                     by = "gene", merge_distance = Inf)

source("new_sequential_lik.R")

setnames(set_cancer_rates, "rate", "metastasis_mu")

cesa$samples <- cesa$samples %>%
  left_join(consistent_samples %>% select(Unique_Patient_Identifier, Sample_type), by = "Unique_Patient_Identifier")

for(comp_ind in 1:length(compound)){
  
  this_comp <- compound[comp_ind, ]
  
  this_gene <- unlist(unique(this_comp$snv_info$genes))
  these_props <- mut_rate_df[mut_rate_df$gene == this_gene,c("p_1","p_2")]
  these_props <- c(these_props$p_1, these_props$p_2)
  
  cesa <- ces_variant(cesa = cesa, variants = this_comp, model = sequential_lik_dev, 
                      ordering_col = "Sample_type",
                                        ordering = c('Primary', 'Metastasis'),
                                        lik_args = list(sequential_mut_prop = these_props), run_name = this_gene)
  
}

plot_effects(effects = cesa$selection$ARID1A, group_by = "variant")
