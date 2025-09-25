library(cancereffectsizeR)
library(data.table)
library(readxl)
library(dplyr)
library(ces.refset.hg19)
library(stringr)
library(ggrepel)
library(MutationalPatterns)
library(RColorBrewer)


# Set Working Directory
setwd("C:/Moein/projects/andrewju_project")

Sample_Information <- read.delim("C:/Moein/projects/andrewju_project/Sample_Information_MR.txt")

### preparing the data

# Paper Data(hg19)
paper_tc_data <- fread("Paper_TC_Data.txt", skip = 1)

# Refset determined using IGV
tc_maf <- preload_maf(maf = paper_tc_data, refset = "ces.refset.hg19", 
                      sample_col = "Sample ID",  start_col = "Position",
                      ref_col = "Reference", tumor_allele_col = "Alternate")

# Loading TCGA Data (hg38)
tcga_maf_file <- "TCGA-THCA.maf.gz"
if (!file.exists(tcga_maf_file)) {
  get_TCGA_project_MAF(project = "THCA", filename = tcga_maf_file)
}

tcga_maf <- preload_maf(maf = tcga_maf_file, 
                        chain_file = "hg38ToHg19.over.chain", 
                        refset = "ces.refset.hg19")

# Loading Genie Data ----
genie <- fread("GENIE_Mutation_Data.txt")
genie <- genie[, .(Hugo_Symbol, Chromosome, Start_Position, End_Position, Variant_Classification, 
                   Variant_Type, Reference_Allele, Tumor_Seq_Allele1, Tumor_Seq_Allele2, 
                   dbSNP_RS, Tumor_Sample_Barcode)]
genie <- genie %>% filter(`Variant_Type` == "SNP")
genie_maf <- preload_maf(maf = genie, refset = "ces.refset.hg19")

genie_maf <- genie_maf %>%
  inner_join(Sample_Information, by = "Unique_Patient_Identifier") %>%
  select(all_of(names(genie_maf)))


#Loading 2nd paper data(hg19):
paper_data_2nd <- fread("2nd_paper_Data.txt")

# Refset determined using paper:
paper_2nd_maf <- preload_maf(maf = paper_data_2nd, refset = "ces.refset.hg19", 
                             sample_col = "Unique_Patient_Identifier",  start_col = "Start", chr_col = "Chr",
                             ref_col = "Ref", tumor_allele_col = "Alt")



# Keep samples where column Problem is equal to NA:
tc_maf <- tc_maf %>% filter(is.na(problem))
tcga_maf <- tcga_maf %>% filter(is.na(problem))
genie_maf <- genie_maf %>% filter(is.na(problem))
paper_2nd_maf <- paper_2nd_maf %>% filter(is.na(problem))

# keeping only samples that do not occur at germline variant sites:
tc_maf <- tc_maf %>% filter (`germline_variant_site` == FALSE)
tcga_maf <- tcga_maf %>% filter (`germline_variant_site` == FALSE)
genie_maf <- genie_maf %>% filter (`germline_variant_site` == FALSE)
paper_2nd_maf <- paper_2nd_maf %>% filter (`germline_variant_site` == FALSE)

# keeping only samples that do not occur in repetitive regions 
tc_maf <- tc_maf %>% filter (`repetitive_region` == FALSE | cosmic_site_tier %in% 1:3)
tcga_maf <- tcga_maf %>% filter (`repetitive_region` == FALSE | cosmic_site_tier %in% 1:3)
genie_maf <- genie_maf %>% filter (`repetitive_region` == FALSE | cosmic_site_tier %in% 1:3)
paper_2nd_maf <- paper_2nd_maf %>% filter (`repetitive_region` == FALSE | cosmic_site_tier %in% 1:3)

# keeping snv:
tc_maf <- subset(tc_maf, variant_type == "snv")
tcga_maf <- subset(tcga_maf, variant_type == "snv")
genie_maf <- subset(genie_maf, variant_type == "snv")
paper_2nd_maf <- subset(paper_2nd_maf, variant_type == "snv")

# 3. CESAnalysis Creation and General Results ----
# Creating CESAnalysis
cesa <- CESAnalysis(refset = "ces.refset.hg19")

# Filter Variants
top_tgs_genes <- c("TP53", "PIK3CA", "TERT", "NF1", "NRAS", "BRAF", "CDKN2A", "CDKN2B", 
                   "NKX2-1","RET", "KMT2C", "BCOR", "TBX3", "EIF1AX", "RBM10", 
                   "ATM", "ARID1A")

tgs_coverage <- ces.refset.hg19$gr_genes[ces.refset.hg19$gr_genes$names %in% top_tgs_genes]

# Loading MAF into CESAnalysis
cesa <- load_maf(cesa = cesa, maf = tc_maf, coverage = "targeted", maf_name = "THCA",
                 covered_regions = tgs_coverage, covered_regions_name = "top_genes",
                 covered_regions_padding = 10)

cesa <- load_maf(cesa = cesa, maf = tcga_maf, coverage = "exome", maf_name = "TCGA_THCA")
cesa <- load_maf(cesa, maf = genie_maf, maf_name = "Genie_THCA", coverage = "targeted",
                 covered_regions = tgs_coverage, covered_regions_name = "top_genes",
                 covered_regions_padding = 10)

cesa <- load_maf(cesa = cesa, maf = paper_2nd_maf, coverage = "exome", maf_name = "paper_2nd")


cesa <- load_sample_data(cesa, Sample_Information)


#Defining group
primary_groups_all <- cesa$samples[Sample_type == "Primary", unique(Unique_Patient_Identifier)]
meta_groups_all <- cesa$samples[Sample_type == "Metastasis", unique(Unique_Patient_Identifier)]

#defining groups for gene mutation rate using exome:
primary_groups <- cesa$samples[Sample_type == "Primary" & coverage == "exome", unique(Unique_Patient_Identifier)]
meta_groups <- cesa$samples[Sample_type == "Metastasis" & coverage == "exome", unique(Unique_Patient_Identifier)]

cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa, covariates = "THCA", samples = primary_groups, save_all_dndscv_output = T)
cesa_samples_by_groups <- gene_mutation_rates(cesa = cesa_samples_by_groups, covariates = "THCA", samples = meta_groups, save_all_dndscv_output = T)


RefCDS = ces.refset.hg19$RefCDS
dndscv_gene_names <- cesa_samples_by_groups$gene_rates$gene
nsyn_sites = sapply(RefCDS[dndscv_gene_names], function(x) colSums(x[["L"]])[1])


# selecting mutation rate data for samples in primary_groups
samples_in_primary_groups <- length(unique(cesa_samples_by_groups$dNdScv_results$rate_grp_1$annotmuts$sampleID))

# selecting mutation rate data for samples in metastasis_groups
samples_in_metastasis_groups <- length(unique(cesa_samples_by_groups$dNdScv_results$rate_grp_2$annotmuts$sampleID))

library(tidyverse)

### creating a data frame with mutation rate data for primary_groups and metastasis_groups
mut_rate_df <- tibble(gene = cesa_samples_by_groups$dNdScv_results$rate_grp_1$genemuts$gene_name,
                      exp_primary_mu = cesa_samples_by_groups$dNdScv_results$rate_grp_1$genemuts$exp_syn_cv,
                      exp_metastasis_mu = cesa_samples_by_groups$dNdScv_results$rate_grp_2$genemuts$exp_syn_cv)

# Add n_syn_sites column to mut_rate_df:
mut_rate_df$n_syn_sites = nsyn_sites[mut_rate_df$gene]

mut_rate_df %>% 
  mutate(primary_mu = (exp_primary_mu / n_syn_sites) / samples_in_primary_groups) %>%
  mutate(metastasis_mu = (exp_metastasis_mu / n_syn_sites) / samples_in_metastasis_groups) %>%
  mutate(cancer_greater = metastasis_mu > primary_mu) -> 
  mut_rate_df


# saving gene mutation rates into separate data frame:
primary_rate <- mut_rate_df %>%
  select(gene, rate = primary_mu) %>%
  data.table::setDT()
Meta_rate <- mut_rate_df %>%
  select(gene, rate = metastasis_mu) %>%
  data.table::setDT()


# clear the gene rates in the cesa object 
cesa_samples_by_groups <- clear_gene_rates(cesa = cesa_samples_by_groups)

# setting gene rates for primary and metastasis:
cesa_samples_by_groups <- set_gene_rates(cesa = cesa_samples_by_groups, rates = primary_rate, missing_genes_take_nearest = T, samples = cesa$samples[Sample_type=="Primary"]) 
cesa_samples_by_groups <- set_gene_rates(cesa = cesa_samples_by_groups, rates = Meta_rate, missing_genes_take_nearest = T, samples = cesa$samples[Sample_type=="Metastasis"]) 

# infer trinculeotide-context-specific relative rates of SNV mutation from a mutational signature analysis
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "THCA", 
                                                            treatment_naive = TRUE)

# estimating trinucleotide mutation rates
cesa_samples_by_groups <- trinuc_mutation_rates(cesa = cesa_samples_by_groups, signature_set = "COSMIC_v3.2", 
                                                signature_exclusions = signature_exclusions)

twostage_final <- cesa_samples_by_groups
saveRDS(twostage_final, file="twostage_final.rds")

mut_rate_final <- data.frame(gene=twostage_final@mutrates$gene,
                                 primary_rate=twostage_final@mutrates$rate_grp_1,
                                 metastasis_rate=twostage_final@mutrates$rate_grp_2)

#Mann–Whitney U test:
library(rcompanion)

wilcox_test_result <- wilcox.test(mut_rate_final$primary_rate, mut_rate_final$metastasis_rate)
print(wilcox_test_result)

combined_rates <- c(mut_rate_final$primary_rate, mut_rate_final$metastasis_rate)
groups <- c(rep("pri", length(mut_rate_final$primary_rate)), rep("meta", length(mut_rate_final$metastasis_rate)))
rank_biserial_pri_meta <- wilcoxonR(x = combined_rates, g = groups)
cat("Rank-biserial correlation (pri vs meta):", rank_biserial_pri_meta, "\n")


# Subset and rename columns
mut_rate_subset <- mut_rate_final[mut_rate_final$gene %in% top_tgs_genes, 
                                  c("gene","primary_rate","metastasis_rate")]


write.csv(mut_rate_subset, "mutation_rates_primary_vs_metastasis.csv", row.names = FALSE)



#End
