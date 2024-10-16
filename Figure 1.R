library(cancereffectsizeR)
library(data.table)
library(readxl)
library(dplyr)
library(ces.refset.hg19)
library(stringr)
library(ggrepel)

# Set Working Directory
setwd("C:/Moein/projects/andrewju_project")

Sample_Information <- read.delim("C:/Moein/projects/andrewju_project/Sample_Information_MR.txt")

###Checking for duplicates:
# Keep only rows with "Primary" and "Metastasis" in the "Sample_type" column
filtered_samples <- Sample_Information %>% 
  filter(Sample_type %in% c("Primary", "Metastasis"))

# Remove duplicates based on "Unique_Patient_Identifier"
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

# View the results
print("Consistent samples:")
print(consistent_samples)

print("Inconsistent samples:")
print(inconsistent_samples)



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


#Loading 2nd paper data(hg19):
paper_data_2nd <- fread("2nd_paper_Data.txt")

# Refset determined using IGV
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
top_tgs_genes <- c("TP53", "PIK3CA", "TERT", "NF1", "NF2", "NRAS", "BRAF", "CDKN2A", "CDKN2B", 
                   "NKX2-1","RET", "KMT2C", "KMT2D", "BCOR", "TBX3", "PTEN", "EIF1AX", "RBM10", 
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



# We'll use all suggested exclusions (TCGA primary tumors are treatment-naive)
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "THCA",
                                                            treatment_naive = TRUE)

# Adding information about snv_counts, raw_attributions, biological_weights and trinuc_rates
# to the CESAnalysis
cesa <- trinuc_mutation_rates(cesa,
                              signature_set = ces.refset.hg19$signatures$COSMIC_v3.2,
                              signature_exclusions = signature_exclusions)

##Figure_1:
snv_counts <- cesa$mutational_signatures$snv_counts

summed_snv_by_group <- data.table()
receptor_groups <- unique(na.omit(cesa$samples$Sample_type))
samples_with_snvs <- cesa$samples[colnames(snv_counts), on = "Unique_Patient_Identifier"]
for (grp in receptor_groups) {
  curr_samples <- samples_with_snvs[grp, Unique_Patient_Identifier, on = "Sample_type"]
  curr_snv_sum <- rowSums(snv_counts[, curr_samples])
  summed_snv_by_group[, (grp) := curr_snv_sum]
}
summed_snv_by_group <- as.matrix(summed_snv_by_group)
colnames(summed_snv_by_group)[c(1, 2)] <- c("Metastases", "Primary")
summed_snv_by_group <- summed_snv_by_group[, c("Primary", "Metastases")]
rownames(summed_snv_by_group) <- rownames(snv_counts)
Figure_1 <- MutationalPatterns::plot_96_profile(summed_snv_by_group, ymax = 0.6)
ggsave("Figure_1.png", width = 8, height = 6, dpi = 600)


**********************************************************************************************************************************************

###old code:
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

# Load more Metastasis Data ----
metastasis_1 <- read_excel("Supp_Table3.xlsx", skip = 28)
metastasis_1$Chromosome <- str_sub(metastasis_1$Chr, 4)
meta1_maf <- preload_maf(maf = metastasis_1, refset = "ces.refset.hg19", 
                         sample_col = "#CaseID",  start_col = "Start",
                         ref_col = "Ref", tumor_allele_col = "Alt",
                         keep_extra_columns =  TRUE)

# Keep samples where column Problem is equal to NA:
tc_maf <- tc_maf %>% filter(is.na(problem))
tcga_maf <- tcga_maf %>% filter(is.na(problem))
genie_maf <- genie_maf %>% filter(is.na(problem))
meta1_maf <- meta1_maf %>% filter(is.na(problem))

# keeping only samples that do not occur at germline variant sites:
tc_maf <- tc_maf %>% filter (`germline_variant_site` == FALSE)
tcga_maf <- tcga_maf %>% filter (`germline_variant_site` == FALSE)
genie_maf <- genie_maf %>% filter (`germline_variant_site` == FALSE)
meta1_maf <- meta1_maf %>% filter (`germline_variant_site` == FALSE)

# keeping only samples that do not occur in repetitive regions 
tc_maf <- tc_maf %>% filter (`repetitive_region` == FALSE | cosmic_site_tier %in% 1:3)
tcga_maf <- tcga_maf %>% filter (`repetitive_region` == FALSE | cosmic_site_tier %in% 1:3)
genie_maf <- genie_maf %>% filter (`repetitive_region` == FALSE | cosmic_site_tier %in% 1:3)
meta1_maf <- meta1_maf %>% filter (`repetitive_region` == FALSE | cosmic_site_tier %in% 1:3)

# keeping snv:
tc_maf <- subset(tc_maf, variant_type == "snv")
tcga_maf <- subset(tcga_maf, variant_type == "snv")
genie_maf <- subset(genie_maf, variant_type == "snv")
meta1_maf <- subset(meta1_maf, variant_type == "snv")

# 3. CESAnalysis Creation and General Results ----
# Creating CESAnalysis
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
cesa <- load_maf(cesa, maf = meta1_maf, coverage = "exome", maf_name = "Meta1")

# Loading Clinical Data
cesa <- load_sample_data(cesa, tcga_clinical)
cesa <- load_sample_data(cesa, genie_clinical)
cesa <- load_sample_data(cesa, Sample_Information)

# We'll use all suggested exclusions (TCGA primary tumors are treatment-naive)
signature_exclusions <- suggest_cosmic_signature_exclusions(cancer_type = "THCA",
                                                            treatment_naive = TRUE)

# Adding information about snv_counts, raw_attributions, biological_weights and trinuc_rates
# to the CESAnalysis
cesa <- trinuc_mutation_rates(cesa,
                              signature_set = ces.refset.hg19$signatures$COSMIC_v3.2,
                              signature_exclusions = signature_exclusions)

##Figure:

snv_counts <- cesa$mutational_signatures$snv_counts

summed_snv_by_group <- data.table()
receptor_groups <- unique(na.omit(cesa$samples$Sample_type))
samples_with_snvs <- cesa$samples[colnames(snv_counts), on = "Unique_Patient_Identifier"]
for (grp in receptor_groups) {
  curr_samples <- samples_with_snvs[grp, Unique_Patient_Identifier, on = "Sample_type"]
  curr_snv_sum <- rowSums(snv_counts[, curr_samples])
  summed_snv_by_group[, (grp) := curr_snv_sum]
}
summed_snv_by_group <- as.matrix(summed_snv_by_group)
colnames(summed_snv_by_group)[c(1, 2)] <- c("Metastases", "Primary")
summed_snv_by_group <- summed_snv_by_group[, c("Primary", "Metastases")]
rownames(summed_snv_by_group) <- rownames(snv_counts)
Figure_1 <- MutationalPatterns::plot_96_profile(summed_snv_by_group, ymax = 0.4)
ggsave("Figure.png", width = 8, height = 6, dpi = 600)

#End
