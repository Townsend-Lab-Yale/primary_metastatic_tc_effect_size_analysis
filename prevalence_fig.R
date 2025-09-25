# Load necessary libraries
library(cancereffectsizeR)
library(data.table)
library(readxl)
library(dplyr)
library(ces.refset.hg19)
library(stringr)
library(ggrepel)
library(ggpubr)

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
top_tgs_genes <- c("TP53", "PIK3CA", "TERT", "NF1", "NRAS", "BRAF", "CDKN2B", 
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



# Extract prevalence data
prevalence <- cesa@maf %>%
  inner_join(Sample_Information, by = "Unique_Patient_Identifier")

# Classify data into different Sample_type 
primary <- prevalence %>% filter(Sample_type == "Primary")
metastasis <- prevalence %>% filter(Sample_type == "Metastasis")

#Count the number of samples in each group:
primary_samples <- primary %>%
  summarise(Unique_Count = n_distinct(Unique_Patient_Identifier))

metastasis_samples <- metastasis %>%
  summarise(Unique_Count = n_distinct(Unique_Patient_Identifier))


# Define genes of interest with correct order (Removed "MUC16")
genes_of_interest <- c("TP53", "PIK3CA", "TERT", "NF1", "NRAS", "BRAF", "CDKN2B", 
                   "NKX2-1","RET", "KMT2C", "BCOR", "TBX3", "EIF1AX", "RBM10", 
                   "ATM", "ARID1A")


# Filter data for genes of interest
filtered_primary <- primary %>%
  filter(genes %in% genes_of_interest) %>%
  distinct(Unique_Patient_Identifier, genes)

filtered_metastasis <- metastasis %>%
  filter(genes %in% genes_of_interest) %>%
  distinct(Unique_Patient_Identifier, genes)

# Count unique patients per gene and normalize frequency
primary_gene_frequencies <- filtered_primary %>%
  group_by(genes) %>%
  summarise(Unique_Patient_Count = n(), .groups = "drop") %>%
  mutate(Frequency_Percentage = (Unique_Patient_Count / 1601) * 100, Category = "Primary")

metastasis_gene_frequencies <- filtered_metastasis %>%
  group_by(genes) %>%
  summarise(Unique_Patient_Count = n(), .groups = "drop") %>%
  mutate(Frequency_Percentage = (Unique_Patient_Count / 984) * 100, Category = "Metastasis")

# Combine datasets
combined_gene_frequencies <- bind_rows(primary_gene_frequencies,
                                       metastasis_gene_frequencies)

# Convert Category to a factor for correct grouping
combined_gene_frequencies$Category <- factor(combined_gene_frequencies$Category,
                                             levels = c("Primary", "Metastasis"))
# Calculate total frequency per gene for ordering
gene_order <- combined_gene_frequencies %>%
  group_by(genes) %>%
  summarise(Total_Frequency = sum(Frequency_Percentage)) %>%
  arrange(desc(Total_Frequency)) %>%
  pull(genes)

# Update the genes factor levels based on total frequency
combined_gene_frequencies$genes <- factor(combined_gene_frequencies$genes, levels = gene_order)


# Define new distinct color scheme for each tumor stage
stage_colors <- c("Primary" = "#1b9e77",  # Green
                  "Metastasis" = "#7570b3")   # Blue

# Generate one-panel grouped bar chart
ggplot(combined_gene_frequencies, aes(x = genes, y = Frequency_Percentage, fill = Category)) +
  geom_bar(stat = "identity", color = "black", width = 0.7, position = position_dodge(width = 0.8)) +  
  labs(x = "Genes", y = "Mutation Frequency (%)") +  
  scale_fill_manual(values = stage_colors) +  # Apply distinct colors per stage
  theme_pubr(base_size = 16) +  # Professional journal-quality theme
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16, color = "black"),  
    axis.text.y = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 20, face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 16, face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("Gene_Mutation_Frequency_Grouped.png", width = 12, height = 6, dpi = 600)

#End