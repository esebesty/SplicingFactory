library(tidyverse)
library(TCGAbiolinks)
library(SplicingFactory)

#downloaded files will come here
#set for your own preference
location <- "/disk/work/users/tp1/projects/transcriptome-noise-in-cc/data/TCGAbiolinks"
#depth of your location path (depends on number of sub-directories)
pathDepth = 16

#location of downloaded kgXref table relative to 'location'
#downloaded from here: http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/kgXrefOld5.txt.gz
annLocation <- "/annotations/kgXrefOld5.txt"

#Read in annotation
annotation <- read_delim(paste0(location, annLocation),
                         "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE) %>%
  rename(gene_name = X5) %>%
  # Manual curation due to biological phenomena (1 transcript, 2 gene)
  mutate(isoform_id = ifelse(!str_detect(X1, "uc010nxr"), 
                             str_split_fixed(X1, "\\.", 2)[, 1], X1))

# TCGA data gathering
query <- GDCquery(project = "TCGA-BRCA",
                  legacy = TRUE,
                  data.category = "Gene expression",
                  data.type = "Isoform expression quantification",
                  experimental.strategy = "RNA-Seq",
                  sample.type = c("Primary Tumor","Solid Tissue Normal"))

#Download data files
#GDCdownload(query, directory = paste0(location, "/data"))

#select breast cancer samples
dataSubt <- TCGAquery_subtype(tumor = "brca") %>%
  select(patient, BRCA_Subtype_PAM50)
#select Luminal A subtype patients
LumAPatients <- dataSubt %>%
  filter(BRCA_Subtype_PAM50 == "LumA") %>%
  select(patient)

#First 20 patients with both pathological and control sample
#We keep only 40 samples to make run time short
result_table <- query[[1]][[1]] %>%
  mutate(tags = as.character(tags)) %>%
  filter(str_detect(tags, "unnormalized")) %>%
  mutate(patient = str_extract(cases, "^.{12}")) %>%
  select(patient, cases, sample_type, file_name) %>%
  group_by(patient) %>%
  filter(length(patient) != 1) %>%
  filter(patient %in% LumAPatients$patient) %>%
  ungroup() %>%
  arrange(patient) %>% slice_head(n=40)

# Data read-in
samples <- as.data.frame(list.files(paste0(location,
                                           "/data/TCGA-BRCA/legacy/Gene_expression/Isoform_expression_quantification"),
                                    recursive = TRUE, full.names = TRUE))

colnames(samples) <- "file"
#create a vector of samples tissue type (Tumor/Normal)
samples <- samples %>%
  mutate(file_name = str_split_fixed(file, "\\/", pathDepth)[, pathDepth]) %>%
  filter(file_name %in% result_table$file_name)%>%
  left_join(result_table, by = "file_name") %>%
  arrange(patient)

#Function takes vector of strings with path to data files as arguments
#Gives back list of tibbles generated from the data files with 3 columns:
#isoform ID; Raw count and scaled estimate - we do not use the latter
sampleRead <- function(x) {
  xd <- read_delim(x, "\t", escape_double = FALSE,
                   trim_ws = TRUE)
  return(xd)
}


samplesS <- sapply(as.character(samples$file), FUN = sampleRead, simplify = FALSE,
                   USE.NAMES = TRUE)

# Preprocess
data <- bind_rows(samplesS, .id = "id") %>%
  mutate(id = str_split_fixed(id, "\\/", pathDepth)[, pathDepth],
         isoform_id = ifelse(!str_detect(isoform_id, "uc010nxr"), 
                             str_split_fixed(isoform_id, "\\.", 2)[, 1], isoform_id)) %>%
  left_join(select(result_table, file_name, patient, sample_type), by = c("id" = "file_name"))%>%
  arrange(isoform_id) %>%
  mutate(patient = ifelse(sample_type == "Solid Tissue Normal", paste0(patient, "_N"), paste0(patient, "_T")))


table <- select(data, isoform_id, raw_count, scaled_estimate, patient) %>%
  pivot_wider(id_cols = isoform_id, names_from = patient, values_from = c(raw_count, scaled_estimate)) %>%
  left_join(select(annotation, isoform_id, gene_name), by = "isoform_id")

genes <- table$gene_name

samples <- select(data, patient, sample_type) %>%
  distinct()
samples <- as.character(samples$sample_type)

#create table with readcounts; each column represents a sample
count_table <- select(data, isoform_id, raw_count, patient) %>%
  filter(isoform_id %in% table$isoform_id) %>% 
  pivot_wider(id_cols = isoform_id, names_from = patient, values_from = raw_count) %>%
  select(-isoform_id) %>% 
  as.data.frame()

# Entropy calculation
Laplace_diversity <- calculate_diversity(count_table, genes, method = "laplace")
Laplace_readcount_Wilcox <- calculate_difference(Laplace_diversity, samples, control = "Solid Tissue Normal",
                                                 method = "mean", test = "wilcoxon")%>%
  arrange(desc(abs(log2_fold_change)))


#add gene names to the rows/observations
count_table <- cbind(genes, count_table)

#select top 100 most diverse genes
#to make sure that sample dataset will work
top_genes <- Laplace_readcount_Wilcox[1:100,]

#add further 200 genes to the sample dataset
set.seed(42)
random_genes <- Laplace_readcount_Wilcox[sample(nrow(Laplace_readcount_Wilcox[100:nrow(Laplace_readcount_Wilcox),]), 200), ]
geneset <- rbind(top_genes, random_genes) %>% select(genes) 
tcga_brca_luma_dataset <- count_table %>% filter(genes %in% geneset$genes)

save(tcga_brca_luma_dataset, file="tcga_brca_luma_dataset.RData")
