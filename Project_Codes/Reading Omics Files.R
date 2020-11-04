library(pacman)
p_load(tidyverse, openxlsx, piggyback)
# devtools::install_github("ropensci/piggyback@87f71e8")
piggyback::pb_track("Project_Datasets/*")
sample_info <- read.csv("./Project_Datasets/sample_info.csv", stringsAsFactors = FALSE)
#Raw from https://depmap.org/portal/download/ downloaded 2/11/2020


#Raw from https://portals.broadinstitute.org/ccle/data CCLE_metabolomics_20190502.csv 2/11/2020
Metabolites <- openxlsx::read.xlsx("./Project_Datasets/CCLE_metabolites_landscape_of_cancer.xlsx",
                         sheet = "1-clean data")[-760,] %>% 
    mutate(X1 = str_match(X1,"^([:graph:]*?)_")[,2]) %>%
    column_to_rownames("X1")%>% t()


#Raw from https://gygi.med.harvard.edu/publications/ccle 2/11/2020
CCLE_proteins <- openxlsx::read.xlsx("./Project_datasets/Table_S2_Protein_Quant_Normalized.xlsx", sheet =2) %>%
   .[,!str_detect(colnames(.),"Peptides|Column")]
colnames(CCLE_proteins) <- str_remove_all(colnames(CCLE_proteins),"_TenPx..")
rownames(CCLE_proteins) <- NULL
CCLE_proteins <- CCLE_proteins %>%
    column_to_rownames(var = "Uniprot_Acc")
CCLE_proteins <- CCLE_proteins[,-c(1:5)]
colnames(CCLE_proteins) <- str_match(colnames(CCLE_proteins),"^([:graph:]*?)_")[,2]
CCLE_proteins <- CCLE_proteins[,which(!duplicated(colnames(CCLE_proteins)))]

#Raw from https://depmap.org/portal/download/ 2/11/2020
RNA_seq <- read.csv("./Project_Datasets/CCLE_expression.csv") 
RNA_seq <- inner_join(sample_info[,1:2],RNA_seq, by = c("DepMap_ID" = "X"))[,-1] %>%
    column_to_rownames(var = "stripped_cell_line_name") %>% t() %>%
    magrittr::set_rownames(str_match(rownames(.), "^([:graph:]*?)\\.")[,2])

#Raw from https://depmap.org/portal/download/ 2/11/2020
Achilles <- read.csv("./Project_Datasets/Achilles_gene_effect.csv")
Achilles <- inner_join(sample_info[,1:2],Achilles, by = "DepMap_ID")[,-1] %>%
    column_to_rownames(var = "stripped_cell_line_name") %>% t() %>%
    magrittr::set_rownames(str_match(rownames(.), "^([:graph:]*?)\\.")[,2])

#Raw from https://depmap.org/portal/download/ 2/11/2020
RNAi <- read.csv("./Project_Datasets/D2_combined_gene_dep_scores.csv")
RNAi <- RNAi %>% mutate(Gene_name = str_match(Gene_name,"^([:graph:]*?) ")[,2]) %>%
    column_to_rownames(var = "Gene_name") %>% 
    setNames(str_match(colnames(.), "^([:graph:]*?)_")[,2])


#Raw from  https://www.nature.com/articles/s41467-019-09695-9#Sec28 2/11/2020
NCI_60_metabolites <- read.xlsx("./Project_Datasets/41467_2019_9695_MOESM2_ESM.xlsx", sheet = 3, startRow = 4, na.strings = "NaN") %>%
    .[str_detect(.$`Annotation.ID`,"H|C"),-c(1:3,5)] %>% setNames(str_remove_all(colnames(.), "^[:graph:]*?_"))  %>%
    na.omit() %>% remove_rownames() %>% column_to_rownames("Annotation.ID")

#Raw from https://www.sciencedirect.com/science/article/pii/S2589004219304407#mmc2 3/11/2020 SWATH paper 2019
NCI_60_proteins <- read.xlsx("./Project_Datasets/1-s2.0-S2589004219304407-mmc2.xlsx", sheet =6) %>%
    .[,-c(2:8)] %>% remove_rownames() %>%
    column_to_rownames(var = "protein.accession.number") %>%
    setNames(str_remove(colnames(.), "^[:graph:]*?_"))

#Raw from https://discover.nci.nih.gov/cellminer/loadDownload.do 3/11/2020 - RNAseq Composite expression
NCI_60_RNA <- readxl::read_xls("./Project_Datasets/RNA__RNA_seq_composite_expression.xls", skip = 10) %>%
    .[,-c(2:6)] %>% setNames(str_remove_all(colnames(.), "^[:graph:]*:|-| ")) %>%
    column_to_rownames("Genenamed")

# Mutations <- read.xlsx("./Project_Datasets/CCLE_mutations.xlsx", sheet = 2)[,c("Hugo_Symbol", "DepMap_ID")] %>%
#     left_join(sample_info[,1:2]) %>%
#     mutate(Value  = 1)
# 
# Mutations<- Mutations[,c(1,3:4)] %>%
#     group_by(Hugo_Symbol, stripped_cell_line_name) %>%
#     dplyr::summarise(abundance = sum(Value)) %>%
#     tidyr::pivot_wider(names_from = stripped_cell_line_name, values_from = abundance, 
#                        values_fill = 0) %>%
#     column_to_rownames("Hugo_Symbol")

#Raw from https://github.com/J-Andy/Protein-expression-in-human-cancer/tree/master/data 3/11/2020
Vizcaino_Proteome <- read.delim("./Project_Datasets/proteinGroups_cellLines_ppbNorm_min50validValues_svdImputeMissing_removeBatchEffect.txt")

save(RNA_seq, 
     RNAi, 
     CCLE_proteins, 
     sample_info,
     Metabolites,
     Achilles,
     file = "./Project_Datasets/CCLE_OMICS.rda")
     
save(NCI_60_proteins,
     NCI_60_metabolites,
     NCI_60_RNA,
     Vizcaino_Proteome,
     file = "./Project_Datasets/NCI60_OMICS.rda")

pb_new_release("Skourtis/multiOmics", "v1")
pb_track() %>% pb_upload(repo = "Skourtis/multiOmics")

