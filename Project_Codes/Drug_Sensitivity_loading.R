### Preparing Drug sensitivity ###
pacman::p_load("stringr", "readr","tidyverse","renv")
renv::init()
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.12")
renv::snapshot()
###downloaded from https://dgidb.org/downloads 12/12/2020
DGIdb <- read_tsv(here::here("Project_datasets","interactions.tsv"))

Metabolic_drugs <- openxlsx::read.xlsx(here::here("Project_Datasets",
                                                  "metabolic compounds.xlsx"),
                                       sheet =2)

CCLE_Drug_screening_info <- readr::read_csv(here::here("Project_Datasets",
                                                       "primary-screen-replicate-treatment-info.csv"))

CCLE_metabolic_drug <- intersect(Metabolic_drugs$Product.Name,str_to_title(CCLE_Drug_screening_info$name))
Purine_Drug <- CCLE_Drug_screening_info$name[str_detect(CCLE_Drug_screening_info$moa,"purine")] %>% 
    na.omit() %>% as.vector()
CCLE_Drug_screening_original <- read_csv(here::here("Project_Datasets",
                                                    "primary-screen-replicate-collapsed-logfold-change.csv"))
drug_targets <- CCLE_Drug_screening_info %>% 
    dplyr::select(broad_id,target) %>% distinct()%>%
    na.omit() %>% separate_rows(target, sep = ", ") %>% 
    #left_join(Master_pathways[,-2], by = c("target" = "ID")) %>%
    left_join(humankegg_df[,c("Gene_name","Reactome_Pathway")], by = c("target" = "Gene_name")) %>%
    mutate(Reactome_Pathway = if_else(is.na(Reactome_Pathway),"Not_Found_reactome",Reactome_Pathway),
           broad_id = str_match(broad_id,"^BRD-([:graph:]*?)-")[,2])
drug_targets_2 <- CCLE_Drug_screening_info %>% 
    dplyr::select(broad_id,target) %>% distinct()%>%
    na.omit() %>% separate_rows(target, sep = ", ") %>% 
    left_join(Master_pathways[,-2], by = c("target" = "ID")) %>%
    #left_join(SMPDB_pathways[,c("Gene_Name","Pathway")], by = c("target" = "Gene_Name")) %>%
    mutate(Pathway = if_else(is.na(Pathway),"Non_Metabolic",Pathway),
           broad_id = str_match(broad_id,"^BRD-([:graph:]*?)-")[,2])


Detecting_dominant_target <- function(x,y,z){
    Pathway <- subset(y, broad_id == x) %>%
        pull({{z}}) %>% table() %>% names() %>% .[1]
    data.frame(Drug = x) %>%
        mutate({{z}} := Pathway)
}

drug_targets <- map_dfr(unique(drug_targets$broad_id),
                        ~Detecting_dominant_target(x = .x, y = drug_targets,z = "Reactome_Pathway")) %>%
    left_join(map_dfr(unique(drug_targets_2$broad_id),
                      ~Detecting_dominant_target(x = .x, y = drug_targets_2, z = "Pathway")), by = "Drug")

Sanger_Drug_Screening_original <- read_csv(here::here("Project_Datasets",
                                                      "sanger-viability.csv"))
sample_info <- read_csv(here::here("Project_Datasets","sample_info.csv"))
Sanger_Drug_Screening <- inner_join(Sanger_Drug_Screening_original[,c(7,10,11,12)], 
                                    sample_info[,1:2],
                                    by = c("ARXSPAN_ID" = "DepMap_ID"))[,-2] 
#Sanger_metabolic_drug <- intersect(toupper(Metabolic_drugs$Product.Name),toupper(Sanger_Drug_Screening$DRUG_NAME))
#Sanger_Metabolic <- CCLE_Drug_screening_info[str_to_title(CCLE_Drug_screening_info$name) %in% str_to_title(Metabolic_drugs$Product.Name),c("broad_id","name","moa")] %>% unique()
#Sanger_Metabolic$broad_id <-  str_match(CCLE_Metabolic$broad_id,"^BRD-([:graph:]*?)-")[,2]
#Sanger_Drug_Screening <- pivot_longer(Sanger_Drug_Screening, -stripped_cell_line_name, names_to = "broad_id", values_to = "Sensitivity")
#Sanger_Drug_Screening$BROAD_ID <- str_match(Sanger_Drug_Screening$BROAD_ID,"^BRD-([:graph:]*)")[,2]
#Sanger_Drug_Screening <- Sanger_Drug_Screening %>% inner_join(unique(CCLE_Metabolic))
#Sanger_Drug_Screening <- Sanger_Drug_Screening[Sanger_Drug_Screening$DRUG_NAME %in% Sanger_metabolic_drug,] #%>% right_join(CCLE_Metabolic)



#CCLE_Drug_screening_original <- read.csv("./Project_Datasets/secondary-screen-replicate-collapsed-logfold-change.csv")
CCLE_Drug_screening <- inner_join(CCLE_Drug_screening_original,
                                  sample_info[,1:2],
                                  by = c("X1" = "DepMap_ID"))[,-1]

#colnames(CCLE_Drug_screening_test) <- str_match(colnames(CCLE_Drug_screening_test),"^BRD\\.([:graph:]*?)\\.")[,2]
CCLE_Metabolic <- CCLE_Drug_screening_info[str_to_title(CCLE_Drug_screening_info$name) %in% str_to_title(Metabolic_drugs$Product.Name),c("broad_id","name","moa")] %>% unique()
CCLE_Metabolic$broad_id <-  str_match(CCLE_Metabolic$broad_id,"^BRD-([:graph:]*?)-")[,2]

CCLE_Drug_screening <- pivot_longer(CCLE_Drug_screening, -stripped_cell_line_name, names_to = "broad_id", values_to = "Sensitivity")
CCLE_Drug_screening$broad_id <- str_match(CCLE_Drug_screening$broad_id,"^BRD-([:graph:]*?)-")[,2]
#CCLE_Drug_screening <- CCLE_Drug_screening %>% inner_join(unique(CCLE_Metabolic))
#CCLE_Drug_screening <- left_join(groups[,1:2], CCLE_Drug_screening, by = c("Cell_line" = "stripped_cell_line_name"))
#CCLE_Drug_screening <- CCLE_Drug_screening%>% right_join(CCLE_Metabolic)
CCLE_Drug_screening_1 <- CCLE_Drug_screening %>% 
    group_by(broad_id) %>% 
    #subset(broad_id %in% c("A00077618","A00100033","A00147595")) %>% 
    mutate(value_z = scale(Sensitivity)) %>%
    mutate(Classification = case_when(
        value_z<(-0.68) ~ "Sensitive",
        value_z>(0.68) ~ "Resistant",
        between(value_z,-0.68,0.68)~ "Vague"))


#CCLE_Drug_screening_metab <- CCLE_Drug_screening[,colnames(CCLE_Drug_screening) %in% c("X",as.character(CCLE_Metabolic))]

#CCLE_Metabolic_drugs <- Metabolic_drugs[Metabolic_drugs$Product.Name %in% CCLE_metabolic_drug,]
Sanger_Drug_Screening$DRUG_NAME <- str_replace_all(Sanger_Drug_Screening$DRUG_NAME,'-..-',"-")
Sanger_Drug_Screening$DRUG_NAME <- str_replace_all(Sanger_Drug_Screening$DRUG_NAME,'/',"-")
Sanger_Drug_Screening <- Sanger_Drug_Screening %>% 
    mutate(BROAD_ID = str_match(BROAD_ID,"^BRD-([:graph:]*?)$")[,2],
           Classification = case_when(
               Z_SCORE_PUBLISHED<(-0.68) ~ "Sensitive",
               Z_SCORE_PUBLISHED>(0.68) ~ "Resistant",
               between(Z_SCORE_PUBLISHED,-0.68,0.68)~ "Vague")) %>%
    left_join(drug_targets, by = c("BROAD_ID"= "Drug"))
#Sanger_Drug_Screening <- Sanger_Drug_Screening$DRUG_NAME %>% na.omit()


CancerDAP <- read_tsv(here::here("Project_Datasets","SubpathwayEntry.txt"), col_names = F) %>%
    setNames(c("Drug","Status","Pathway","Pathway_ID","N_gene","Score","p_val","Gene","Gene_ID")) %>% 
    select(-c(Status,Pathway_ID,Gene_ID)) %>%
    mutate(Drug = str_remove(Drug," rep."))

save(CCLE_Drug_screening,Sanger_Drug_Screening,CancerDAP,CCLE_Drug_screening_1, file = here::here("Project_Datasets","Drug_sensitivity.RData"))
