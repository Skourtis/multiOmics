pacman::p_load(org.Hs.eg.db)
x <- org.Hs.egUNIPROT
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
Entrez = tibble(Entrez = names(xx),
                    Uniprot = xx)
Entrez <- unnest(Entrez, cols = c(Uniprot))

pacman::p_load(graphite)
humankegg <- pathways("hsapiens", "reactome")

humankegg_df <- data.frame(SubPathway = NULL,
                           ID = NULL)
for(i in names(humankegg)){
    print(i)
    if(length(nodes(humankegg[[i]]))>0){
        humankegg_df <- rbind(humankegg_df,data.frame(SubPathway = i,
                                  ID = nodes(humankegg[[i]])))
    }
    
}
humankegg_df <- humankegg_df %>% subset(!duplicated(ID)) %>%
    mutate(ID = str_remove_all(ID,"ENTREZID:")) %>% 
    left_join(Entrez, by = c("ID" = "Entrez")) %>% na.omit() %>%
    dplyr::select(-ID)





HUMAN_9606_idmapping <- read_tsv("./Project_Datasets/HUMAN_9606_idmapping.dat",
                                 col_names = FALSE) %>% 
    setNames(c("Uniprot", "Type", "ID"))

HUMAN_9606_idmapping_hsa <- HUMAN_9606_idmapping[str_detect(HUMAN_9606_idmapping$ID, "hsa:"),-2] %>% 
    left_join(HUMAN_9606_idmapping[HUMAN_9606_idmapping$Type == "Gene_Name", -2], by = "Uniprot") %>% 
    na.omit() %>% .[,-1] %>% setNames(c("KEGG_ID", "Gene.names"))


#Raw from

Master_pathways_metabolites <- openxlsx::read.xlsx("./Project_Datasets/JCI71180sd2.xlsx", 
                                                   sheet = 2, startRow = 5, 
                                                   cols = c(2,6,10)) %>%
    pivot_longer(-SUPER_PATHWAY, names_to = "Type", values_to = "ID") %>%
    setNames(c("Pathway","Type","ID")) %>%
    mutate(Type = "Metabolite") %>%
    subset(!is.na(ID)) %>% unique()

#Raw from https://www.sciencedirect.com/science/article/pii/S2211124718304388?via%3Dihub#mmc2
Master_pathways <- map_df(1:7, ~openxlsx::read.xlsx("./Project_Datasets/1-s2.0-S2211124718304388-mmc2.xlsx", 
                                          startRow = 2, 
                                          sheet = .x)) %>% 
    left_join(HUMAN_9606_idmapping[HUMAN_9606_idmapping$Type == "Gene_Name",-2], by = c("Genes" = "ID")) %>%
    pivot_longer(!Pathway, names_to = "Type", values_to="ID") %>% na.omit() %>% unique() %>%
    rbind(Master_pathways_metabolites)

genes_reactions <- data.frame(Gene_id = NULL,
                              Reaction = NULL, 
                              pathway = NULL,
                              stringsAsFactors = F)
# compound_reactions <- data.frame(Reaction = NULL,
#                                  Substrate = NULL,
#                                  Product = NULL,
#                                  Direction = NULL,
#                                  pathway = NULL,
#                                  stringsAsFactors = F)


Retrieve_Genes_KEGG <- function(pathway, Gene_df =  genes_reactions){
    # pathway = pathways[1]
    # Gene_df = genes_reactions
    # Compound_df = compound_reactions
    KEGGgraph::retrieveKGML(pathway, organism="hsa", destfile=tmp, method="auto", quiet=TRUE)
    KGMLGraph <- KEGGgraph::parseKGML2Graph(tmp)
    KGML_n <- KEGGgraph::parseKGML(tmp)
    # mapkG2 <- KEGGpathway2Graph(KGML_n, expandGenes=TRUE)
    # set.seed(124)
    # randomNodes <- sample(nodes(mapkG2), 25)
    # mapkGsub <- subGraph(randomNodes, mapkG2)
    # plot(mapkGsub)
    
    for (i in 1:length(KGMLGraph@nodeData@defaults[["KEGGNode"]][["nodes"]])){
        df <- data.frame(Gene_id = KGMLGraph@nodeData@defaults[["KEGGNode"]][["nodes"]][[i]]@entryID,
                         Reaction = KGMLGraph@nodeData@defaults[["KEGGNode"]][["nodes"]][[i]]@reaction,
                         pathway = pathway,
                         stringsAsFactors = FALSE)
        Gene_df <- rbind(Gene_df,df)
        
    }
    
    
    # for (i in 1:length(KGML_n@reactions)){
    #     df <- data.frame(Reaction = KGML_n@reactions[[i]]@name,
    #                      Substrate = KGML_n@reactions[[i]]@substrateName,
    #                      Product = KGML_n@reactions[[i]]@productName,
    #                      Direction = KGML_n@reactions[[i]]@type,
    #                      pathway = pathway,
    #                      stringsAsFactors = FALSE)
    #     Compound_df <- rbind(Compound_df,df)
    #     
    # }
    
    Gene_df
}

tmp <- tempfile()
pathways <- c(#"01110", # Biosynthesis of secondary metabolites
    #"01120", # Microbial metabolism in diverse environments
    "01200", # Carbon metabolism
    "01210", # 2-Oxocarboxylic acid metabolism
    "01212", # Fatty acid metabolism
    "01230") # Biosynthesis of amino acids
#"01220") # Degradation of aromatic compounds
#KEGGgraph::retrieveKGML(pathways, organism="hsa", destfile=tmp, method="auto", quiet=TRUE)
pathways <- c("03010",# Ribosome [PATH:hsa03010]
              "00970", #Aminoacyl-tRNA biosynthesis [PATH:hsa00970]
              "03013", #RNA transport [PATH:hsa03013]
              "03015", #mRNA surveillance pathway [PATH:hsa03015]
              "03008",# Ribosome biogenesis in eukaryotes [PATH:hsa03008]"
              "03060", # Protein export [PATH:hsa03060]
              "04141", # Protein processing in endoplasmic reticulum [PATH:hsa04141]
              "04130", # SNARE interactions in vesicular transport [PATH:hsa04130]
              "04120", # Ubiquitin mediated proteolysis [PATH:hsa04120]
              "04122", # Sulfur relay system [PATH:hsa04122]
              "03050", # Proteasome [PATH:hsa03050]
              "03018") # RNA degradation [PATH:hsa03018])




Gene_pathways <- map_df(pathways,Retrieve_Genes_KEGG) %>% 
    left_join(HUMAN_9606_idmapping_hsa, by= c("Gene_id" = "KEGG_ID")) %>%
    subset(!duplicated(Gene.names)) %>% 
    left_join(HUMAN_9606_idmapping[HUMAN_9606_idmapping$Type == "Gene_Name",-2], 
              by = c("Gene.names" = "ID"))

Gene_pathways <- Gene_pathways %>%
    mutate(SubPathway = case_when(pathway %in% pathways[1:5] ~ "Ribosomal",
                                  pathway %in% pathways[-c(1:5)] ~ "Proteosome"))

humankegg_df <- rbind(humankegg_df, Gene_pathways[,c("SubPathway","Uniprot")])

#Raw from http://tko.ccbr.utoronto.ca/ 4/11/2020
Essential_genes <- openxlsx::read.xlsx("./Project_Datasets/reference_essentials_and_nonessentials_sym_hgnc_entrez.xlsx")$Gene
SMPDB_pathways <- read.csv("./Project_Datasets/smpdb_pp_metabolism.proteins.csv", header = F) %>%
    setNames(c("ID","Pathway","Type", "Uniprot","Name","HMDB","Other","mRNA","Gene_Name","Location"))

save(Master_pathways,
     Essential_genes,
     SMPDB_pathways,
     Gene_pathways,
     humankegg_df,
     file = "./Project_Datasets/Metabolic_Pathway_annotations.RData")
