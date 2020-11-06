HUMAN_9606_idmapping <- read_tsv("./Project_Datasets/HUMAN_9606_idmapping.dat",
                                 col_names = FALSE) %>% 
    setNames(c("Uniprot", "Type", "ID"))

HUMAN_9606_idmapping_hsa <- HUMAN_9606_idmapping[str_detect(HUMAN_9606_idmapping$ID, "hsa:"),-2] %>% 
    left_join(HUMAN_9606_idmapping[HUMAN_9606_idmapping$Type == "Gene_Name", -2], by = "Uniprot") %>% 
    na.omit() %>% .[,-1] %>% setNames(c("KEGG_ID", "Gene.names"))

#Raw from https://www.sciencedirect.com/science/article/pii/S2211124718304388?via%3Dihub#mmc2
Master_pathways <- map_df(1:7, ~read.xlsx("./Project_Datasets/1-s2.0-S2211124718304388-mmc2.xlsx", 
                                          startRow = 2, 
                                          sheet = .x)) %>% 
    left_join(HUMAN_9606_idmapping[HUMAN_9606_idmapping$Type == "Gene_Name",-2], by = c("Genes" = "ID"))

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
    testing <- KEGGgraph::parseKGML2Graph(tmp)
    test_folate <- KEGGgraph::parseKGML(tmp)
    # mapkG2 <- KEGGpathway2Graph(test_folate, expandGenes=TRUE)
    # set.seed(124)
    # randomNodes <- sample(nodes(mapkG2), 25)
    # mapkGsub <- subGraph(randomNodes, mapkG2)
    # plot(mapkGsub)
    
    for (i in 1:length(testing@nodeData@defaults[["KEGGNode"]][["nodes"]])){
        df <- data.frame(Gene_id = testing@nodeData@defaults[["KEGGNode"]][["nodes"]][[i]]@entryID,
                         Reaction = testing@nodeData@defaults[["KEGGNode"]][["nodes"]][[i]]@reaction,
                         pathway = pathway,
                         stringsAsFactors = FALSE)
        Gene_df <- rbind(Gene_df,df)
        
    }
    
    
    # for (i in 1:length(test_folate@reactions)){
    #     df <- data.frame(Reaction = test_folate@reactions[[i]]@name,
    #                      Substrate = test_folate@reactions[[i]]@substrateName,
    #                      Product = test_folate@reactions[[i]]@productName,
    #                      Direction = test_folate@reactions[[i]]@type,
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

Gene_pathways <- map_df(pathways,Retrieve_Genes_KEGG) %>% 
    left_join(HUMAN_9606_idmapping_hsa, by= c("Gene_id" = "KEGG_ID")) %>%
    subset(!duplicated(Gene.names))

#Raw from http://tko.ccbr.utoronto.ca/ 4/11/2020
Essential_genes <- openxlsx::read.xlsx("./Project_Datasets/reference_essentials_and_nonessentials_sym_hgnc_entrez.xlsx")$Gene
SMPDB_pathways <- read.csv("./Project_Datasets/smpdb_pp_metabolism.proteins.csv", header = F) %>%
    setNames(c("ID","Pathway","Type", "Uniprot","Name","HMDB","Other","mRNA","Gene_Name","Location"))

save(Master_pathways,
     Essential_genes,
     SMPDB_pathways,
     Gene_pathways,
     file = "./Project_Datasets/Metabolic_Pathway_annotations.rda")
