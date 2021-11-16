#This script is part of supplementary documents of "Global Characterization of Fungal Mitogenomes: New Insights on Genomic Diversity and Dynamism of Coding Genes and Accessory Elements"
#published in Frontiers in Microbiology, research topic "Mitochondrial Genomes and Mitochondrion Related Gene Insights to Fungal Evolution"
#DOI: 10.3389/fmicb.2021.787283
#Authors:  Paula L. C. Fonseca, Ruth B. De-Paula, Daniel S. Araújo, Luiz Marcelo R. Tomé, Thairine Mendes-Pereira, Wenderson F. Rodrigues, 
#Luiz-Eduardo Del-Bem, Eric R. Aguiar, and Aristóteles Góes-Neto
#External help: Dener Eduardo Bortolini

#This script uses the results from BLASTn to build similarity networks.  

#******************************************************************************#
#                            Run the code in R		                           #
#******************************************************************************#

library(tidygraph)
library(data.table)
library(ggraph)
library(dplyr)
library(igraph)
library(viridis)

acc <- fread("accessions.m.txt")
dup <- fread("duplicated.blastn.r.txt")
dup <- dup %>% dplyr::group_by(V1,V2) %>% slice(which.max(V4))  #keep largest alignment 
dup <- dup[,1:3]
colnames(dup) <- c("from","to","weight")
dup <- subset(dup, dup$from != dup$to)
dup <- subset(dup, dup$from != "NC_003060.1")
dup <- subset(dup, dup$from != "NC_003061.1")
dup <- left_join(dup, acc, by=c("from"="Accession_number"))
colnames(dup) <- c("from","to","weight","Saccharomyces")
NAs <- subset(dup, is.na(dup$Saccharomyces)) 
NAs <- unique(as.vector(NAs$from))
for (item in NAs){
  dup <- subset(dup, to!=item)
  dup <- subset(dup, from!=item)
}

tbl_links <- dup[,1:3]
tbl_nodes <- distinct(dup[,c(1,4)])
graph2 <- graph_from_data_frame(tbl_links, tbl_nodes, directed = FALSE)

ggraph(graph2, layout = 'fr', weights = weight) + 
  geom_edge_link(aes(color=weight)) + 
  geom_node_voronoi(aes(color=Saccharomyces, fill=Saccharomyces), max.radius = 0.01, expand = unit(-0.5, 'mm')) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4) +
  theme_void() + scale_edge_colour_viridis()

ggsave("duplicated.mitogenomes.blastn.network.connectedonly.sacc.edgeviridis-link-reducednode.pdf",width=7,height=5)

acc <- fread("filos.txt")
all <- fread("general.blastn.r.txt")
all <- all %>% dplyr::group_by(V1,V2) %>% slice(which.max(V4)) #keep largest alignment 
all <- all[,1:3]
colnames(all) <- c("from","to","weight")
all <- subset(all, all$from != all$to)
all <- subset(all, all$from != "NC_003060.1")
all <- subset(all, all$from != "NC_003061.1")

all <- left_join(all, acc, by=c("from"="Accession_number"))
colnames(all) <- c("from","to","weight","Phylum")
NAs <- subset(all, is.na(all$Phylum)) 
NAs <- unique(as.vector(NAs$from))

for (item in NAs){
  all <- subset(all, to!=item)
  all <- subset(all, from!=item)
}

tbl_links <- all[,1:3]
tbl_nodes <- distinct(all[,c(1,4)])
graph2 <- graph_from_data_frame(tbl_links, tbl_nodes, directed = FALSE)

ggraph(graph2, layout = 'fr', weights = weight) + 
  geom_edge_link(aes(color=weight)) + 
  geom_node_point(aes(color=Phylum)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 4) +
  theme_void() + scale_edge_colour_viridis()

ggsave("general.mitogenomes.blastn.network.connectedonly.noSpiz.byfilo_v2.edgeviridis.pdf",width=7,height=5)
