# BiocManager::install(c("TreeAndLeaf","RedeR"))
# install.packages(c("igraph","RColorBrewer"))

#This R code generate the figure 05 of the manuscript: Viral chronotypes and their role in shaping seasonal viral dynamics in the Western English Channel

#Load packages
library("TreeAndLeaf")
library("RedeR")
library("igraph")
library("ape")
library("ggtree")
library("dendextend")
library("dplyr")
library("ggplot2")
library("RColorBrewer")
library("tidytree")
library("cowplot")
library("ggnewscale")

#################
#####PRIMASE#####

#read the two files: metadata and newick tree
#Read metadata, this file was specifically tailored for the contigs 
primase_md <- read.table("primase_genomes_metadata_wSAR11.txt", header=T, check.names=F, sep="\t")

#Read Vcontact2 taxonomic classification 
peak<-read.table("3090ctgs_metadata_rafah_peak.txt", header=T, check.names=F, row.names = NULL,sep= "\t")

#subset contig and host_for_tree
ctg_viral_family<-peak[,c(2,27)]

#merge by the contig column
primase_md_fam<- merge(primase_md, ctg_viral_family, by = "contig")

#Read tree
tree_primase<-ape::read.tree("pharokka3090_primase_ChronoSAR11.aln.fasta.treefile")
#reformat to treedata
tree_primase_data <- as.treedata(tree_primase)

#Subset to create a new dataframe 
merged_df_reduced<- data.frame(label=primase_md_fam$label, group=primase_md_fam$final_kluster, seasonal=primase_md_fam$seasonal_pattern, family=primase_md_fam$host_for_tree)

#join tree with metadata 
phylo_tree_prim <- dplyr::full_join(tree_primase_data,merged_df_reduced, by='label')

#Primase panel #layout equal_angle for aesthetics and circular to double check the clusterings
p2<-ggtree(phylo_tree_prim, layout="equal_angle", aes(color=family))+ scale_color_manual(values=c("Pelagibacter"="orangered3","other"="gray50","PelagibacterVC2"="purple3"))#+theme(legend.position="none")
primase_final<-p2+ new_scale_color()+geom_tippoint(aes(fill=seasonal), shape=21, size=1)+ scale_fill_manual(values=c("seasonal"="ivory","non.seasonal"="black","pelagiphage_reference"="orangered3"))#+theme(legend.position="none")

#########################
######Ferrochelatase#####
#########################

#Read metadata, this file was specifically tailored for the contigs 
Fe_md <- read.table("Fe_genomes_metadata_wSAR11.txt", header=T, check.names=F, sep="\t")

#merge by the contig column
Fe_md_fam<- merge(Fe_md, ctg_viral_family, by = "contig")

#Read tree
tree_Fe<-ape::read.tree("pharokka3090_Fe_ChronoSAR11.aln.fasta.treefile")
#reformat to treedata
tree_Fe_data <- as.treedata(tree_Fe)

#Subset to create a new dataframe 
merged_df_Fe_reduced<- data.frame(label=Fe_md_fam$label, group=Fe_md_fam$final_kluster, seasonal=Fe_md_fam$seasonal_pattern, family=Fe_md_fam$host_for_tree)

#join tree with metadata 
phylo_tree_Fe <- dplyr::full_join(tree_Fe_data,merged_df_Fe_reduced, by='label')

#Ferrochelatase panel
p3<-ggtree(phylo_tree_Fe, layout="equal_angle", aes(color=family))+ scale_color_manual(values=c("Pelagibacter"="orangered3","other"="gray50","PelagibacterVC2"="purple3"))#+theme(legend.position="none")
Fe_final<-p3+ new_scale_color()+geom_tippoint(aes(fill=seasonal), shape=21, size=1)+ scale_fill_manual(values=c("seasonal"="ivory","non.seasonal"="black","pelagiphage_reference"="orangered3"))#+theme(legend.position="none")

#########################
######2OG oxygenase #####
#########################
OG_md <- read.table("2OG_genomes_metadata_wSAR11.txt", header=T, check.names=F, sep="\t")

#merge by the contig column
OG_md_fam<- merge(OG_md, ctg_viral_family, by = "contig")

#Read tree
tree_OG<-ape::read.tree("pharokka3090_2OG_ChronoSAR11.aln.fasta.treefile")
#reformat to treedata
tree_OG_data <- as.treedata(tree_OG)

#Subset to create a new dataframe 
merged_df_OG_reduced<- data.frame(label=OG_md_fam$CDS, group=OG_md_fam$final_kluster, seasonal=OG_md_fam$seasonal_pattern, family=OG_md_fam$host_for_tree)

#join tree with metadata 
phylo_tree_OG <- dplyr::full_join(tree_OG_data,merged_df_OG_reduced, by='label')

#OG panel
p1<-ggtree(phylo_tree_OG, layout="equal_angle", aes(color=family))+ scale_color_manual(values=c("Pelagibacter"="orangered3","other"="gray50","PelagibacterVC2"="purple3"))#+theme(legend.position="none")
OG_final<-p1+ new_scale_color()+geom_tippoint(aes(fill=seasonal), shape=21, size=1)+ scale_fill_manual(values=c("seasonal"="ivory","non.seasonal"="black","pelagiphage_reference"="orangered3"))#+theme(legend.position="none")


##FINAL PLOT
#svg("phylog_3fridwSAR11.svg", height=10)
plot_grid(OG_final,Fe_final,primase_final, align = "v", ncol = 1, rel_heights=c(.5,.35,.2))
#dev.off()

#Print them individually
#svg("/Users/lb808/Documents/MyDrafts/InProgress/01_WEC_Viromes_ShortReads/Figures/04_ProtFamPhylog/phylog_OG.svg", height=12, width=12)
OG_final
#dev.off()

#Print them individually
#svg("/Users/lb808/Documents/MyDrafts/InProgress/01_WEC_Viromes_ShortReads/Figures/04_ProtFamPhylog/phylog_Fe.svg", height=12, width=12)
Fe_final
#dev.off()

#Print them individually
#svg("/Users/lb808/Documents/MyDrafts/InProgress/01_WEC_Viromes_ShortReads/Figures/04_ProtFamPhylog/primase_final.svg", height=12, width=12)
primase_final
#dev.off()

###The END###
