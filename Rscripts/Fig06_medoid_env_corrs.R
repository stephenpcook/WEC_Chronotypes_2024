#This R code generate the analysis and figures used to correlate the chronotype medoids and environmental factors

#In previous codes: cluster_analysis_final.R and cluster_analysis_nonseas_final.R we generated a table with estimated centroids for each chronotype. 
#The files generated contain a list of medoids that were obtain from clustering a chronotype into two groups. Theoretically either of these should represent the chronotypes. 

#That was done using z-scores. What about using the RPKM average to get something closer to a real medioid?
#NOTE I should generate the medioids of the z-score and do the correlations using this normalized value. 

library(dplyr)
library(cluster)
library(ggplot2)
library(reshape2)
library(zoo)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(dtwclust)
library(psych)

#read table
rpkms_3090 <- read.table("3090cat_interpol_rpkms.tab", header=T, row.names=1, check.names=F)

#read environmental metadata
metadata_3090<-read.table("virome_sr_environm_data_EQ_time_nogaps.txt", header=T, row.names=1, check.names=F, sep="\t")
#match the row names to the column names of the rpkms_3090
rownames(metadata_3090)<-paste0("X",1:nrow(metadata_3090),".",rownames(metadata_3090))


#read chronotype metadata
chronotypeMD_3090<-read.table("3090ctgs_metadata_rafah_peakV2_108clsters.txt",header=T, check.names=F, row.names = NULL,sep= "\t")

#cut the first two columns to have a link between final_kluster - contig3090
finalkluster_ctg<-chronotypeMD_3090[,c(1,6)]

#generate a column contig
rpkms_3090$contig<-rownames(rpkms_3090)
rownames(rpkms_3090) <- NULL
rpkms_3090$contig <- gsub("\\.\\.", "||", rpkms_3090$contig)
rpkms_3090$contig <- gsub("\\.", "-", rpkms_3090$contig)
rpkms_3090$contig <- gsub("X", "", rpkms_3090$contig)

# Join the dataframes
joined_data <- inner_join(finalkluster_ctg, rpkms_3090, by = "contig")
  # Initialize an empty dataframe for medoids
  medoids_df <- data.frame(matrix(ncol = ncol(rpkms_3090) - 1, nrow = length(unique(finalkluster_ctg$final_kluster)))) #empty df size of medioid repr (153 as chronotypes) and 32 samples +1 
  colnames(medoids_df) <- colnames(rpkms_3090)[-length(rpkms_3090)]  # Assuming the first column is "Element_column"
  # Initialize a vector to store cluster names
  cluster_names <- character(length = length(unique(finalkluster_ctg$final_kluster)))
 
 #For Loop te generate de medoids of each chronotype
  # Group by Group_column and calculate medoids
  for (i in 1:length(unique(finalkluster_ctg$final_kluster))) {
    group_data <- joined_data %>% filter(final_kluster == unique(finalkluster_ctg$final_kluster)[i])
    # Compute dissimilarity matrix using daisy function from cluster package
    diss_matrix <- daisy(group_data[, -c(1, 2)], metric = "euclidean")
    # Calculate medoid using pam() function from cluster package
    medoid_index <- pam(diss_matrix, diss = TRUE, k = 1)$medoids
    medoids_df[i, ] <- as.numeric(group_data[medoid_index, -c(1, 2)])
    # Store cluster names
    cluster_names[i] <- unique(finalkluster_ctg$final_kluster)[i]}
  # Set row names of the medoids_df to cluster names
  rownames(medoids_df) <- cluster_names

#medoids_df contain the representative of each chronotype. 
#Sanity check: plot them
medoids_df$chronotype<-rownames(medoids_df)
# Create a new column "Seasonality" based on the suffixes in the "Value" column
medoids_df$seasonal_pattern <- ifelse(grepl("ns$", medoids_df$chronotype), "non.seasonal", "seasonal")
#melt
medoids_df_melt<-melt(medoids_df, id.vars=c("chronotype","seasonal_pattern"))

#use zscore 
medoids_df_melt$value_zscore<-zscore(medoids_df_melt$value)

#two panels: seasonal non.seasonal
medoids_df_melt_seasonal<-medoids_df_melt[medoids_df_melt$seasonal_pattern=="seasonal",]
medoids_df_melt_ns<-medoids_df_melt[medoids_df_melt$seasonal_pattern=="non.seasonal",]

######SEASONAL PLOT mediods based on rpkms######
medoids_df_melt_seasonal_no92<-medoids_df_melt_seasonal[medoids_df_melt_seasonal$chronotype != "92s",]
medoid92s<-medoids_df_melt_seasonal[medoids_df_melt_seasonal$chronotype == "92s",]

#plot
plot_all<-ggplot(medoids_df_melt_seasonal,aes(x = variable, y = value_zscore, group=chronotype)) +
geom_line(size=.2,aes( x=variable,y=value_zscore)) +
scale_x_discrete(expand = c(0.01,0))+
geom_vline(xintercept = c(3,15,27), color = "black")+theme_bw()+ facet_wrap(~chronotype)+
guides(color=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.text.x = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x=element_blank())

plot_all_no92<-ggplot(medoids_df_melt_seasonal_no92 ,aes(x = variable, y = value, group=chronotype)) +
geom_line(size=.2,aes( x=variable,y=value)) +
scale_x_discrete(expand = c(0.01,0))+
geom_vline(xintercept = c(3,15,27), color = "black")+theme_bw()+ facet_wrap(~chronotype)+
guides(color=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.text.x = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x=element_blank())

plot_92<-ggplot(medoid92s,aes(x = variable, y = value, group=chronotype)) +
geom_line(size=.2,aes( x=variable,y=value)) +
scale_x_discrete(expand = c(0.01,0))+
geom_vline(xintercept = c(3,15,27), color = "black")+theme_bw()+ facet_wrap(~chronotype)+
guides(color=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.text.x = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x=element_blank())

######Non Seasonal PLOT mediods based on rpkms######
plot_all_ns<-ggplot(medoids_df_melt_ns,aes(x = variable, y = value_zscore, group=chronotype)) +
geom_line(size=.2,aes( x=variable,y=value_zscore)) +
scale_x_discrete(expand = c(0.01,0))+
geom_vline(xintercept = c(3,15,27), color = "black")+theme_bw()+ facet_wrap(~chronotype)+
guides(color=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.text.x = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x=element_blank())

#########################
### Pwise Correlation ###
#########################

#concatenate mediods with env factors
medoids_only_df <- medoids_df[, -c((ncol(medoids_df) - 1):ncol(medoids_df))]
metadata_3090_t<-as.data.frame(t(metadata_3090), stringsAsFactors = FALSE)
medoids_metadata<-rbind(medoids_only_df,metadata_3090_t)

#remove dates from here and interpolate NAs
filtered_medoids_metadata <- medoids_metadata[!(rownames(medoids_metadata) %in% c("Date_wavelet","Date")), , drop = FALSE]
filtered_medoids_metadata[filtered_medoids_metadata == "NA"] <- NA

# Custom function for row-wise interpolation
interpolate_rows <- function(row) {
  interpolated_values <- na.approx(row, na.rm = FALSE)
  return(interpolated_values)
}

# Apply row-wise interpolation while preserving row and column names
filtered_medoids_metadata_interpol <- t(apply(filtered_medoids_metadata, 1, interpolate_rows))
filtered_medoids_metadata_interpol_df <- as.data.frame(filtered_medoids_metadata_interpol, stringsAsFactors = FALSE)

# Restore row and column names
rownames(filtered_medoids_metadata_interpol_df) <- rownames(filtered_medoids_metadata)
colnames(filtered_medoids_metadata_interpol_df) <- colnames(filtered_medoids_metadata)

#UPDATE April 2024 #remove flourescence and density. 
filtered_medoids_metadata_interpol_dfV2 <- filtered_medoids_metadata_interpol_df[!rownames(filtered_medoids_metadata_interpol_df) %in% c("Fluorescence_volts_2.5","Density_kg.m3_2.5"),]

#get the pairwise correlations preserving the row names pair
cor_medoids_metadata <- cor(t(filtered_medoids_metadata_interpol_dfV2), method = "pearson")

#I'll do two plots: only between chronotypes and other chronotypes and env factors
rows_to_split_env_chrono<-row.names(metadata_3090_t)[3:length(row.names(metadata_3090_t))]
rows_to_split_env_chronoV2 <- rows_to_split_env_chrono[!(rows_to_split_env_chrono %in% c("Fluorescence_volts_2.5","Density_kg.m3_2.5"))]

#split th dfs
cor_medoids_metadata_chronos<-cor_medoids_metadata[!(rownames(cor_medoids_metadata) %in% rows_to_split_env_chronoV2), ] #still have the env vs env in the last 10 columns 
cor_medoids_metadata_env<-cor_medoids_metadata[rownames(cor_medoids_metadata) %in% rows_to_split_env_chronoV2, ] #still have the env vs env in the last 10 columns 

#trim the last ten columns
cor_medoids_metadata_env_trim<-cor_medoids_metadata_env[,1:154] #trim the env
cor_enviro_only<-cor_medoids_metadata_env[,155:163]

#get only one small for env factors
cor_enviro_only<-cor_medoids_metadata[rownames(cor_medoids_metadata) %in% rows_to_split_env_chronoV2,colnames(cor_medoids_metadata) %in% rows_to_split_env_chronoV2 ]

###################
#PLOTS with HEATMAP#
###################

md_col<-unique(chronotypeMD_3090[,c(1,7)])#extract kluster, seasonal_pattern, peaking

#unique chronotype
md_col_filtered<- md_col %>% 
  group_by(final_kluster) %>% 
  slice(1)

#get the dat to add the frame
seas_pat<-data.frame(md_col_filtered[,c(1,2)])
rownames(seas_pat)<-seas_pat[,1]
seas_pat <- seas_pat[, -1, drop = FALSE]

#get the gene ratio data into a variable (3 AMGS). First plot the three, then remove primase that I don't care too much. 
#read table
enriched_AMG <- read.table("ratio_enriched_AMGs108.tab", header=T,check.names=F, sep = '\t')

#subset each category
ferro<-enriched_AMG[enriched_AMG$annot=="ferrochelatase",][,c(1,6)]
OG<-enriched_AMG[enriched_AMG$annot=="2OG-Fe(II) oxygenase",] [,c(1,6)]
prim<-enriched_AMG[enriched_AMG$annot=="primase",] [,c(1,6)]

rownames(ferro)<-ferro[,1]
rownames(OG)<-OG[,1]
rownames(prim)<-prim[,1]

ferro<- ferro[, -1, drop = FALSE]
OG <- OG[, -1, drop = FALSE]
prim <- prim[, -1, drop = FALSE]

#Add this info in column
col_seasonal <- HeatmapAnnotation(
  seasonal_pattern = seas_pat$seasonal_pattern,
  #peak=peak_pat$peaking,
  ferrochelatase = anno_barplot(ferro$ratio),
  oxygenase = anno_barplot(OG$ratio),
  col = list(
    seasonal_pattern = c("non.seasonal" = "black", "seasonal"  = "ivory")
   # peak=c("spring"="#4f9710", "summer"="#E08B1A", "autumn"="#937615", "winter"="#255C97", "nd"="black")
)
)

#####CHRONOS#####
# Perform hierarchical clustering
#hc_rows_chronos <- hclust(dist(cor_medoids_metadata_chronos_trim))
#hc_cols_chronos <- hclust(dist(t(cor_medoids_metadata_chronos_trim)))

# Create the row and column dendrograms
#row_dendro_chronos <- as.dendrogram(hc_rows_chronos)
#col_dendro_chronos <- as.dendrogram(hc_cols_chronos)

#heatmap_corr_chronos<-Heatmap(cor_medoids_metadata_chronos_trim,
 # cluster_rows = row_dendro_chronos ,
 # cluster_columns = col_dendro_chronos,
 # column_names_gp = grid::gpar(fontsize = 4),
 # row_names_gp = grid::gpar(fontsize = 4),
 # top_annotation =col_seasonal)

#####CHRONOS-ENV#####
# Perform hierarchical clustering
hc_rows_env <- hclust(dist(cor_medoids_metadata_env_trim))
hc_cols_env <- hclust(dist(t(cor_medoids_metadata_env_trim)))

# Create the row and column dendrograms
row_dendro_env <- as.dendrogram(hc_rows_env)
col_dendro_env <- as.dendrogram(hc_cols_env)

#with z-score -> this might be influencing the correaltions#
#1) join the medioid rpkm and env data
#2) get the z-score
#3) replicate the heatmap
#4) Plot significant correlations with an asterisk

#z-score this one: filtered_medoids_metadata_interpol_df
filtered_medoids_metadata_interpol_df_zscore<- data.frame(zscore(filtered_medoids_metadata_interpol_df))

#get the pairwise correlations preserving the row names pair
cor_medoids_metadata_zscore <- cor(t(filtered_medoids_metadata_interpol_df_zscore), method = "pearson")

# Calculate p-values
p_values <- corr.test(cor_medoids_metadata_zscore)$p #holm correction inputting the matrix of correlations

#hetampa with a function based on the p-values matrix
#also add another constrain for larger 0.5 and lesser than -0.5
heatmap_corr_env<-Heatmap(cor_medoids_metadata_env_trim,
  cluster_rows = row_dendro_env,
  cluster_columns = col_dendro_env,
  column_names_gp = grid::gpar(fontsize = 4),
  row_names_gp = grid::gpar(fontsize = 10),
  top_annotation =col_seasonal,
  cell_fun = function(j, i, x, y, w, h, fill) {
    if(p_values[i, j] < 0.05) {
        grid.text("*", x, y)
    }
  })

#svg("env_chrono_corrs_heatmap_2024V2_108c.svg", width=10, height=4)
heatmap_corr_env
#dev.off()

#From cor_enviro_only generate the cor test
p_values_env<-corr.test(cor_enviro_only)$p #holm correction inputting the matrix of correlation

# Perform hierarchical clustering
hc_rows_env2 <- hclust(dist(cor_enviro_only))
hc_cols_env2 <- hclust(dist(t(cor_enviro_only)))

# Create the row and column dendrograms
row_dendro_env2 <- as.dendrogram(hc_rows_env2)
col_dendro_env2 <- as.dendrogram(hc_cols_env2)

heatmap_corr_env2<-Heatmap(cor_enviro_only,
  cluster_rows = row_dendro_env2,
  cluster_columns = col_dendro_env2,
  column_names_gp = grid::gpar(fontsize = 6),
  row_names_gp = grid::gpar(fontsize = 10),
  cell_fun = function(j, i, x, y, w, h, fill) {
    if(p_values_env[i, j] < 0.05) {
        grid.text("*", x, y)
    }
  })

#svg("env_correls_heatmap_smallpanel2024.svg", width=6, height=6)
heatmap_corr_env2
#dev.off()
