#This R code generate the Figure 04 of Bolanos et al., 2024

#After annotating the vOTUs sequences, these were analyzed within the ChronoClustR framework. 

#Load libraries#
library(dplyr)
library(gplots)
library(tidyr)
library(dendextend)
library(ggplot2)
library(cowplot)
library(ggplotify) 

#Define half violin
# Half violin plot function:

"%||%" <- function(a, b) {
  if (!is.null(a))
    a
  else
    b
}

#function definitions#

geom_flat_violin <-
  function(mapping = NULL,
           data = NULL,
           stat = "ydensity",
           position = "dodge",
           trim = TRUE,
           scale = "area",
           show.legend = NA,
           inherit.aes = TRUE,
           ...) {
    ggplot2::layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomFlatViolin,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(trim = trim,
                    scale = scale,
                    ...)
    )
  }

GeomFlatViolin <-
  ggproto(
    "GeomFlatViolin",
    Geom,
    setup_data = function(data, params) {
      data$width <- data$width %||%
        params$width %||% (resolution(data$x, FALSE) * 0.9)
      
      # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
      data %>%
        dplyr::group_by(.data = ., group) %>%
        dplyr::mutate(
          .data = .,
          ymin = min(y),
          ymax = max(y),
          xmin = x,
          xmax = x + width / 2
        )
    },
    
    draw_group = function(data, panel_scales, coord)
    {
      # Find the points for the line to go all the way around
      data <- base::transform(data,
                              xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
      
      # Make sure it's sorted properly to draw the outline
      newdata <-
        base::rbind(
          dplyr::arrange(.data = base::transform(data, x = xminv), y),
          dplyr::arrange(.data = base::transform(data, x = xmaxv), -y)
        )
      
      # Close the polygon: set first and last point the same
      # Needed for coord_polar and such
      newdata <- rbind(newdata, newdata[1,])
      
      ggplot2:::ggname("geom_flat_violin",
                       GeomPolygon$draw_panel(newdata, panel_scales, coord))
    },
    
    draw_key = draw_key_polygon,
    
    default_aes = ggplot2::aes(
      weight = 1,
      colour = "grey20",
      fill = "white",
      size = 0.5,
      alpha = NA,
      linetype = "solid"
    ),
    
    required_aes = c("x", "y")
  )

#End of half violin#

####################### 
#Read categories table#
#######################

kluster_cat_df<-read.table("kluster_cat2.tsv",header=T, check.names=F, row.names = NULL,sep= "\t",comment.char = "")
kluster_annot_df<-read.table("kluster_anot2.tsv",header=T, check.names=F, row.names = NULL,sep= "\t",comment.char = "")
kluster_phrog_df<-read.table("kluster_phrog2.tsv",header=T, check.names=F, row.names = NULL,sep= "\t",comment.char = "")

##################
####Enrichment####
##################

#####BY ANNOT#####
md3090_df<-read.table("3090ctgs_metadata_rafah_peakV2_108clsters.txt",header=T, check.names=F, row.names = NULL,sep= "\t")

ctg_kluster<-md3090_df[1:6]
colnames(ctg_kluster)<-c("kluster","contig")

#Generate contig counts per chronotype
category_counts <- table(ctg_kluster$kluster)

#Create a new data frame with unique values and their counts
ctg_cluster_df <- data.frame(
  final_kluster = names(category_counts),
  count = as.vector(category_counts)
)

#merge 
merged_annot_df <- merge(kluster_annot_df, ctg_cluster_df , by = "final_kluster", all.x = TRUE)

#get the ratio
merged_annot_df$ratio<-merged_annot_df$count.x/merged_annot_df$count.y #dim 57834x6

#remove all zeros 
merged_annot_df_non0<-merged_annot_df[merged_annot_df$ratio!=0,]

#remove all hypothetical and unknown
merged_annot_df_non0_hyp<-merged_annot_df_non0[merged_annot_df_non0$annot!="hypothetical protein",]
merged_annot_df_non0_hyp_unk<-merged_annot_df_non0_hyp[merged_annot_df_non0_hyp$annot!="unknown function",]

######ORIGINAL RATIO ########
#Get the unique annotation categories in the data
unique_annot_ratio <- unique(merged_annot_df$annot)

#Initialize a vector to store adjusted p-values
adjusted_p_values_ratio <- vector("numeric", length(unique_annot_ratio))

#Loop through each unique anot and perform t-test with Bonferroni correction
for (i in 1:length(unique_annot_ratio)) {
  annot_ratio <- unique_annot_ratio[i]
  
  annot_subset_ratio <- merged_annot_df[merged_annot_df$annot == annot_ratio, ]
  
  seasonal_occurrences_ratio <- annot_subset_ratio[annot_subset_ratio$pattern == "seasonal", "ratio"]
  not_seasonal_occurrences_ratio <- annot_subset_ratio[annot_subset_ratio$pattern == "non_seasonal", "ratio"]
  
  t_test_result_ratio <- t.test(seasonal_occurrences_ratio, not_seasonal_occurrences_ratio)
  
  adjusted_p_values_ratio[i] <- t_test_result_ratio$p.value * length(unique_annot_ratio)  #Bonferroni correction
}

#Check for significant anot after Bonferroni correction
significant_anot_ratio <- unique_annot_ratio[adjusted_p_values_ratio <= 0.2]

number_sign_res<-length(significant_anot_ratio)

#Extract the p-vales from 
sign_ratios_to_add<-sort(adjusted_p_values_ratio)[1:number_sign_res]

merged_annot_significant<-merged_annot_df[merged_annot_df$annot %in% significant_anot_ratio, ] #subset and create a dataframe of the significant different categories

#Add a Z to the new phosphoheptose isomerase panel to keep the order
merged_annot_significant$annot <- gsub("phosphoheptose isomerase", "zphosphoheptose isomerase", merged_annot_significant$annot)

#BOXPLOTS of the ratio distributions of the annotated categories VERSION April 2024#

plot_sing_seas_noseas<-ggplot(merged_annot_significant, aes(x =pattern, y = ratio))+labs(y= "gene observation / chronotype members", x='')+
theme(axis.text = element_text(size = 14, color="black"), strip.text.x = element_text(size = 14, color = "black", face = "bold"),panel.background = element_blank(),axis.text.x=element_text(size=14,color = "black",vjust = .5), legend.title=element_blank(),  axis.text.y=element_text(vjust=0.4, hjust=1,size=14,color = "black"),panel.grid.major.x = element_line(size=.1,color="black"),panel.grid.major.y =element_line(size=.1, color="black") , panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), strip.background = element_blank(),axis.title.y = element_text(size = 14))

boxplot_seas_noseas_Nov2023<-plot_sing_seas_noseas + stat_boxplot(position = position_nudge(x = 0, y = 0), geom = "errorbar", width = 0.15) + 
geom_flat_violin(position = position_nudge(x = 0, y = 0), adjust=2, scale = "width", fill = "lightgray", colour = "gray")+
geom_boxplot(outlier.shape = NA, alpha = 0.3, width = .5, colour = "black")+
geom_point(position = position_jitter(width = .075), size = .75, shape=21)+facet_grid(. ~ annot)

#svg("significant_seas_non_seas_Apr24V2TSCLUSTonly.svg", width=18, height=5)
boxplot_seas_noseas_Nov2023
#dev.off()

####Finally write the table for correlations figure 06####
#write.table(merged_annot_significant, "ratio_enriched_AMGs108.tab", sep="\t", row.names=FALSE, quote=FALSE)
