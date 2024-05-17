#VISUALIZATION in R
library(zoo)
library(gtools)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)
library(gtools)
library(reshape2)
library(ggpubr)
library(colorRamp2)
library(ComplexHeatmap)
library(gplots)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(vegan)
library(reshape2)
library(ggpubr)
library(splines2)
library(ggridges)
library(pacman)
library(psych)
library(fpc)
library(factoextra)
library(gplots)
library(gridExtra)
library(corrplot)
library(dtwclust)
library(TSclust)
library(GeneCycle) #fisher.g.test
library(zoo) #interpolation
library(ggrepel)
library(KneeArrower)

#Figure 3 panels (07 in the final version of the manuscipt)
#(a) (3090 chronotypes original): Bowtie 98%, contig length >10,000 bp, vOTUs
#(b) (4588 chronotypes cross): Bowtie 98%, contig length >10,000 bp, cross-assembly
#(c) (5158 chronotypes cross): Bowtie 98%, contig length >5,000 bp, cross-assembly

#(a) (3090 chronotypes original): Bowtie 98%, contig length >10,000 bp, vOTUs
checkv_2<-read.table("decent_seqs.tsv", header=F,check.names=F, sep ="\t")
snps_2<-read.table("SNPcountspergenome.56int.list", header=F,check.names=F, sep ="\t")

colnames(snps_2)<-c("snps","contig")
length_2<-checkv_2[,c(1,2)]
colnames(length_2)<-c("contig","length")

subset_length_2 <- length_2[length_2$contig %in% snps_2$contig, ]

merged2<-merge(subset_length_2,snps_2, by="contig")
merged2$Perc<-(merged2$snps/merged2$length)*100

size_SNPS_plot_2nd<-ggplot(merged2, aes(x = Perc, y = length)) +
geom_point() +  scale_y_continuous(breaks = seq(10000, 130000, by = 25000),limits = c(10000, 130000),labels = scales::comma)+
scale_x_continuous(breaks = seq(0, .7, by = .1),limits = c(0,0.7))+  
labs(x = "Percentage of variable sites", y = "Contig size")+
theme(axis.text.x=element_text(size=13, color = "black"),axis.text.y=element_text(size=13, color = "black"), panel.background = element_blank(),strip.background = element_blank(), legend.title=element_blank(),panel.grid.major.x = element_blank(),panel.grid.major.y =element_blank(), panel.grid.major=element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.title.x = element_text(size = 13,color = "black"),axis.title.y = element_text(size = 13,color = "black"))

#(b) (4588 chronotypes cross): Bowtie 98%, contig length >10,000 bp, cross-assembly
checkv_3<-read.table("merged_shared_RQH_thirdAttempt/decent_seqs.tsv", header=F,check.names=F, sep ="\t")
snps_3<-read.table("merged_shared_RQH_thirdAttempt/SNPcountspergenome.68int.list", header=F,check.names=F, sep ="\t")

colnames(snps_3)<-c("snps","contig")
length_3<-checkv_3[,c(1,2)]
colnames(length_3)<-c("contig","length")

subset_length_3 <- length_3[length_3$contig %in% snps_3$contig, ]

merged3<-merge(subset_length_3,snps_3, by="contig")
merged3$Perc<-(merged3$snps/merged3$length)*100

size_SNPS_plot_3rd<-ggplot(merged3, aes(x = Perc, y = length)) +
geom_point() +  scale_y_continuous(breaks = seq(10000, 130000, by = 25000),limits = c(10000, 130000),labels = scales::comma)+
scale_x_continuous(breaks = seq(0, .7, by = .1),limits = c(0,0.7))+ 
labs(x = "Percentage of variable sites", y = "Contig size")+
theme(axis.text.x=element_text(size=13, color = "black"),axis.text.y=element_text(size=13, color = "black"), panel.background = element_blank(),strip.background = element_blank(), legend.title=element_blank(),panel.grid.major.x = element_blank(),panel.grid.major.y =element_blank(), panel.grid.major=element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.title.x = element_text(size = 13,color = "black"),axis.title.y = element_text(size = 13,color = "black"))

#(c) (5158 chronotypes cross): Bowtie 98%, contig length >5,000 bp, cross-assembly
checkv_4<-read.table("merged_shared_RQH_fourthAttempt/decent_seqs.tsv", header=F,check.names=F, sep ="\t")
snps_4<-read.table("merged_shared_RQH_fourthAttempt/SNPcountspergenome.48int.list", header=F,check.names=F, sep ="\t")

colnames(snps_4)<-c("snps","contig")
length_4<-checkv_4[,c(1,2)]
colnames(length_4)<-c("contig","length")

subset_length_4 <- length_4[length_4$contig %in% snps_4$contig, ]

merged4<-merge(subset_length_4,snps_4, by="contig")
merged4$Perc<-(merged4$snps/merged4$length)*100

size_SNPS_plot_4th<-ggplot(merged4, aes(x = Perc, y = length)) +
geom_point() +  scale_y_continuous(breaks = seq(10000, 130000, by = 25000),limits = c(10000, 130000),labels = scales::comma)+
scale_x_continuous(breaks = seq(0, .7, by = .1),limits = c(0,0.7))+ 
labs(x = "Percentage of variable sites", y = "Contig size")+
theme(axis.text.x=element_text(size=13, color = "black"),axis.text.y=element_text(size=13, color = "black"), panel.background = element_blank(),strip.background = element_blank(), legend.title=element_blank(),panel.grid.major.x = element_blank(),panel.grid.major.y =element_blank(), panel.grid.major=element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),axis.title.x = element_text(size = 13,color = "black"),axis.title.y = element_text(size = 13,color = "black"))

#svg("3panel_ctg_sites.svg", width=10)
plot_grid(size_SNPS_plot_2nd,size_SNPS_plot_3rd,size_SNPS_plot_4th, nrow=1)
#dev.off()

##########
#Figure 3#

#I'll do two alternatives to this figure. 
#(a) Lines and on top normalized RPKM and coverage
#(b) Heatmap and on top normalized RPKM and coverage
#Note we need to introduce NA -> zero columns in the coverage for the missing samples. for RPKM is fine to use the interpolated. 

#READ the tables from the 98% chronotypes 3090
directory <- "merged_shared_RQH_secAttempt/"

# Get the list of files in the directory
file_list <- list.files(directory, full.names = TRUE)

#grep those ending in txt
txt_files <- grep("\\.txt$", file_list, value = TRUE)

# Create an empty list to store the data frames
data_frames <- list()

# Loop to read the table
for (file in txt_files) {
  # Extract the section of the file name
  section <- gsub(".txt", "", basename(file))
  
  # Read the table from the file
  data <- read.table(file, header=T,check.names=F, sep ="\t")

  # Add a new column to the data frame with the name of the file, this will help us to identify which group of IDS belong to what fixed date
  data$fixed_date <- section
  
  # Assign the data frame to a name
  data_frames[[section]] <- data
}

#generate the master DF 
master_df <- do.call(rbind, data_frames) #dim 4988 29 (27 + contig ID + fixed date ID)

#In order to get the time series correct, I need to add the X[num] in front of the string
#This is not interpolated, I will interpolate it. 

#Remove "-viral-fraction" from the name and generate a numerical series in the column names
new_colnames <- gsub("-viral-fraction", "", colnames(master_df))
#vector with the ordered new colnames
complete_colnames_manual<-c("ID","X2.shared_with_2018-Dec","X1.shared_with_2018-Nov","X6.shared_with_2019-Apr","X10.shared_with_2019-Aug","X14.shared_with_2019-Dec","X4.shared_with_2019-Feb","X3.shared_with_2019-Jan","X9.shared_with_2019-Jul","X8.shared_with_2019-Jun","X5.shared_with_2019-Mar","X7.shared_with_2019-May","X13.shared_with_2019-Nov","X11.shared_with_2019-Sep","X22.shared_with_2020-Aug","X26.shared_with_2020-Dec","X16.shared_with_2020-Feb","X15.shared_with_2020-Jan","X21.shared_with_2020-Jul","X17.shared_with_2020-Mar","X25.shared_with_2020-Nov","X24.shared_with_2020-Oct","X23.shared_with_2020-Sep","X30.shared_with_2021-Apr","X27.shared_with_2021-Jan","X32.shared_with_2021-Jun","X29.shared_with_2021-Mar","X31.shared_with_2021-May","fixed_date")
#replace
colnames(master_df)<-complete_colnames_manual

# names of missing samples 
missing_samples <- c("X12.shared_with_2019-Oct","X18.shared_with_2020-Apr","X19.shared_with_2020-May","X20.shared_with_2020-Jun","X28.shared_with_2021-Feb")

# Add missing samples filled with NA
master_df[, missing_samples] <- NA

#Split master_df into two -> numerical and ordered and IS/Fixed_date
master_df_nums<-master_df[,c(2:28,30:34)]
master_df_strings<-master_df[,c(1,29)]

#sort and interpolate to the closest integer master_df_nums
master_df_nums_sorted <- master_df_nums[, mixedsort(names(master_df_nums))]
#Save the colnames
master_df_nums_sorted_cols<-colnames(master_df_nums_sorted)
#interpolate
master_df_nums_sorted_interp <- as.data.frame(t(apply(master_df_nums_sorted, 1, function(row) na.approx(row, na.rm = FALSE, method = "constant", rule = 2.5))))
#add the columns after interpolating
colnames(master_df_nums_sorted_interp)<-master_df_nums_sorted_cols
#merge numbers and text back
final_masterdf_interpol <- cbind(master_df_nums_sorted_interp , master_df_strings)

#melt the master_df
master_df_melt<-melt(final_masterdf_interpol)

#read RPKM norm
NONSignSeasonal_rpkms.normz_melt<-read.table("NONSignSeasonal_rpkms.normz_melt.txt", header=T,check.names=F, sep ="\t")
SignSeasonal_rpkms.normz_melt<-read.table("SignSeasonal_rpkms.normz_melt.txt", header=T,check.names=F, sep ="\t")

rpkms_all<-rbind(SignSeasonal_rpkms.normz_melt,NONSignSeasonal_rpkms.normz_melt)
rpkms_all_df<-dcast(rpkms_all, contig ~ variable)
#sort the columns so after melting this order stays as a factor.
sorted_column_names_rpkm <- mixedsort(names(rpkms_all_df))
sorted_rpkms_all_df<- rpkms_all_df[, sorted_column_names_rpkm]

#melt
rpkms_all_melted<-melt(sorted_rpkms_all_df)
rpkms_all_melted$variable <- gsub("(?<=\\d)\\.(?=[A-Za-z])", "-", rpkms_all_melted$variable, perl = TRUE)

#reorder 
rpkms_all_melted$variable<-factor(rpkms_all_melted$variable, levels = c("X1.2018-Nov","X2.2018-Dec","X3.2019-Jan","X4.2019-Feb","X5.2019-Mar","X6.2019-Apr","X7.2019-May","X8.2019-Jun","X9.2019-Jul","X10.2019-Aug","X11.2019-Sep","X12.2019-Oct","X13.2019-Nov","X14.2019-Dec","X15.2020-Jan","X16.2020-Feb","X17.2020-Mar","X18.2020-Apr","X19.2020-May","X20.2020-Jun","X21.2020-Jul","X22.2020-Aug","X23.2020-Sep","X24.2020-Oct","X25.2020-Nov","X26.2020-Dec","X27.2021-Jan","X28.2021-Feb","X29.2021-Mar","X30.2021-Apr","X31.2021-May","X32.2021-Jun"))

#READ the tables #second attempt 3090 decent and 26k all
directorycov <- "cov_files"
# Get the list of files in the directorycov
file_listcov <- list.files(directorycov, full.names = TRUE)
#grep those ending in txt
txt_filescov <- grep("\\.mod$", file_listcov, value = TRUE)
# Create an empty list to store the data frames
data_framescov <- list()
# Loop to read the table
for (filecov in txt_filescov) {
  # Extract the section of the file name
  sectioncov <- gsub("coverage_table.98filtered.tab.filtJed.", "", basename(filecov))
    # Read the table from the file
  datacov <- read.table(filecov, header=T,check.names=F, sep ="\t")
  # Assign the data frame to a name
  data_framescov[[sectioncov]] <- datacov
}
# Check column names for all data frames. At the end these colnames need to be in the form: X1.2018.Nov so they match with the normalized rpkms
col_names_listcov <- lapply(data_framescov, colnames)
# Merge data frames based on the common column
merged_df_cov <- Reduce(function(x, y) merge(x, y, by = "rname", all = TRUE), data_framescov)
#manually change colnames so I can add the other columns
complete_colnames_manual_cov<-c("contig","X2.2018-Dec","X1.2018-Nov","X6.2019-Apr","X10.2019-Aug","X14.2019-Dec","X4.2019-Feb","X3.2019-Jan","X9.2019-Jul","X8.2019-Jun","X5.2019-Mar","X7.2019-May","X13.2019-Nov","X11.2019-Sep","X22.2020-Aug","X26.2020-Dec","X16.2020-Feb","X15.2020-Jan","X21.2020-Jul","X17.2020-Mar","X25.2020-Nov","X24.2020-Oct","X23.2020-Sep","X30.2021-Apr","X27.2021-Jan","X32.2021-Jun","X29.2021-Mar","X31.2021-May")
colnames(merged_df_cov)<-complete_colnames_manual_cov

# names of missing samples 
missing_samples <- c("X12.2019-Oct","X18.2020-Apr","X19.2020-May","X20.2020-Jun","X28.2021-Feb")
# Add missing samples filled with NA
merged_df_cov[, missing_samples] <- NA

#turn the NAs into 0s
merged_df_cov[is.na(merged_df_cov)] <- 0

#sort the columns so after melting this order stays as a factor.
sorted_column_names <- mixedsort(names(merged_df_cov))
sorted_merged_df_cov<- merged_df_cov[, sorted_column_names]

#melt
merged_df_cov_melted<-melt(sorted_merged_df_cov)


###############
#####PLOTS#####
#Add the coverage by month

#PLOT Function (x~y)
function_plotsTS <- function(vOTU,DataFrame){
#Here subset the values
cluster_df_zscore<-DataFrame[DataFrame$ID == vOTU, ]
#PLOT
ggplot(cluster_df_zscore,aes(x = variable, y = value, group=fixed_date)) +
geom_line(size=.2,aes( x=variable,y=value)) +
scale_x_discrete(expand = c(0.01,0))+
ylim(0, max(cluster_df_zscore$value))+
geom_vline(xintercept = c(3,15,27), color = "black")+
labs(title = unique(cluster_df_zscore$ID))+
theme_bw()+
guides(color=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.text.x = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(size = 7))
  }

function_plotsTS_heatmap <- function(vOTU,DataFrame){
#Here subset the values
cluster_df_zscore<-DataFrame[DataFrame$ID == vOTU, ]
#PLOT
ggplot(cluster_df_zscore, aes(x = variable, y = fixed_date, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +  # You can adjust colors as needed
  scale_x_discrete(expand = c(0.01, 0)) +
  geom_vline(xintercept = c(3, 15, 27), color = "black") +
  labs(title = unique(cluster_df_zscore$ID)) +
  theme_bw() +
  guides(fill = FALSE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 7)
  )
 }

function_plotsTS_cov <- function(vOTU,DataFrame){
#Here subset the values
cluster_df_zscore<-DataFrame[DataFrame$contig == vOTU, ]
#PLOT
ggplot(cluster_df_zscore,aes(x = variable, y = value, group=contig)) +
geom_line(size=1,aes( x=variable,y=value), color="darkorange") +
scale_x_discrete(expand = c(0.01,0))+
ylim(0, max(cluster_df_zscore$value))+
geom_vline(xintercept = c(3,15,27), color = "black")+
labs(title = unique(cluster_df_zscore$contig))+
theme_bw()+
guides(color=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.text.x = element_blank() ,axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(size = 7))
  }

function_plotsTS_rpkm <- function(vOTU,DataFrame){
#Here subset the values
cluster_df_zscore<-DataFrame[DataFrame$contig == vOTU, ]
#PLOT
ggplot(cluster_df_zscore,aes(x = variable, y = value, group=contig)) +
geom_line(size=1,aes( x=variable,y=value), color="navy") +
scale_x_discrete(expand = c(0.01,0))+
ylim(min(cluster_df_zscore$value), max(cluster_df_zscore$value))+
geom_vline(xintercept = c(3,15,27), color = "black")+
labs(title = unique(cluster_df_zscore$contig))+
theme_bw()+
guides(color=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.text.x = element_blank() ,axis.title.y = element_blank(),axis.title.x = element_blank(),plot.title = element_text(size = 7))
  }


#PLOT Function (heatmap, instead of overlaped cov and rpkm put it on top)
fixed_dates<-unique(master_df_melt$fixed_date) #extract the fixed dates
#generate ordered numbers to be attached to the fixed_dates names
to_attach<-c("X2","X1","X6","X10","X13","X4","X3","X9","X8","X5","X7","X12","X11","X18","X22","X15","X14","X17","X16","X21","X20","X19","X25","X23","X27","X24","X26")
#paste them into a df
temp_df <- data.frame(fixed_dates = fixed_dates, order = to_attach)
temp_df$new_fixed_dates<-paste(temp_df$order, temp_df$fixed_dates, sep = ".")

#new dataframe
master_df_melt_orderedHM<-master_df_melt

#Substitute, now we can plot y and x axis in a chronological order.
master_df_melt_orderedHM$fixed_date <- temp_df$new_fixed_dates[match(master_df_melt_orderedHM$fixed_date, temp_df$fixed_dates)]

master_df_melt_orderedHM$fixed_date <- factor(master_df_melt_orderedHM$fixed_date, levels = c("X1.merged_shared_final_2018-Nov","X2.merged_shared_final_2018-Dec","X3.merged_shared_final_2019-Jan","X4.merged_shared_final_2019-Feb", "X5.merged_shared_final_2019-Mar","X6.merged_shared_final_2019-Apr", "X7.merged_shared_final_2019-May", "X8.merged_shared_final_2019-Jun","X9.merged_shared_final_2019-Jul", "X10.merged_shared_final_2019-Aug","X11.merged_shared_final_2019-Sep","X12.merged_shared_final_2019-Nov","X13.merged_shared_final_2019-Dec", "X14.merged_shared_final_2020-Jan","X15.merged_shared_final_2020-Feb","X16.merged_shared_final_2020-Mar", "X17.merged_shared_final_2020-Jul", "X18.merged_shared_final_2020-Aug", "X19.merged_shared_final_2020-Sep","X20.merged_shared_final_2020-Oct","X21.merged_shared_final_2020-Nov","X22.merged_shared_final_2020-Dec","X23.merged_shared_final_2021-Jan","X24.merged_shared_final_2021-Mar","X25.merged_shared_final_2021-Apr","X26.merged_shared_final_2021-May","X27.merged_shared_final_2021-Jun"))

###############################################
#Select the contigs with SNPs greater than 10k#

tmp1<-master_df_melt[,c(1,4)] #extract ID-contig and value (SNPs)

#summarize SNPs per contig
counts_by_ID<-data.frame(tmp1 %>%
group_by(ID) %>%
summarize(Value1_sum = sum(value, na.rm = TRUE)))%>%
  arrange(Value1_sum)

#extract contigs greater than 10k (these are 9)
list10K<-counts_by_ID[counts_by_ID$Value1_sum>10000,]$ID

plot_list <- list()

for (id in list10K){
plot<-function_plotsTS(id,master_df_melt)
plot_list[[id]] <- plot}

#svg("topsubs10_Ts.svg",width=12,height=12)
plot_grob <- arrangeGrob(grobs=plot_list)
grid.arrange(plot_grob)
#dev.off()

###############################################
###Coverage plot of contigs greater than 10k###
plot_listcov <- list()

for (id in list10K){
plotcov<-function_plotsTS_cov(id,merged_df_cov_melted)
plot_listcov[[id]] <- plotcov}

#svg("topsubs10_cov.svg",width=12,height=12)
plot_grob_cov <- arrangeGrob(grobs=plot_listcov)
grid.arrange(plot_grob_cov)
#dev.off()

###############################################
###RPKM plot of contigs greater than 10k###
plot_listrpkm <- list()

for (id in list10K){
plotrpkm<-function_plotsTS_rpkm(id,rpkms_all_melted)
plot_listrpkm[[id]] <- plotrpkm}

#svg("topsubs10_rpkm.svg",width=12,height=12)
plot_grob_rpkm <- arrangeGrob(grobs=plot_listrpkm)
grid.arrange(plot_grob_rpkm)
#dev.off()

###3 examples for Figure 07###
#These have high competness and interesting pattern.
#(a) Most complete:"2019-Nov-viral-fraction__viral-spades-short__as-yet-unknown__56697bp__000078||full ", 
#(b) interesting pattern:"2019-Mar-viral-fraction__viral-spades-short__as-yet-unknown__48511bp__000161||full",
#(c) interesting host (w most completness):2020-Jan-viral-fraction__viral-spades-short__as-yet-unknown__126324bp__000005||full 

#grep from rafah_peak2_metadata the list of the top10
md3090_df<-read.table("3090ctgs_metadata_rafah_peakV2_108clsters.txt",header=T, check.names=F, row.names = NULL,sep= "\t")

#3 to plot:
plot3list<-c("2019-Nov-viral-fraction__viral-spades-short__as-yet-unknown__56697bp__000078||full","2019-Mar-viral-fraction__viral-spades-short__as-yet-unknown__48511bp__000161||full","2020-Jan-viral-fraction__viral-spades-short__as-yet-unknown__126324bp__000005||full") 
md3090_df_top10ksnps<-md3090_df[md3090_df$contig %in% plot3list,]

##
a_rpkm<-function_plotsTS_rpkm("2019-Nov-viral-fraction__viral-spades-short__as-yet-unknown__56697bp__000078||full",rpkms_all_melted)
a_cov<-function_plotsTS_cov("2019-Nov-viral-fraction__viral-spades-short__as-yet-unknown__56697bp__000078||full",merged_df_cov_melted)
a_snps<-function_plotsTS("2019-Nov-viral-fraction__viral-spades-short__as-yet-unknown__56697bp__000078||full",master_df_melt)
a3<-plot_grid(a_rpkm, a_cov, a_snps, align = "hv", ncol=1,rel_heights=c(.15,.15,.7))

b_rpkm<-function_plotsTS_rpkm("2019-Mar-viral-fraction__viral-spades-short__as-yet-unknown__48511bp__000161||full",rpkms_all_melted)
b_cov<-function_plotsTS_cov("2019-Mar-viral-fraction__viral-spades-short__as-yet-unknown__48511bp__000161||full",merged_df_cov_melted)
b_snps<-function_plotsTS("2019-Mar-viral-fraction__viral-spades-short__as-yet-unknown__48511bp__000161||full",master_df_melt)
b3<-plot_grid(b_rpkm, b_cov, b_snps, align = "hv", ncol=1,rel_heights=c(.15,.15,.7))

c_rpkm<-function_plotsTS_rpkm("2020-Jan-viral-fraction__viral-spades-short__as-yet-unknown__126324bp__000005||full",rpkms_all_melted)
c_cov<-function_plotsTS_cov("2020-Jan-viral-fraction__viral-spades-short__as-yet-unknown__126324bp__000005||full",merged_df_cov_melted)
c_snps<-function_plotsTS("2020-Jan-viral-fraction__viral-spades-short__as-yet-unknown__126324bp__000005||full",master_df_melt)
c3<-plot_grid(c_rpkm, c_cov, c_snps, align = "hv", ncol=1,rel_heights=c(.15,.15,.7))

#svg("selection3.svg",width=12,height=11)
plot_grid(a3,b3,c3, align = "v", nrow=1)
#dev.off()

#correlation coverage and rpkm -> Sup con los 3 paneles 

#Generate a dataframe of the 3090 cov, snps, length, etc
#Variant Density= Number of Variants/ (Genome Length×Coverage)
#sorted_merged_df_cov = ~23k coverage

list3090<-md3090_df$contig

sorted_merged_df_cov_3090<-sorted_merged_df_cov[sorted_merged_df_cov$contig %in% list3090, ]

sorted_merged_df_cov_3090_melted<-melt(sorted_merged_df_cov_3090)

#Generate vcf profiles and then merge all the info to generate a variant density
#READ the tables (second attempt 3090 decent and 26k all)
directoryvcf <- "vcf_files"
# Get the list of files in the directorycov
file_listvcf <- list.files(directoryvcf, full.names = TRUE)
#grep those ending in txt
txt_filesvcf <- grep("\\.vcf$", file_listvcf, value = TRUE)
# Create an empty list to store the data frames
data_framesvcf <- list()
# Loop to read the table
for (filevcf in txt_filesvcf) {
  datavcf <- read.table(filevcf, header=T,check.names=F)
  data_framesvcf[[filevcf]] <- datavcf
}

# Check column names for all data frames. At the end these colnames need to be in the form: X1.2018.Nov so they match with the normalized rpkms
col_names_listvcf <- lapply(data_framesvcf, colnames)
# Merge data frames based on the common column
merged_df_vcf <- Reduce(function(x, y) merge(x, y, by = "contig", all = TRUE), data_framesvcf)
#manual colnames so I can add the other columns
complete_colnames_manual_vcf<-c("contig","X2.2018-Dec","X1.2018-Nov","X6.2019-Apr","X10.2019-Aug","X14.2019-Dec","X4.2019-Feb","X3.2019-Jan","X9.2019-Jul","X8.2019-Jun","X5.2019-Mar","X7.2019-May","X13.2019-Nov","X11.2019-Sep","X22.2020-Aug","X26.2020-Dec","X16.2020-Feb","X15.2020-Jan","X21.2020-Jul","X17.2020-Mar","X25.2020-Nov","X24.2020-Oct","X23.2020-Sep","X30.2021-Apr","X27.2021-Jan","X32.2021-Jun","X29.2021-Mar","X31.2021-May")
colnames(merged_df_vcf)<-complete_colnames_manual_vcf

#Subset from the 3090
merged_df_vcf_3090<-merged_df_vcf[merged_df_vcf$contig %in% list3090,] #3058/3090
#turn NA into 0
merged_df_vcf_3090[is.na(merged_df_vcf_3090)] <- 0

#melt
merged_df_vcf_3090_melted<-melt(merged_df_vcf_3090)
colnames(merged_df_vcf_3090_melted)<-c("contig","variable","variants")

#Merge melted vcf,cov,length to generate a variant index:
#remove the interpolated values of missing samples from coverage: 
sorted_merged_df_cov_3090_melted_clean<-sorted_merged_df_cov_3090_melted[! (sorted_merged_df_cov_3090_melted$variable %in% missing_samples), ]

merge_vcf_cov<-merge(merged_df_vcf_3090_melted, sorted_merged_df_cov_3090_melted_clean, by = c("contig", "variable")) #2967 contigs
colnames(merge_vcf_cov)<-c("contig","variable","variants","cov")

md3090_df_reduced<-md3090_df[c(6,5,8,7)]

#remove the interpolated values of missing samples from rpkms:
rpkms_all_melted_reduced<-rpkms_all_melted[! (rpkms_all_melted$variable %in% missing_samples), ]
# Replace second period with hyphen
rpkms_all_melted_reduced$variable <- gsub("(?<=\\d)\\.(?=[A-Za-z])", "-", rpkms_all_melted_reduced$variable, perl = TRUE)

merged_allvalues1<-merge(merge_vcf_cov,md3090_df_reduced, by ="contig")

merged_allvalues<-merge(merged_allvalues1,rpkms_all_melted_reduced,by = c("contig","variable"))

#variant density index
merged_allvalues$denom<-(merged_allvalues$contig_length*merged_allvalues$cov) #56213 has zero values and 23888 nonzero
merged_allvalues$variants_index<-(-log(merged_allvalues$variants/merged_allvalues$denom))
is.na(merged_allvalues)<-sapply(merged_allvalues, is.infinite) 
merged_allvalues[is.na(merged_allvalues)]<-0 #inf into 0

merged_allvalues$variable<-factor(merged_allvalues$variable, levels = c("X1.2018-Nov","X2.2018-Dec","X3.2019-Jan","X4.2019-Feb","X5.2019-Mar","X6.2019-Apr","X7.2019-May","X8.2019-Jun","X9.2019-Jul","X10.2019-Aug","X11.2019-Sep","X12.2019-Oct","X13.2019-Nov","X14.2019-Dec","X15.2020-Jan","X16.2020-Feb","X17.2020-Mar","X18.2020-Apr","X19.2020-May","X20.2020-Jun","X21.2020-Jul","X22.2020-Aug","X23.2020-Sep","X24.2020-Oct","X25.2020-Nov","X26.2020-Dec","X27.2021-Jan","X28.2021-Feb","X29.2021-Mar","X30.2021-Apr","X31.2021-May","X32.2021-Jun"))

vec_no_zeros <- merged_allvalues$variants_index[merged_allvalues$variants_index != 0]

# Find the maximum and minimum values of the column
max_value <- max(vec_no_zeros) #19.03126
min_value <- min(vec_no_zeros) #7.67842

#heatmap Function
function_plots_SNPsHM<-function(vOTU,DataFrame){
#Here subset the values
cluster_df_zscore<-DataFrame[DataFrame$klusters_tree_TsClust108 == vOTU, ]
ggplot(cluster_df_zscore, aes(x = variable, y = contig, fill = variants_index)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = c("white", "white", "white", "lightyellow", "darkorange", "red", "darkred"),  # Colors for the gradient
    values = c(0, 0.16, 0.33, 0.5, 0.67, 0.83, 1),  # Positions of colors along the gradient
    limits = c(0, 20),  # Min and max values of the scale
    breaks = c(0, 7.6, 20),  # Breaks for the color scale
    labels = c("0", "7.6", "20")  # Labels for the breaks
  )  +  # 
  scale_x_discrete(expand = c(0.01, 0)) +
  geom_vline(xintercept = c(2.5,13.5,22.5,27.5), color = "black") +
  theme_bw() +
  #guides(fill = FALSE) + #comment or uncomment depending whether the legend is needed. 
  labs(title = vOTU)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,size=6),
    axis.text.y =element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 6),
    legend.key.size = unit(0.2, "cm"),
    legend.title=element_text(size=3),
    legend.text=element_text(size=5)
  )
}

#x~y correlation snvs and rpkms - All big plot 

#Clear multiseasonal peaking chronotypes
#To study diversity - population from one month to another without being the most abundant. 

#Seasonal chronotypes with defined peaks, how similar are the snv profiles month to month
clearpeaks<-c("16s","20s","24s","25s","51s","52s","54s","62s","65s","70s","72s","79s","83s","87s", "98s", "100s") #16 4X4

merged_allvalues_nonzeroRPKM <- merged_allvalues[merged_allvalues$variants_index != 0, ] 

#remove contigs with zero sum of variability
merged_allvalues_nonzero_SUM <- merged_allvalues %>%
  group_by(contig) %>%
  summarise(total_sum = sum(variants_index)) %>%
  filter(total_sum != 0)

#turn it into a df
merged_allvalues_nonzero_SUM_df<-as.data.frame(merged_allvalues_nonzero_SUM) #from 2967 to 2596
non_zero_sum_contigs<-merged_allvalues_nonzero_SUM_df$contig #get the names to keep

#subset this list from merged_allvalues
merged_allvalues_Divnonzerosum<-merged_allvalues[merged_allvalues$contig %in% non_zero_sum_contigs,]

plot_list <- list()
for (id in clearpeaks){
plot<-function_plots_SNPsHM(id,merged_allvalues_Divnonzerosum)
plot_list[[id]] <- plot}
plot_grob_div <- arrangeGrob(grobs=plot_list)
grid.arrange(plot_grob_div)

#ggsave("heatmaps_divdensity_seas.svg", plot=grid.arrange(plot_grob_div),width=12, height=12)


##x~y correlation snvs index and rpkms for each of the above #use values de rpkm # Local polynomial regression (LOESS)
merged_allvalues_nonzeroRPKM <- merged_allvalues[merged_allvalues$variants_index != 0, ] #remove zero values for variance index

corr_snvs_rpkms<-ggplot(merged_allvalues_nonzeroRPKM, aes(x=variants_index, y=value)) + geom_point(aes(colour=seasonal_pattern),shape=21,size = .1)+scale_fill_manual(values=c("black","white"))+scale_color_manual(values=c("black","black"))+
labs(y= "normalized RPKM", x="variance density")+geom_smooth(method = 'lm')+
theme(axis.text = element_text(size = 12, color="black"),
    legend.position = "none",panel.background = element_blank(),strip.background = element_blank(),
    axis.text.x=element_text(size=10,color = "black"), legend.title=element_blank(),  
    axis.text.y=element_text(size=10,color = "black"),panel.grid.major.x = element_line(size=.1, color="black"),
    panel.grid.major.y =element_line(size=.1, color="black") , panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 12),axis.title.x = element_text(size = 14),
    strip.text.x = element_text(size = 14, color = "black", face = "bold"))+facet_wrap(~seasonal_pattern)+ stat_cor(method = "spearman", label.x = 7.5, label.y = -2.5)

corr_snvs_covs<-ggplot(merged_allvalues_nonzeroRPKM, aes(x=variants_index, y=cov)) + geom_point(aes(colour=seasonal_pattern),shape=21,size = .1)+scale_fill_manual(values=c("black","white"))+scale_color_manual(values=c("black","black"))+
labs(x="variance density",y= "coverage")+geom_smooth(method = 'lm')+
theme(axis.text = element_text(size = 12, color="black"),
    legend.position = "none",panel.background = element_blank(),strip.background = element_blank(),
    axis.text.x=element_text(size=10,color = "black"), legend.title=element_blank(),  
    axis.text.y=element_text(size=10,color = "black"),panel.grid.major.x = element_line(size=.1, color="black"),
    panel.grid.major.y =element_line(size=.1, color="black") , panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 12),axis.title.x = element_text(size = 14),
    strip.text.x = element_text(size = 14, color = "black", face = "bold"))+facet_wrap(~seasonal_pattern)+ stat_cor(method = "spearman", label.x = 7, label.y = 7500)

#svg("xy_snpindex_cov_rpkm.svg", width=4, height=5)
plot_grid(corr_snvs_rpkms,corr_snvs_covs, ncol=1, align='hv')
#dev.off()

######
#Function of cov_snvs and rpkm_snvs to plot the same chronotypes in the HM and see whether covs and rpkms covariate similarly for each group. (more cov - more polymorphisms)

#Correlations plot functions
function_plots_SNPsRPKMS<-function(vOTU,DataFrame){ #vOTU one or a list of chronotype ids
  cluster_df_zscore<-DataFrame[DataFrame$klusters_tree_TsClust108 == vOTU, ]
  cluster_df_zscore <- cluster_df_zscore[cluster_df_zscore$variants_index != 0, ] #remove zero values for variance index
  ggplot(cluster_df_zscore, aes(x=variants_index, y=value)) + geom_point(aes(colour=seasonal_pattern),shape=21,size = .1)+scale_fill_manual(values=c("black","white"))+scale_color_manual(values=c("black","black"))+
labs(y= "normalized RPKM", x="variance density", title=vOTU)+geom_smooth(method = 'lm')+
theme(axis.text = element_text(size = 12, color="black"),
    legend.position = "none",panel.background = element_blank(),strip.background = element_blank(),
    axis.text.x=element_text(size=10,color = "black"), legend.title=element_blank(),  
    axis.text.y=element_text(size=10,color = "black"),panel.grid.major.x = element_line(size=.1, color="black"),
    panel.grid.major.y =element_line(size=.1, color="black") , panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 12),axis.title.x = element_text(size = 14),
    strip.text.x = element_text(size = 14, color = "black", face = "bold"))+ stat_cor(method = "spearman", label.x = 7.5, label.y = -2.5)
}

plot_list_snp_rpkms <- list()
for (id in clearpeaks){
plot<-function_plots_SNPsRPKMS(id,merged_allvalues_Divnonzerosum)
plot_list_snp_rpkms[[id]] <- plot}
plot_grob_snp_rpkms <- arrangeGrob(grobs=plot_list_snp_rpkms)
#svg("corrsRPKM_matching_heatmap.svg", width=12, height=12)
grid.arrange(plot_grob_snp_rpkms)
#dev.off()

function_plots_SNPsCOV<-function(vOTU,DataFrame){ #vOTU one or a list of chronotype ids
  cluster_df_zscore<-DataFrame[DataFrame$klusters_tree_TsClust108 == vOTU, ]
  cluster_df_zscore <- cluster_df_zscore[cluster_df_zscore$variants_index != 0, ] #remove zero values for variance index
  ggplot(cluster_df_zscore, aes(x=variants_index, y=cov)) + geom_point(aes(colour=seasonal_pattern),shape=21,size = .1)+scale_fill_manual(values=c("black","white"))+scale_color_manual(values=c("black","black"))+
labs(y= "coverage", x="variance density", title=vOTU)+geom_smooth(method = 'lm')+
theme(axis.text = element_text(size = 12, color="black"),
    legend.position = "none",panel.background = element_blank(),strip.background = element_blank(),
    axis.text.x=element_text(size=10,color = "black"), legend.title=element_blank(),  
    axis.text.y=element_text(size=10,color = "black"),panel.grid.major.x = element_line(size=.1, color="black"),
    panel.grid.major.y =element_line(size=.1, color="black") , panel.grid.minor = element_blank(), 
    axis.line = element_line(colour = "black"), axis.title.y = element_text(size = 12),axis.title.x = element_text(size = 14),
    strip.text.x = element_text(size = 14, color = "black", face = "bold"))+ stat_cor(method = "spearman", label.x = 7.5, label.y = -2.5)
}

plot_list_snp_cov <- list()
for (id in clearpeaks){
plot<-function_plots_SNPsCOV(id,merged_allvalues_Divnonzerosum)
plot_list_snp_cov[[id]] <- plot}
plot_grob_snp_cov <- arrangeGrob(grobs=plot_list_snp_cov)
#svg("corrsCOV_matching_heatmap.svg", width=12, height=12)
grid.arrange(plot_grob_snp_cov)
#dev.off()

#double check for rpkm chronotypes
maxzscore<-ceiling(max(SignSeasonal_rpkms.normz_melt$value)) #set max 6
minzscore<-floor(min(SignSeasonal_rpkms.normz_melt$value))

function_plotsTS <- function(plots,DataFrame){
cluster<- DataFrame[DataFrame$klusters_tree_list3==plots,]$obs
#Here subset the zscores 
cluster_df_zscore<-rpkms_all_melted[rpkms_all_melted$contig %in% cluster, ]
#PLOT
ggplot(cluster_df_zscore,aes(x = variable, y = value, group=contig)) +
geom_line(size=.2,aes( x=variable,y=value)) +
scale_x_discrete(expand = c(0.01,0))+
scale_y_continuous(limits = c(minzscore, maxzscore))+
labs(title=plots)+
geom_vline(xintercept = c(3,15,27), color = "black")+theme_bw()+
guides(color=FALSE)+theme(plot.title = element_text(size = 5),panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.text.x = element_text(size = 5, angle=90),axis.title.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x=element_blank())
  }

dendro_best_number_contigcluster_df<-read.table("seas_TsClust_contig_memb_opt3.txt",header=T,sep ="\t")

merged_df_vcf_pos<-read.table("vcf_files_positions/vcfpos_merged.txt", header=T,sep ="\t")

#Create a column that matches to the contig list
merged_df_vcf_pos$contig <- gsub("_\\d+$", "", merged_df_vcf_pos$SNP) #2813 contigs.

#subset the metadata of the clearpeaks chronotype. 
clearpeaks_metadata<-md3090_df_reduced[md3090_df_reduced$klusters_tree_TsClust108 %in% clearpeaks,]

#####################
######FIGURE ########

#APRIL 19
#Evaluate whether the second biggest number in a Row is a month before - after or one year later

#(A) extract the 3090 dataset from merged_df_vcf_pos
md3090_contiglist<-md3090_df$contig
merged_df_vcf_pos3090<-merged_df_vcf_pos[merged_df_vcf_pos$contig %in% md3090_contiglist,]

#############################
#PanelA Seasonal Clear Peaks#

#(B) only extract clear peaks
clearpeaks_contigs<-md3090_df[md3090_df$klusters_tree_TsClust108  %in% clearpeaks,]$contig
merged_df_vcf_clearpeaks<-merged_df_vcf_pos[merged_df_vcf_pos$contig %in% clearpeaks_contigs,]

non_contig_columns <- colnames(merged_df_vcf_clearpeaks)[!colnames(merged_df_vcf_clearpeaks) %in% c("contig", "SNP")]

# Initialize a list
results <- list()

# Loop through each column
for(col in non_contig_columns) {
  # Perform comparison of the column with itself
  comparison <- ifelse(merged_df_vcf_clearpeaks[[col]] == 1 & merged_df_vcf_clearpeaks[[col]] == 1, 1, 0)
  result <- data.frame(col1 = col, col2 = col, comparison, contig = merged_df_vcf_clearpeaks$contig)
  results <- append(results, list(result))
}

# Loop through each pair of non-contig columns
for(i in 1:(length(non_contig_columns) - 1)) {
  for(j in (i+1):length(non_contig_columns)) {
    # Get the names of the current pair of columns
    col1 <- non_contig_columns[i]
    col2 <- non_contig_columns[j]
    
    # Perform pairwise comparison and store the result
    comparison <- ifelse(merged_df_vcf_clearpeaks[[col1]] == 1 & merged_df_vcf_clearpeaks[[col2]] == 1, 1, 0)
    result <- data.frame(col1, col2, comparison, contig = merged_df_vcf_clearpeaks$contig)
    results <- append(results, list(result))
    
    # Swap the order of columns for the second comparison
    comparison <- ifelse(merged_df_vcf_clearpeaks[[col2]] == 1 & merged_df_vcf_clearpeaks[[col1]] == 1, 1, 0)
    result <- data.frame(col1 = col2, col2 = col1, comparison, contig = merged_df_vcf_clearpeaks$contig)
    results <- append(results, list(result))
  }
}

# Combine the results into a single data frame
pairwise_comp_df <- do.call(rbind, results)

shared_vcf_bycontig_pairwise <- pairwise_comp_df %>%
  group_by(col1, col2, contig) %>%
  summarise(total_comparison = sum(comparison))

shared_vcf_bycontig_pairwise_df<-data.frame(shared_vcf_bycontig_pairwise)

ranks_time_vcf_clear_peaks<-shared_vcf_bycontig_pairwise_df

#convert to dates
ranks_time_vcf_clear_peaks$col1[] <- lapply(ranks_time_vcf_clear_peaks$col1[], function(x) sub("Jan", "01", x))
ranks_time_vcf_clear_peaks$col1[] <- lapply(ranks_time_vcf_clear_peaks$col1[], function(x) sub("Feb", "02", x))
ranks_time_vcf_clear_peaks$col1[] <- lapply(ranks_time_vcf_clear_peaks$col1[], function(x) sub("Mar", "03", x))
ranks_time_vcf_clear_peaks$col1[] <- lapply(ranks_time_vcf_clear_peaks$col1[], function(x) sub("Apr", "04", x))
ranks_time_vcf_clear_peaks$col1[] <- lapply(ranks_time_vcf_clear_peaks$col1[], function(x) sub("May", "05", x))
ranks_time_vcf_clear_peaks$col1[] <- lapply(ranks_time_vcf_clear_peaks$col1[], function(x) sub("Jun", "06", x))
ranks_time_vcf_clear_peaks$col1[] <- lapply(ranks_time_vcf_clear_peaks$col1[], function(x) sub("Jul", "07", x))
ranks_time_vcf_clear_peaks$col1[] <- lapply(ranks_time_vcf_clear_peaks$col1[], function(x) sub("Aug", "08", x))
ranks_time_vcf_clear_peaks$col1[] <- lapply(ranks_time_vcf_clear_peaks$col1[], function(x) sub("Sep", "09", x))
ranks_time_vcf_clear_peaks$col1[] <- lapply(ranks_time_vcf_clear_peaks$col1[], function(x) sub("Oct", "10", x))
ranks_time_vcf_clear_peaks$col1[] <- lapply(ranks_time_vcf_clear_peaks$col1[], function(x) sub("Nov", "11", x))
ranks_time_vcf_clear_peaks$col1[] <- lapply(ranks_time_vcf_clear_peaks$col1[], function(x) sub("Dec", "12", x))

ranks_time_vcf_clear_peaks$col2[] <- lapply(ranks_time_vcf_clear_peaks$col2[], function(x) sub("Jan", "01", x))
ranks_time_vcf_clear_peaks$col2[] <- lapply(ranks_time_vcf_clear_peaks$col2[], function(x) sub("Feb", "02", x))
ranks_time_vcf_clear_peaks$col2[] <- lapply(ranks_time_vcf_clear_peaks$col2[], function(x) sub("Mar", "03", x))
ranks_time_vcf_clear_peaks$col2[] <- lapply(ranks_time_vcf_clear_peaks$col2[], function(x) sub("Apr", "04", x))
ranks_time_vcf_clear_peaks$col2[] <- lapply(ranks_time_vcf_clear_peaks$col2[], function(x) sub("May", "05", x))
ranks_time_vcf_clear_peaks$col2[] <- lapply(ranks_time_vcf_clear_peaks$col2[], function(x) sub("Jun", "06", x))
ranks_time_vcf_clear_peaks$col2[] <- lapply(ranks_time_vcf_clear_peaks$col2[], function(x) sub("Jul", "07", x))
ranks_time_vcf_clear_peaks$col2[] <- lapply(ranks_time_vcf_clear_peaks$col2[], function(x) sub("Aug", "08", x))
ranks_time_vcf_clear_peaks$col2[] <- lapply(ranks_time_vcf_clear_peaks$col2[], function(x) sub("Sep", "09", x))
ranks_time_vcf_clear_peaks$col2[] <- lapply(ranks_time_vcf_clear_peaks$col2[], function(x) sub("Oct", "10", x))
ranks_time_vcf_clear_peaks$col2[] <- lapply(ranks_time_vcf_clear_peaks$col2[], function(x) sub("Nov", "11", x))
ranks_time_vcf_clear_peaks$col2[] <- lapply(ranks_time_vcf_clear_peaks$col2[], function(x) sub("Dec", "12", x))

ranks_time_vcf_clear_peaks$col1 <-  gsub("^[^.]+\\.", "", ranks_time_vcf_clear_peaks$col1) #remove X[digit]
ranks_time_vcf_clear_peaks$col1 <- gsub("\\.","-",ranks_time_vcf_clear_peaks$col1) #convert . to -

ranks_time_vcf_clear_peaks$col2 <-  gsub("^[^.]+\\.", "", ranks_time_vcf_clear_peaks$col2) #remove X[digit]
ranks_time_vcf_clear_peaks$col2 <- gsub("\\.","-",ranks_time_vcf_clear_peaks$col2) #convert . to -

#add -01 
ranks_time_vcf_clear_peaks$col1 <- paste0(ranks_time_vcf_clear_peaks$col1, "-01")
ranks_time_vcf_clear_peaks$col2 <- paste0(ranks_time_vcf_clear_peaks$col2, "-01")

#generate new date columns
ranks_time_vcf_clear_peaks$date1<-as.Date(ranks_time_vcf_clear_peaks$col1)
ranks_time_vcf_clear_peaks$date2<-as.Date(ranks_time_vcf_clear_peaks$col2)

ranks_time_vcf_clear_peaks$Month_Difference <- round(as.numeric(difftime(ranks_time_vcf_clear_peaks$date2,ranks_time_vcf_clear_peaks$date1, units = "weeks") / 4.34524))

####plot the positive values to see the patterns
ranks_time_vcf_clear_peaks_sub<-ranks_time_vcf_clear_peaks[ranks_time_vcf_clear_peaks$Month_Difference> -1,]
ranks_time_vcf_clear_peaks_sub$Month_Difference <-as.character(ranks_time_vcf_clear_peaks_sub$Month_Difference)
ranks_time_vcf_clear_peaks_sub$total_comparison<-as.numeric(ranks_time_vcf_clear_peaks_sub$total_comparison)
ranks_time_vcf_clear_peaks_sub$Month_Difference<-factor(ranks_time_vcf_clear_peaks_sub$Month_Difference, levels=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"))

ranks_time_vcf_clear_peaks_sub_mean <- ranks_time_vcf_clear_peaks_sub %>% group_by(col1,contig) %>%
  mutate(total_value = max(total_comparison)) %>%
  ungroup() %>%
  # Calculate the percentage
  mutate(percentage = (total_comparison / total_value) * 100)

ranks_time_vcf_clear_peaks_sub_mean_df<-data.frame(ranks_time_vcf_clear_peaks_sub_mean)
ranks_time_vcf_clear_peaks_sub_mean_df[is.na(ranks_time_vcf_clear_peaks_sub_mean_df)] <- 0
ranks_time_vcf_clear_peaks_sub_mean_df$Month_Difference<-as.character(ranks_time_vcf_clear_peaks_sub_mean_df$Month_Difference)
ranks_time_vcf_clear_peaks_sub_mean_df$Month_Difference<-factor(ranks_time_vcf_clear_peaks_sub_mean_df$Month_Difference, levels=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"))

#correlation between month difference and shared polymorphisms. Just confirm the sign is correct in the month difference.

####PLOT####
function_plots_hists<-function(DataFrame){
ggplot(DataFrame, aes(x = Month_Difference, y = percentage, group=Month_Difference)) + geom_boxplot(outlier.shape = NA) + #outliers are not shown to better show the quartile distribution (boxplots). 
labs(y = "Percentage of shared SNPs", x = "Time shift (months)") +
  #geom_vline(xintercept = c(2.5,13.5,22.5,27.5), color = "black") +
  theme_bw() +
  #guides(fill = FALSE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size=10, color="black"),
    axis.text.x = element_text(size=7,color="black"),
    axis.title.y = element_text(size=12,color="black"),
    axis.title.x = element_text(size=12,color="black"),
    plot.title = element_blank())}

#Recover those who have a 12 month signal
#Check the distribution of these values and set a threshold, those with a smaller difference may have a recurrent SNP profile behavoir. 
ranks_time_vcf_clear_peaks_sub_mean_difs_df<-data.frame(ranks_time_vcf_clear_peaks_sub_mean)
ranks_time_vcf_clear_peaks_sub_mean_difs_df[is.na(ranks_time_vcf_clear_peaks_sub_mean_difs_df)] <- 0

#exctract the mean of the percentage per contig per tim-shift and do a dendrogram clustering. #contig/time diference/mean
ranks_new<-data.frame(ranks_time_vcf_clear_peaks_sub_mean_difs_df %>% group_by(contig, Month_Difference) %>%
  summarise(mean_percentage = mean(percentage, na.rm = TRUE)))

dates_sorted<-unique(sort(ranks_time_vcf_clear_peaks_sub_mean_difs_df$date1))
dates_completed<-c(dates_sorted[1:11], "2019-10-01",dates_sorted[12:16],"2020-04-01","2020-05-01","2020-06-01",dates_sorted[17:23],"2021-02-01",dates_sorted[24:27])
month_dif<-as.vector(unique(sort(ranks_time_vcf_clear_peaks_sub_mean_difs_df$Month_Difference)))

df_temporal1 <- data.frame(date_wavelet = dates_completed, Month_Difference= month_dif)

TS_analisis<-merge(ranks_new,df_temporal1 , by="Month_Difference")[2:4]
TS_analisis$contig<-gsub("\\|\\|","..",TS_analisis$contig)
TS_analisis$contig<-gsub("-",".",TS_analisis$contig)

#FUNCTION 2 check for seasonality 
ts_conversor<- function(df, s=c(2018, 11), e=c(2021,6), f=32 ){
  ts(df$mean_percentage, start = s, end=e, frequency = f) 
}

check_seasonality_fisher <- function(ls){
  ls %>% 
    map(~ts_conversor(.x)) %>% 
    map(~data.frame(per = periodogram(.)$freq[
      which.max(periodogram(.)$spec)],
      whole = periodogram(.),
      pval = fisher.g.test(.))) %>% 
  bind_rows( .id = "contig") 
}

#Run the seasonal check
TS_analisis_check <-TS_analisis %>% 
  split(.$contig,drop = T) %>% 
  check_seasonality_fisher()

#select unique rows
TS_analisis_check_uniq<-TS_analisis_check[!duplicated(TS_analisis_check$contig),]

#select contigs names of unique less than 1X10-10
seasonal_TS<-TS_analisis_check_uniq[ which(TS_analisis_check_uniq$pval < 0.01), ][,1]

seasonal_TS<-gsub("\\.\\.","||",seasonal_TS)
seasonal_TS<-gsub("\\.","-",seasonal_TS)

ranks_time_vcf_clear_peaks_SEAS<- ranks_time_vcf_clear_peaks_sub_mean_df[ranks_time_vcf_clear_peaks_sub_mean_df$contig %in% seasonal_TS,]
ranks_time_vcf_clear_peaks_NONS<- ranks_time_vcf_clear_peaks_sub_mean_df[!ranks_time_vcf_clear_peaks_sub_mean_df$contig %in% seasonal_TS,]

mean_variants_all_plot<-function_plots_hists(ranks_time_vcf_clear_peaks_sub_mean_df)
mean_variants_seas_plot<-function_plots_hists(ranks_time_vcf_clear_peaks_SEAS)
mean_variants_nonseas_plot<-function_plots_hists(ranks_time_vcf_clear_peaks_NONS)

peak_final_plot<-plot_grid(mean_variants_all_plot,mean_variants_seas_plot,mean_variants_nonseas_plot, ncol=1, align='hv')

#############################
#PanelB Esporadic CLEAR#
#Seasonal chronotypes with a clear ephemeral pattern
clearephemeral<-c("1ns","2ns","5ns","12ns","16ns","18ns","21ns","22ns","24ns","29ns","30ns","34ns","35ns","38ns","45ns","46ns") #16 4X4
#(B) only extract clear ephemeral
#(B) only extract clear ephemeral
clearephemeral_contigs<-md3090_df[md3090_df$final_kluster %in% clearephemeral,]$contig
merged_df_vcf_clearephemeral<-merged_df_vcf_pos[merged_df_vcf_pos$contig %in% clearephemeral_contigs,]

non_contig_columns_ephe <- colnames(merged_df_vcf_clearephemeral)[!colnames(merged_df_vcf_clearephemeral) %in% c("contig", "SNP")]

# Initialize a list
results_ephe <- list()

# Loop through each column
for(col in non_contig_columns_ephe) {
  # Perform comparison of the column with itself
  comparison <- ifelse(merged_df_vcf_clearephemeral[[col]] == 1 & merged_df_vcf_clearephemeral[[col]] == 1, 1, 0)
  result_ephe<- data.frame(col1 = col, col2 = col, comparison, contig = merged_df_vcf_clearephemeral$contig)
  results_ephe <- append(results_ephe, list(result_ephe))
}

# Loop through each pair of non-contig columns
for(i in 1:(length(non_contig_columns_ephe) - 1)) {
  for(j in (i+1):length(non_contig_columns_ephe)) {
    # Get the names of the current pair of columns
    col1 <- non_contig_columns_ephe[i]
    col2 <- non_contig_columns_ephe[j]
    
    # Perform pairwise comparison and store the result
    comparison <- ifelse(merged_df_vcf_clearephemeral[[col1]] == 1 & merged_df_vcf_clearephemeral[[col2]] == 1, 1, 0)
    result_ephe<- data.frame(col1, col2, comparison, contig = merged_df_vcf_clearephemeral$contig)
    results_ephe <- append(results_ephe, list(result_ephe))
    
    # Swap the order of columns for the second comparison
    comparison <- ifelse(merged_df_vcf_clearephemeral[[col2]] == 1 & merged_df_vcf_clearephemeral[[col1]] == 1, 1, 0)
    result_ephe<- data.frame(col1 = col2, col2 = col1, comparison, contig = merged_df_vcf_clearephemeral$contig)
    results_ephe <- append(results_ephe, list(result_ephe))
  }
}

# Combine the results_ephe into a single data frame
pairwise_comp_df_ephe <- do.call(rbind, results_ephe)

shared_vcf_bycontig_pairwise_ephe <- pairwise_comp_df_ephe %>%
  group_by(col1, col2, contig) %>%
  summarise(total_comparison = sum(comparison))

shared_vcf_bycontig_pairwise_df_ephe<-data.frame(shared_vcf_bycontig_pairwise_ephe)

ranks_time_vcf_clear_ephemeral<-shared_vcf_bycontig_pairwise_df_ephe

#convert to dates
ranks_time_vcf_clear_ephemeral$col1[] <- lapply(ranks_time_vcf_clear_ephemeral$col1[], function(x) sub("Jan", "01", x))
ranks_time_vcf_clear_ephemeral$col1[] <- lapply(ranks_time_vcf_clear_ephemeral$col1[], function(x) sub("Feb", "02", x))
ranks_time_vcf_clear_ephemeral$col1[] <- lapply(ranks_time_vcf_clear_ephemeral$col1[], function(x) sub("Mar", "03", x))
ranks_time_vcf_clear_ephemeral$col1[] <- lapply(ranks_time_vcf_clear_ephemeral$col1[], function(x) sub("Apr", "04", x))
ranks_time_vcf_clear_ephemeral$col1[] <- lapply(ranks_time_vcf_clear_ephemeral$col1[], function(x) sub("May", "05", x))
ranks_time_vcf_clear_ephemeral$col1[] <- lapply(ranks_time_vcf_clear_ephemeral$col1[], function(x) sub("Jun", "06", x))
ranks_time_vcf_clear_ephemeral$col1[] <- lapply(ranks_time_vcf_clear_ephemeral$col1[], function(x) sub("Jul", "07", x))
ranks_time_vcf_clear_ephemeral$col1[] <- lapply(ranks_time_vcf_clear_ephemeral$col1[], function(x) sub("Aug", "08", x))
ranks_time_vcf_clear_ephemeral$col1[] <- lapply(ranks_time_vcf_clear_ephemeral$col1[], function(x) sub("Sep", "09", x))
ranks_time_vcf_clear_ephemeral$col1[] <- lapply(ranks_time_vcf_clear_ephemeral$col1[], function(x) sub("Oct", "10", x))
ranks_time_vcf_clear_ephemeral$col1[] <- lapply(ranks_time_vcf_clear_ephemeral$col1[], function(x) sub("Nov", "11", x))
ranks_time_vcf_clear_ephemeral$col1[] <- lapply(ranks_time_vcf_clear_ephemeral$col1[], function(x) sub("Dec", "12", x))

ranks_time_vcf_clear_ephemeral$col2[] <- lapply(ranks_time_vcf_clear_ephemeral$col2[], function(x) sub("Jan", "01", x))
ranks_time_vcf_clear_ephemeral$col2[] <- lapply(ranks_time_vcf_clear_ephemeral$col2[], function(x) sub("Feb", "02", x))
ranks_time_vcf_clear_ephemeral$col2[] <- lapply(ranks_time_vcf_clear_ephemeral$col2[], function(x) sub("Mar", "03", x))
ranks_time_vcf_clear_ephemeral$col2[] <- lapply(ranks_time_vcf_clear_ephemeral$col2[], function(x) sub("Apr", "04", x))
ranks_time_vcf_clear_ephemeral$col2[] <- lapply(ranks_time_vcf_clear_ephemeral$col2[], function(x) sub("May", "05", x))
ranks_time_vcf_clear_ephemeral$col2[] <- lapply(ranks_time_vcf_clear_ephemeral$col2[], function(x) sub("Jun", "06", x))
ranks_time_vcf_clear_ephemeral$col2[] <- lapply(ranks_time_vcf_clear_ephemeral$col2[], function(x) sub("Jul", "07", x))
ranks_time_vcf_clear_ephemeral$col2[] <- lapply(ranks_time_vcf_clear_ephemeral$col2[], function(x) sub("Aug", "08", x))
ranks_time_vcf_clear_ephemeral$col2[] <- lapply(ranks_time_vcf_clear_ephemeral$col2[], function(x) sub("Sep", "09", x))
ranks_time_vcf_clear_ephemeral$col2[] <- lapply(ranks_time_vcf_clear_ephemeral$col2[], function(x) sub("Oct", "10", x))
ranks_time_vcf_clear_ephemeral$col2[] <- lapply(ranks_time_vcf_clear_ephemeral$col2[], function(x) sub("Nov", "11", x))
ranks_time_vcf_clear_ephemeral$col2[] <- lapply(ranks_time_vcf_clear_ephemeral$col2[], function(x) sub("Dec", "12", x))

ranks_time_vcf_clear_ephemeral$col1 <-  gsub("^[^.]+\\.", "", ranks_time_vcf_clear_ephemeral$col1) #remove X[digit]
ranks_time_vcf_clear_ephemeral$col1 <- gsub("\\.","-",ranks_time_vcf_clear_ephemeral$col1) #convert . to -

ranks_time_vcf_clear_ephemeral$col2 <-  gsub("^[^.]+\\.", "", ranks_time_vcf_clear_ephemeral$col2) #remove X[digit]
ranks_time_vcf_clear_ephemeral$col2 <- gsub("\\.","-",ranks_time_vcf_clear_ephemeral$col2) #convert . to -

#add -01 
ranks_time_vcf_clear_ephemeral$col1 <- paste0(ranks_time_vcf_clear_ephemeral$col1, "-01")
ranks_time_vcf_clear_ephemeral$col2 <- paste0(ranks_time_vcf_clear_ephemeral$col2, "-01")

#generate new date columns
ranks_time_vcf_clear_ephemeral$date1<-as.Date(ranks_time_vcf_clear_ephemeral$col1)
ranks_time_vcf_clear_ephemeral$date2<-as.Date(ranks_time_vcf_clear_ephemeral$col2)

ranks_time_vcf_clear_ephemeral$Month_Difference <- round(as.numeric(difftime(ranks_time_vcf_clear_ephemeral$date2,ranks_time_vcf_clear_ephemeral$date1, units = "weeks") / 4.34524))

####If you think about this one, it is a complete matrix as we did the pairwise comparisons both ways, it is only neccessary to plot the positive values to see the patterns
ranks_time_vcf_clear_ephemeral_sub<-ranks_time_vcf_clear_ephemeral[ranks_time_vcf_clear_ephemeral$Month_Difference> -1,]
ranks_time_vcf_clear_ephemeral_sub$Month_Difference <-as.character(ranks_time_vcf_clear_ephemeral_sub$Month_Difference)
ranks_time_vcf_clear_ephemeral_sub$total_comparison<-as.numeric(ranks_time_vcf_clear_ephemeral_sub$total_comparison)
ranks_time_vcf_clear_ephemeral_sub$Month_Difference<-factor(ranks_time_vcf_clear_ephemeral_sub$Month_Difference, levels=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"))

ranks_time_vcf_clear_ephemeral_sub_mean <- ranks_time_vcf_clear_ephemeral_sub %>% group_by(col1,contig) %>%
  mutate(total_value = max(total_comparison)) %>%
  ungroup() %>%
  # Calculate the percentage
  mutate(percentage = (total_comparison / total_value) * 100)

ranks_time_vcf_clear_ephemeral_sub_mean_df<-data.frame(ranks_time_vcf_clear_ephemeral_sub_mean)
ranks_time_vcf_clear_ephemeral_sub_mean_df[is.na(ranks_time_vcf_clear_ephemeral_sub_mean_df)] <- 0
ranks_time_vcf_clear_ephemeral_sub_mean_df$Month_Difference<-as.character(ranks_time_vcf_clear_ephemeral_sub_mean_df$Month_Difference)
ranks_time_vcf_clear_ephemeral_sub_mean_df$Month_Difference<-factor(ranks_time_vcf_clear_ephemeral_sub_mean_df$Month_Difference, levels=c("0","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31"))

#correlation between month difference and shared polymorphisms. Just confirm the sign is correct in the month difference.

#Recover those who have a 12 month signal
#Check the distribution of these values and set a threshold, those with a smaller difference may have a recurrent SNP profile behavoir. 
ranks_time_vcf_clear_ephemeral_sub_mean_difs_df<-data.frame(ranks_time_vcf_clear_ephemeral_sub_mean)
ranks_time_vcf_clear_ephemeral_sub_mean_difs_df[is.na(ranks_time_vcf_clear_ephemeral_sub_mean_difs_df)] <- 0

#exctract the mean of the percentage per contig per tim-shift and do a dendrogram clustering. #contig/time diference/mean
ranks_new_ephe<-data.frame(ranks_time_vcf_clear_ephemeral_sub_mean_difs_df %>% group_by(contig, Month_Difference) %>%
  summarise(mean_percentage = mean(percentage, na.rm = TRUE)))

dates_sorted_ephe<-unique(sort(ranks_time_vcf_clear_ephemeral_sub_mean_difs_df$date1))
dates_completed<-c(dates_sorted_ephe[1:11], "2019-10-01",dates_sorted_ephe[12:16],"2020-04-01","2020-05-01","2020-06-01",dates_sorted_ephe[17:23],"2021-02-01",dates_sorted_ephe[24:27])
month_dif<-as.vector(unique(sort(ranks_time_vcf_clear_ephemeral_sub_mean_difs_df$Month_Difference)))

df_temporal1_ephe <- data.frame(date_wavelet = dates_completed, Month_Difference= month_dif)

TS_analisis_ephe<-merge(ranks_new_ephe,df_temporal1_ephe , by="Month_Difference")[2:4]
TS_analisis_ephe$contig<-gsub("\\|\\|","..",TS_analisis_ephe$contig)
TS_analisis_ephe$contig<-gsub("-",".",TS_analisis_ephe$contig)

#Run the seasonal check
TS_analisis_ephe_check <-TS_analisis_ephe %>% 
  split(.$contig,drop = T) %>% 
  check_seasonality_fisher()

#select unique rows
TS_analisis_ephe_check_uniq<-TS_analisis_ephe_check[!duplicated(TS_analisis_ephe_check$contig),]

#select contigs names of unique less than 1X10-10
seasonal_TS_ephe<-TS_analisis_ephe_check_uniq[ which(TS_analisis_ephe_check_uniq$pval < 0.01), ][,1]

seasonal_TS_ephe<-gsub("\\.\\.","||",seasonal_TS_ephe)
seasonal_TS_ephe<-gsub("\\.","-",seasonal_TS_ephe)

ranks_time_vcf_clear_ephemeral_SEAS<- ranks_time_vcf_clear_ephemeral_sub_mean_df[ranks_time_vcf_clear_ephemeral_sub_mean_df$contig %in% seasonal_TS_ephe,]
ranks_time_vcf_clear_ephemeral_NONS<- ranks_time_vcf_clear_ephemeral_sub_mean_df[!ranks_time_vcf_clear_ephemeral_sub_mean_df$contig %in% seasonal_TS_ephe,]

mean_variants_all_plot_ephe<-function_plots_hists(ranks_time_vcf_clear_ephemeral_sub_mean_df)
mean_variants_seas_plot_ephe<-function_plots_hists(ranks_time_vcf_clear_ephemeral_SEAS)
mean_variants_nonseas_plot_ephe<-function_plots_hists(ranks_time_vcf_clear_ephemeral_NONS)

ephe_final_plot<-plot_grid(mean_variants_all_plot_ephe,mean_variants_seas_plot_ephe,mean_variants_nonseas_plot_ephe, ncol=1, align='hv')

#FINAL FIGURE:
#svg("FigureFinal_sharedSNP_hists.svg", width=10,height=12)
plot_grid(peak_final_plot,ephe_final_plot, ncol=2, align = 'hv')
#dev.off()

ranks_time_vcf_clear_peaks_SEAS_DF<-data.frame(contig=unique(ranks_time_vcf_clear_peaks_SEAS$contig), membership="seas.signal")
ranks_time_vcf_clear_peaks_NONS_DF<-data.frame(contig=unique(ranks_time_vcf_clear_peaks_NONS$contig), membership="nonseas.signal")

mt_toadd<-rbind(ranks_time_vcf_clear_peaks_SEAS_DF,ranks_time_vcf_clear_peaks_NONS_DF)

clearpeaks_metadatafinal<-merge(clearpeaks_metadata,mt_toadd, by="contig",all = TRUE)
clearpeaks_metadatafinal_clean <- na.omit(clearpeaks_metadatafinal)


################
###HEAT MAPS####
################

#Based on the previous results, let plot an example 

count_memebership_byChrono <-data.frame(clearpeaks_metadatafinal_clean %>%
  group_by(klusters_tree_TsClust108, membership) %>%
  summarise(count = n()))

function_plots_SNPs_contig<-function(vOTU,DataFrame){
#Here subset the values
cluster_df_zscore<-DataFrame[DataFrame$contig == vOTU, ]
cluster_df_zscore_num<-cluster_df_zscore[2:28]
sum_matrix <- outer(cluster_df_zscore_num, cluster_df_zscore_num, Vectorize(sum_matches))
max_value <- max(sum_matrix)
# Normalize the matrix by dividing each value in the matrix by the max value. 
#This will scale the values between 0 and 1, providing a normalized representation of the counts of matches. OR use pearson corr. 
normalized_matrix <- sum_matrix / max_value
melted <- melt(normalized_matrix)

#be carefull here, we need to reorder 
melted$Var1<-factor(melted$Var1, levels = c("X1.2018.Nov","X2.2018.Dec","X3.2019.Jan","X4.2019.Feb","X5.2019.Mar","X6.2019.Apr","X7.2019.May","X8.2019.Jun","X9.2019.Jul","X10.2019.Aug","X11.2019.Sep","X13.2019.Nov","X14.2019.Dec","X15.2020.Jan","X16.2020.Feb","X17.2020.Mar","X21.2020.Jul","X22.2020.Aug","X23.2020.Sep","X24.2020.Oct","X25.2020.Nov","X26.2020.Dec","X27.2021.Jan","X29.2021.Mar","X30.2021.Apr","X31.2021.May","X32.2021.Jun"))
melted$Var2<-factor(melted$Var2, levels = c("X1.2018.Nov","X2.2018.Dec","X3.2019.Jan","X4.2019.Feb","X5.2019.Mar","X6.2019.Apr","X7.2019.May","X8.2019.Jun","X9.2019.Jul","X10.2019.Aug","X11.2019.Sep","X13.2019.Nov","X14.2019.Dec","X15.2020.Jan","X16.2020.Feb","X17.2020.Mar","X21.2020.Jul","X22.2020.Aug","X23.2020.Sep","X24.2020.Oct","X25.2020.Nov","X26.2020.Dec","X27.2021.Jan","X29.2021.Mar","X30.2021.Apr","X31.2021.May","X32.2021.Jun"))

ggplot(melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") + 
  scale_x_discrete(expand = c(0.01, 0)) +
  #geom_vline(xintercept = c(3, 15, 27), color = "black") +
  labs(title = gsub("viral-fraction__viral-spades-short__as-yet-unknown__","",vOTU)) +
  #geom_vline(xintercept = c(2.5,13.5,22.5,27.5), color = "black") +
  theme_bw() +
  guides(fill = FALSE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_blank(),#element_text(size=6,angle=90),
    axis.text.y = element_blank(),#element_text(size=6),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 4.5))}

#md3090_df_reduced and merged_df_vcf_pos as inputs
######PLOT for the total of chronotypes
function_plots_SNPs_chrono<-function(chronotype,DataFrame1,DataFrame2){ 
#Here subset the values
chronotype_list<-DataFrame1[DataFrame1$klusters_tree_TsClust108==chronotype,]$contig
chronotype_vcf<-DataFrame2[DataFrame2$contig %in% chronotype_list,]

cluster_df_zscore_num<-chronotype_vcf[2:28]
sum_matrix <- outer(cluster_df_zscore_num, cluster_df_zscore_num, Vectorize(sum_matches))
max_value <- max(sum_matrix)
# Normalize the matrix by dividing each value by the max value. 
#This will scale the values between 0 and 1, providing a normalized representation of the counts of matches. OR use pearson corr. 
normalized_matrix <- sum_matrix / max_value
melted <- melt(normalized_matrix)

#be carefull here, we need to reorder 
melted$Var1<-factor(melted$Var1, levels = c("X1.2018.Nov","X2.2018.Dec","X3.2019.Jan","X4.2019.Feb","X5.2019.Mar","X6.2019.Apr","X7.2019.May","X8.2019.Jun","X9.2019.Jul","X10.2019.Aug","X11.2019.Sep","X13.2019.Nov","X14.2019.Dec","X15.2020.Jan","X16.2020.Feb","X17.2020.Mar","X21.2020.Jul","X22.2020.Aug","X23.2020.Sep","X24.2020.Oct","X25.2020.Nov","X26.2020.Dec","X27.2021.Jan","X29.2021.Mar","X30.2021.Apr","X31.2021.May","X32.2021.Jun"))
melted$Var2<-factor(melted$Var2, levels = c("X1.2018.Nov","X2.2018.Dec","X3.2019.Jan","X4.2019.Feb","X5.2019.Mar","X6.2019.Apr","X7.2019.May","X8.2019.Jun","X9.2019.Jul","X10.2019.Aug","X11.2019.Sep","X13.2019.Nov","X14.2019.Dec","X15.2020.Jan","X16.2020.Feb","X17.2020.Mar","X21.2020.Jul","X22.2020.Aug","X23.2020.Sep","X24.2020.Oct","X25.2020.Nov","X26.2020.Dec","X27.2021.Jan","X29.2021.Mar","X30.2021.Apr","X31.2021.May","X32.2021.Jun"))

ggplot(melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "navy") + 
  scale_x_discrete(expand = c(0.01, 0)) +
  #geom_vline(xintercept = c(3, 15, 27), color = "black") +
  labs(title = chronotype) +
  #geom_vline(xintercept = c(2.5,13.5,22.5,27.5), color = "black") +
  theme_bw() +
  #guides(fill = FALSE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(size=6,angle=90),
    axis.text.y = element_blank(),#element_text(size=6),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 12))}

#PLOT all the clearpeak
plot_list_clearpeaks <- list()
for (id in clearpeaks){
plotpeaks<-function_plots_SNPs_chrono(id,md3090_df_reduced,merged_df_vcf_pos)
plot_list_clearpeaks[[id]]  <- plotpeaks}
plot_grob_snp_peaks<- arrangeGrob(grobs=plot_list_clearpeaks)
#svg("clearpeaks_merged.svg", width=12, height=10)
grid.arrange(plot_grob_snp_peaks)
dev.off()

#with ggsave
#ggsave("clearpeaks_merged_ggsave.svg", plot=grid.arrange(plot_grob_snp_peaks),width=12, height=10)

#Hacer 20S, 25s, 65s chronotypes. 
#s20<-clearpeaks_metadata[clearpeaks_metadata$klusters_tree_TsClust108=="20s",]$contig
#s25<-clearpeaks_metadata[clearpeaks_metadata$klusters_tree_TsClust108=="25s",]$contig
#s65<-clearpeaks_metadata[clearpeaks_metadata$klusters_tree_TsClust108=="65s",]$contig
#s54<-clearpeaks_metadata[clearpeaks_metadata$klusters_tree_TsClust108=="54s",]$contig
#s100<-clearpeaks_metadata[clearpeaks_metadata$klusters_tree_TsClust108=="100s",]$contig

plot_list_s65 <- list()
for (id in s65){
plot65<-function_plots_SNPs_contig(id,merged_df_vcf_pos)
plot_list_s65[[id]]  <- plot65}
plot_grob_snp_cov65 <- arrangeGrob(grobs=plot_list_s65)
#ggsave("s65_merged_ggsave.svg", plot=grid.arrange(plot_grob_snp_cov65),width=12, height=10)
grid.arrange(plot_grob_snp_cov65)
#dev.off()

plot_list_s20 <- list()
for (id in s20){
plot20<-function_plots_SNPs_contig(id,merged_df_vcf_pos)
plot_list_s20[[id]]  <- plot20}
plot_grob_snp_cov20 <- arrangeGrob(grobs=plot_list_s20)
#ggsave("s20_merged_ggsave.svg", plot=grid.arrange(plot_grob_snp_cov20),width=10, height=8)
grid.arrange(plot_grob_snp_cov20)
#dev.off()

plot_list_s25 <- list()
for (id in s25){
plot25<-function_plots_SNPs_contig(id,merged_df_vcf_pos)
plot_list_s25[[id]]  <- plot25}
plot_grob_snp_cov25 <- arrangeGrob(grobs=plot_list_s25)
##ggsave("s25_merged_ggsave.svg", plot=grid.arrange(plot_grob_snp_cov25),width=12, height=10)
grid.arrange(plot_grob_snp_cov25)
#dev.off()

plot_list_s100 <- list()
for (id in s100){
plot100<-function_plots_SNPs_contig(id,merged_df_vcf_pos)
plot_list_s100[[id]]  <- plot100}
plot_grob_snp_cov100 <- arrangeGrob(grobs=plot_list_s100)
##ggsave("s25_merged_ggsave.svg", plot=grid.arrange(plot_grob_snp_cov25),width=12, height=10)
grid.arrange(plot_grob_snp_cov100)
#dev.off()

#54,70,98
plot_list_s54 <- list()
for (id in s54){
plot54<-function_plots_SNPs_contig(id,merged_df_vcf_pos)
plot_list_s54[[id]]  <- plot54}
plot_grob_snp_cov54 <- arrangeGrob(grobs=plot_list_s54)
##ggsave("s25_merged_ggsave.svg", plot=grid.arrange(plot_grob_snp_cov25),width=12, height=10)
grid.arrange(plot_grob_snp_cov54)
#dev.off()

#####WHISKERPLOT SUP#####
#Time_shift in x axis -10,-9,-8,-7,-6, 0,6,7,10 and y percentage of shared SNPs
#Data are all-versus-all population profiles at each month compared with previous (negative numbers) and subsequent (positive numbers) months
# In master_df_melt add a column with the month difference. 

#Convert to dates 
master_df_melt2<-master_df_melt

master_df_melt2[] <- lapply(master_df_melt2, function(x) sub("merged_shared_final_", "", x))
master_df_melt2[] <- lapply(master_df_melt2, function(x) gsub(".*shared_with_(\\d{4}-\\w+)", "\\1", x))
master_df_melt2[] <- lapply(master_df_melt2, function(x) sub("Jan", "01", x))
master_df_melt2[] <- lapply(master_df_melt2, function(x) sub("Feb", "02", x))
master_df_melt2[] <- lapply(master_df_melt2, function(x) sub("Mar", "03", x))
master_df_melt2[] <- lapply(master_df_melt2, function(x) sub("Apr", "04", x))
master_df_melt2[] <- lapply(master_df_melt2, function(x) sub("May", "05", x))
master_df_melt2[] <- lapply(master_df_melt2, function(x) sub("Jun", "06", x))
master_df_melt2[] <- lapply(master_df_melt2, function(x) sub("Jul", "07", x))
master_df_melt2[] <- lapply(master_df_melt2, function(x) sub("Aug", "08", x))
master_df_melt2[] <- lapply(master_df_melt2, function(x) sub("Sep", "09", x))
master_df_melt2[] <- lapply(master_df_melt2, function(x) sub("Oct", "10", x))
master_df_melt2[] <- lapply(master_df_melt2, function(x) sub("Nov", "11", x))
master_df_melt2[] <- lapply(master_df_melt2, function(x) sub("Dec", "12", x))

#add -01 
master_df_melt2$fixed_date <- paste0(master_df_melt2$fixed_date, "-01")
master_df_melt2$variable <- paste0(master_df_melt2$variable, "-01")


#generate new date columns
master_df_melt2$fixed_date2<-as.Date(master_df_melt2$fixed_date)
master_df_melt2$variable2<-as.Date(master_df_melt2$variable)

#estimate the difference
master_df_melt2$Month_Difference <- round(as.numeric(difftime(master_df_melt2$variable2,master_df_melt2$fixed_date2, units = "weeks") / 4.34524))
#grep all the rows which values in column X are defined between -12 and 12
master_df_melt2subs <- master_df_melt2[master_df_melt2$Month_Difference >= -12 & master_df_melt2$Month_Difference <= 12, ]
master_df_melt2subs$ID <- as.factor(master_df_melt2subs$ID)
master_df_melt2subs$value <- as.numeric(master_df_melt2subs$value)

master_df_melt2subs_mean <- master_df_melt2subs %>% group_by(ID,fixed_date) %>%
  mutate(total_value = max(value)) %>%
  ungroup() %>%
  # Calculate the percentage
  mutate(percentage = (value / total_value) * 100)

master_df_melt2subs_mean_df<-data.frame(master_df_melt2subs_mean)
master_df_melt2subs_mean_df$Month_Difference<-as.character(master_df_melt2subs_mean_df$Month_Difference)

whiskerplot<-ggplot(master_df_melt2subs_mean, aes(x = Month_Difference, y = percentage, group=Month_Difference)) + geom_boxplot(coef = 2, outlier.shape = NA) +
labs(y = "Percentage of shared SNPs", x = "Time shift (months)") +
  #geom_vline(xintercept = c(2.5,13.5,22.5,27.5), color = "black") +
  theme_bw() +
  #guides(fill = FALSE) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(size=9),
    axis.text.x = element_text(size=9),
    axis.title.y = element_text(size=12),
    axis.title.x = element_text(size=12),
    plot.title = element_blank())

#svg("whiskerplot.svg", width=8, height=6)
whiskerplot
dev.off()
