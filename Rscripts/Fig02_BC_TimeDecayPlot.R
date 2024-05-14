#Description: Figure 02 top panel and Sup Figures 01 and 02 shows multiple pannels depicting Bray-Curtis dissimilarity time-decay plots. 
#We interrogated whether just a few dominant viral populations are driving the annual community turnover oberserved in this study and previous. 

library(phyloseq)
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
library(gridExtra)
library(corrplot)
library(ape)
library(BBmisc)
library(dtwclust)
library(TSclust)
library(zoo)
library(GeneCycle) 

#PLOT 1 BRAY CURTIS TIME DECAY PLOT

#Generate a Phyloseq object
countwec <- read.table("rpkms_populations.txt", header=T, row.names=1, check.names=F)
row.names(countwec) <- sub("^", "X", row.names(countwec))
row.names(countwec) <- sub("[||]+", "..",row.names(countwec) )
row.names(countwec)  <- gsub("-", ".", row.names(countwec) )

env_sr_vir<-read.table("virome_sr_environm_data.txt",header=T, row.names=1, check.names=F, sep ="\t")

#Use as date
env_sr_vir$Date_f <- as.Date(env_sr_vir$Date, "%d/%m/%Y")

###For top panel of figure 02 we will use only the 3090 contigs to be consistent with the overall analysis. 

countwec3090<-read.table("rpkms_populations_3090decent.txt", header=T, row.names=1, check.names=F)
row.names(countwec3090) <- sub("^", "X", row.names(countwec3090))
row.names(countwec3090) <- sub("[||]+", "..",row.names(countwec3090) )
row.names(countwec3090)  <- gsub("-", ".", row.names(countwec3090) )

OTU<-otu_table(countwec, taxa_are_rows=TRUE)
OTU2<-otu_table(countwec3090, taxa_are_rows=TRUE)
SAM<-sample_data(env_sr_vir)
PHYSEQ<-phyloseq(OTU2,SAM) #switch between OTU or OTU2 depending of which dataset will be plotted

#####Split the PHYSEQ by abundances

#GENERATE the distribution of the row sums to determine the quartiles from the distribution of the total abundances
df_sorted<-data.frame(sort(taxa_sums(PHYSEQ)))
colnames(df_sorted)<-"RPKM"

#extract the thresholds
perc25<-quantile(df_sorted$RPKM)[[2]]
perc75<-quantile(df_sorted$RPKM)[[4]] #remove the name (i.e. 0%) and extract the value to use it to subset the PHYSEQ

#subset a vector of names based on the threshelds <25%, <=25% & <=50%, >50%
high_RPKMS<-rownames(subset(df_sorted, RPKM < perc25))
medium_RPKMS<-rownames(subset(df_sorted, RPKM >= perc25 & RPKM <= perc75))
low_RPKMS<-rownames(subset(df_sorted, RPKM > perc75))

#sanity check: sum of elements of the three subsets should be equal to the total in the dataframe.
length(high_RPKMS)+length(medium_RPKMS)+length(low_RPKMS) #26851
dim(df_sorted) #26851

#Subset the generated lists from the phyloseq object
high_subset <- subset(otu_table(PHYSEQ), rownames(otu_table(PHYSEQ)) %in% high_RPKMS)
high_physeq <- merge_phyloseq(high_subset , sample_data(PHYSEQ))

med_subset <- subset(otu_table(PHYSEQ), rownames(otu_table(PHYSEQ)) %in% medium_RPKMS)
med_physeq <- merge_phyloseq(med_subset , sample_data(PHYSEQ))

low_subset <- subset(otu_table(PHYSEQ), rownames(otu_table(PHYSEQ)) %in% low_RPKMS)
low_physeq <- merge_phyloseq(low_subset , sample_data(PHYSEQ))

##Each of these Phyloseq objects I neet to evaluate them to determine the seasonal and not seasonal contigs
#FUNCTION 1: convert our datasets into a time-series object
ts_conversor<- function(df, s=c(2018, 11), e=c(2021,6), f=32 ){
  ts(df$value, start = s, end=e, frequency = f) 
}

#FUNCTION 2 check for seasonality 
check_seasonality_fisher <- function(ls){
  ls %>% 
    map(~ts_conversor(.x)) %>% 
    map(~data.frame(per = periodogram(.)$freq[
      which.max(periodogram(.)$spec)],
      whole = periodogram(.),
      pval = fisher.g.test(.))) %>% 
  bind_rows( .id = "variable") 
}

#READ INPUT FILES with interpolated rpkms and the environmental data with dates fixed to be equidistant
#Exctract whih ones are found in the seasonal and which ones in the non seasonal fraction. 

#All contigs
countwec <- read.table("rpkms_populations.txt", header=T, row.names=1, check.names=F)
countwect<-t(countwec)

#Decent contigs
countwec3090<-read.table("rpkms_populations_3090decent.txt", header=T, row.names=1, check.names=F)
countwec3090t<-t(countwec3090)

VIR_cl_env1<-read.table("virome_sr_environm_data_EQ_time.txt", header=T, row.names=1, check.names=F, sep ="\t")
VIR_cl_env1$sample<-rownames(VIR_cl_env1)
VIR_cl_envdate<-VIR_cl_env1[ , c("sample","Date_wavelet")]
VIR_cl_envdate$sample<-paste0("X",1:nrow(VIR_cl_envdate),".",VIR_cl_envdate$sample)

monthly_RPKM<-merge(VIR_cl_envdate,countwec3090t,by=0, all = TRUE) ######I need to interpolate this
monthly_RPKM$Date_wavelet<-as.Date(monthly_RPKM$Date_wavelet,"%d/%m/%Y")
monthly_RPKM_sorted<-(monthly_RPKM[order(monthly_RPKM$Date_wavelet),]) #dim 32 26854 or 32 3093 for the decent subset

#Interpolation of the numerical columns. I;m skipping 1:3 becuase these are other alphanumerical strings (Row.names-sample-Date_wavelet)
monthly_RPKM_sort_trim_int<- data.frame(monthly_RPKM_sorted[1:3], na.approx(monthly_RPKM_sorted[4:length(monthly_RPKM_sorted)])) #na.approx fn -> interpolation
monthly_RPKM_sort_trim_int_cl<-monthly_RPKM_sort_trim_int[,c(3:length(monthly_RPKM_sorted))]
monthly_RPKM_sort_trim_int_cl_melt<-melt(monthly_RPKM_sort_trim_int_cl, id.vars=c("Date_wavelet"))

#Run the seasonal check
seasonal.check.RPKM <-monthly_RPKM_sort_trim_int_cl_melt %>% 
  split(.$variable,drop = T) %>% 
  check_seasonality_fisher()

#select unique rows
seasonal.check.RPKM_uniq<-seasonal.check.RPKM[!duplicated(seasonal.check.RPKM$variable),]

#select contigs names of unique less than 1X10-10
highly_seasonaldf<-seasonal.check.RPKM_uniq[ which(seasonal.check.RPKM_uniq$pval < 0.05), ][,1]
non_seasonaldf<-seasonal.check.RPKM_uniq[ which(seasonal.check.RPKM_uniq$pval >= 0.05), ][,1]

#Now we need to find the common elements between seaonal and the quartiles and seasonal and the non quartiles. 
#A: high RPKMS 6713 elements -> split them into seasonal and not seasonal 
high_RPKMS_seas<-intersect(high_RPKMS,highly_seasonaldf) #5022
high_RPKMS_noseas<-intersect(high_RPKMS,non_seasonaldf) #1691 -> total 6713

#B) medium RPKMS 13425
med_RPKMS_seas<-intersect(medium_RPKMS,highly_seasonaldf) #11781
med_RPKMS_noseas<-intersect(medium_RPKMS,non_seasonaldf) #1644 -> total 13425

#C) low RPKMs 6713
low_RPKMS_seas<-intersect(low_RPKMS,highly_seasonaldf) #5861
low_RPKMS_noseas<-intersect(low_RPKMS,non_seasonaldf) #852 -> total 6713

#Create the PHYLOSEQs to plot these ones:
high_RPKMS_seas_subset <- subset(otu_table(PHYSEQ), rownames(otu_table(PHYSEQ)) %in% high_RPKMS_seas)
high_RPKMS_seas_physeq <- merge_phyloseq(high_RPKMS_seas_subset , sample_data(PHYSEQ))

high_RPKMS_noseas_subset <- subset(otu_table(PHYSEQ), rownames(otu_table(PHYSEQ)) %in% high_RPKMS_noseas)
high_RPKMS_noseas_physeq <- merge_phyloseq(high_RPKMS_noseas_subset , sample_data(PHYSEQ))

med_RPKMS_seas_subset <- subset(otu_table(PHYSEQ), rownames(otu_table(PHYSEQ)) %in% med_RPKMS_seas)
med_RPKMS_seas_physeq <- merge_phyloseq(high_RPKMS_seas_subset , sample_data(PHYSEQ))

med_RPKMS_noseas_subset <- subset(otu_table(PHYSEQ), rownames(otu_table(PHYSEQ)) %in% med_RPKMS_noseas)
med_RPKMS_noseas_physeq <- merge_phyloseq(med_RPKMS_noseas_subset , sample_data(PHYSEQ))

low_RPKMS_seas_subset <- subset(otu_table(PHYSEQ), rownames(otu_table(PHYSEQ)) %in% low_RPKMS_seas)
low_RPKMS_seas_physeq <- merge_phyloseq(low_RPKMS_seas_subset , sample_data(PHYSEQ))

low_RPKMS_noseas_subset <- subset(otu_table(PHYSEQ), rownames(otu_table(PHYSEQ)) %in% low_RPKMS_noseas)
low_RPKMS_noseas_physeq <- merge_phyloseq(low_RPKMS_noseas_subset , sample_data(PHYSEQ))

ALL_RPKMS_seas_subset <- subset(otu_table(PHYSEQ), rownames(otu_table(PHYSEQ)) %in% highly_seasonaldf) #22664
ALL_RPKMS_seas_physeq <- merge_phyloseq(ALL_RPKMS_seas_subset , sample_data(PHYSEQ))

ALL_RPKMS_noseas_subset <- subset(otu_table(PHYSEQ), rownames(otu_table(PHYSEQ)) %in% non_seasonaldf) #4187
ALL_RPKMS_noseas_physeq <- merge_phyloseq(ALL_RPKMS_noseas_subset , sample_data(PHYSEQ))


###ALL###
#Generate a Bray-Curtis decay plot function: Parameter to pass in the function is PHYSEQ.
#FUNCTION 3

function_plot_BC_decay <- function(PHYSEQ_OBJECT){
PHYSEQ_bray <- phyloseq::distance(physeq = PHYSEQ_OBJECT, method = "bray")
df_BC <- melt(as.matrix(PHYSEQ_bray), varnames = c("row_BC", "col_BC"))
names(df_BC)[3] <- "Bray_Curtis"
rownames(df_BC)<-paste(df_BC$row, df_BC$col, sep="_")
#Generate the days distances between samples to match bray Curtis
env_sr_vir<-data.frame(sample_data(PHYSEQ_OBJECT))
env_sr_vir$rnames<-rownames(env_sr_vir)
#Create a Matrix of days between samples
days_s2010 <- as.vector(difftime(env_sr_vir$Date_f,as.Date("2018-11-20","%Y-%m-%d"),units="days")) #I fixed 2018-11-20 as it was the first sample collected
dist_days <- as.matrix(dist(days_s2010,diag=TRUE,upper=TRUE))
rownames(dist_days) <- env_sr_vir$rnames; colnames(dist_days) <- env_sr_vir$rnames

dist_days[upper.tri(dist_days)] <- NA
dist_days <- melt(as.matrix(dist_days), varnames = c("row_day", "col_day"))
names(dist_days)[3] <- "days"
rownames(dist_days)<-paste(dist_days$row, dist_days$col, sep="_")

#merge the 2 generated dfs
days_BC<-merge(dist_days, df_BC, by=0, all=TRUE)
days_BC_lowt<-na.omit(days_BC) #Remove NAs
days_BC_lowt_final<-days_BC_lowt[,c(2:7)]
days_BC_lowt_final_no0<-days_BC_lowt_final[days_BC_lowt_final$days != 0, ]
#add a column with consecutive numbers so we can label 
days_BC_lowt_final_no0$comparison <- 1:nrow(days_BC_lowt_final_no0)

##Group by days and get the mean of each day category
days_BC_lowt_final_no0_aggr<-aggregate(days_BC_lowt_final_no0$Bray_Curtis, list(days_BC_lowt_final_no0$days), FUN=mean)
colnames(days_BC_lowt_final_no0_aggr)<-c("days","Bray_Curtis") #re attach the correct names

#In the following IF, we are going to plot the harmonic linear regression sinusoidal line if it has a significant amplitude, if not, just the scattered points
#pspline_fit <- lm(Bray_Curtis ~ mSpline(x = days, df = 4, periodic = TRUE, Boundary.knots = c(0,366)), data = days_BC_lowt_final_no0_aggr)
#df_WEC <- cbind(days_BC_lowt_final_no0_aggr, as.data.frame(predict(pspline_fit, interval = "prediction")))

#create the linear model
data_model<-days_BC_lowt_final_no0[,c("days","Bray_Curtis")]

model_BC<-lm(data_model$Bray_Curtis ~ sin(2 * pi * data_model$days / 365.25) + cos(2 * pi * data_model$days / 365.25))

predicted_df <- data.frame(days =data_model$days, pred_BC = predict(model_BC, data_model))

#calculate the stats
amplitude_estimate_BC <- unname(coef(model_BC)[2])
sin_coef <- unname(model_BC$coefficients[grepl("sin", names(model_BC$coefficients))])
cos_coef <- unname(model_BC$coefficients[grepl("cos", names(model_BC$coefficients))])
phase_estimate_BC <- atan2(sin_coef, cos_coef)
angular_frequency_BC <- 2 * pi / 365
amplitude_pvalue_BC <- summary(model_BC)$coef[2, 4] #print this
phase_pvalue_BC <- summary(model_BC)$coef[3, 4] #print this

total_variance <- sum((data_model$Bray_Curtis - mean(data_model$Bray_Curtis))^2)
explained_variance <- sum((predicted_df$pred_BC - mean(data_model$Bray_Curtis))^2)
residual_variance <- sum((data_model$Bray_Curtis - predicted_df$pred_BC)^2)
r_squared <- explained_variance / total_variance

#Evaluate and plot based wehther there was a significant amplitude or not
if (amplitude_pvalue_BC<0.05){
WEC_BC_stat<-ggscatter(days_BC_lowt_final_no0, x = "days", y = "Bray_Curtis",color = "black", shape = 21, size = 1.2)+ scale_x_continuous(breaks=seq(0, 1000, by = 365.25),limits=(c(0,1000)),labels=c("0","1","2")) + 
  geom_line(size=1.25,color='blue',data = predicted_df, aes(x=days, y=pred_BC))+
  xlab("time between samples (years)")+ylab("Bray-Curtis dissimilarity")+ylim(c(0,1))+
  theme(panel.background = element_blank(),legend.position = "none",strip.background = element_blank(),axis.text.y=element_text(colour = 'black',size=14), axis.text.x=element_text(colour = 'black',size=14),axis.line = element_line(colour = "black"),axis.title.x=element_text(colour = 'black',size=14),axis.title.y=element_text(colour = 'black',size=12) )+
  annotate("text", x = 150, y = .3, label = paste0("R-squared = ", round(r_squared, 2)), hjust = 0, vjust = 0, size = 5)+
  annotate("text", x = 150, y = .15, label = paste0("amplitude p-value = ", format(round(amplitude_pvalue_BC, 5), nsmall = 5)), hjust = 0, vjust = 0, size = 5)
}else{
 WEC_BC_stat<-ggscatter(days_BC_lowt_final_no0, x = "days", y = "Bray_Curtis",color = "black", shape = 21, size = 1.2)+ scale_x_continuous(breaks=seq(0, 1000, by = 365.25),limits=(c(0,1000)),labels=c("0","1","2")) + 
  #geom_line(size=1.25,color='blue',data = predicted_df, aes(x=days, y=pred_BC))+
  xlab("time between samples (years)")+ylab("Bray-Curtis dissimilarity")+ylim(c(0,1))+
  theme(panel.background = element_blank(),legend.position = "none",strip.background = element_blank(),axis.text.y=element_text(colour = 'black',size=14), axis.text.x=element_text(colour = 'black',size=14),axis.line = element_line(colour = "black"),axis.title.x=element_text(colour = 'black',size=14),axis.title.y=element_text(colour = 'black',size=12) )+
  annotate("text", x = 150, y = .3, label = paste0("R-squared = ", round(r_squared, 2)), hjust = 0, vjust = 0, size = 5)+
  annotate("text", x = 150, y = .15, label = paste0("amplitude p-value = ", format(round(amplitude_pvalue_BC, 4), nsmall = 4)), hjust = 0, vjust = 0, size = 5)
}
WEC_BC_stat #Final PLOT to retrieve
}

plot_ALL<-function_plot_BC_decay(PHYSEQ)
plot_high<-function_plot_BC_decay(high_physeq)
plot_med<-function_plot_BC_decay(med_physeq)
plot_low<-function_plot_BC_decay(low_physeq)
plot_ALLseas<-function_plot_BC_decay(ALL_RPKMS_seas_physeq)
plot_ALLnoseas<-function_plot_BC_decay(ALL_RPKMS_noseas_physeq)
plot_highseas<-function_plot_BC_decay(high_RPKMS_seas_physeq)
plot_highnoseas<-function_plot_BC_decay(high_RPKMS_noseas_physeq)
plot_medseas<-function_plot_BC_decay(med_RPKMS_seas_physeq )
plot_mednoseas<-function_plot_BC_decay(med_RPKMS_noseas_physeq)
plot_lowseas<-function_plot_BC_decay(low_RPKMS_seas_physeq)
plot_lownoseas<-function_plot_BC_decay(low_RPKMS_noseas_physeq)

#svg("FigS1_3090V2.svg", width=16,height=12)
plot_grid(plot_ALL,plot_high,plot_med,plot_low,plot_ALLseas,plot_highseas,plot_medseas,plot_lowseas,plot_ALLnoseas,plot_highnoseas,plot_mednoseas,plot_lownoseas, align="h", nrow=3)
#dev.off()

#svg("FigS2_extendedV2.svg", width=16,height=12) -> USE OTU in line 49 to construct the phyloseq object using the 26k dataset and not the high-quality 3090
plot_grid(plot_ALL,plot_high,plot_med,plot_low,plot_ALLseas,plot_highseas,plot_medseas,plot_lowseas,plot_ALLnoseas,plot_highnoseas,plot_mednoseas,plot_lownoseas, align="h", nrow=3)
#dev.off()

#Figure 02 PANEL A
#svg("New_Fig2_panelA.svg", height=5, width=10)
plot_grid(plot_ALL)
#dev.off()
