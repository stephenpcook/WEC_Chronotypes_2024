#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

countwec <- read.table(args[1],header=T, row.names=1, check.names=F, sep ="\t")
VIR_cl_env <- read.table(args[2],header=T, row.names=1, check.names=F, sep ="\t")
num_iterations <- as.integer(args[3]) #Bootstrap replicates
size_subsample <- as.integer(args[4]) #specify the size of the subsampling datasets
output_path <- args[5] # Output path for writing results
#start_year <- as.integer(args[6]) #TSconversor parameters
#start_month <- as.integer(args[7]) #TSconversor parameters
#end_year <- as.integer(args[8]) #TSconversor parameters
#end_month <- as.integer(args[9]) #TSconversor parameters
#frequency_all <- as.integer(args[10]) #TSconversor parameters

#TSconversor parameeters

# Check that output path ends with a slash, if not add it
if (!grepl("/$", output_path)) {
    output_path <- paste0(output_path, "/")
}

# Check if the output directory exists, if not, create it
if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE)
}

rm(args)

### Load libraries ### 
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
library(GeneCycle) #this library performs the fisher.g.test
library(zoo) #this library performs the interpolation
library(ggrepel)
library(KneeArrower)

###DEFINE FUNCTIONS###

#FUNCTION 1: convert our datasets into a time-series object
#Parameters are inputted by the user depending on the start/end and frequency of its datasets
#This could be adapted to any type of frequency: hours, days, months, years, etc. 

 # ts_conversor<- function(df, s=c(start_year, start_month), e=c(end_year,end_month), f=frequency_all ){
 #ts(df$value, start = s, end=e, frequency = f) 
 #}

#FUNCTION 2 Seasonality check 
#check_seasonality_fisher <- function(ls){
  #ls %>% 
    #map(~ts_conversor(.x)) %>% 
    #map(~data.frame(per = periodogram(.)$freq[
      #which.max(periodogram(.)$spec)],
      #whole = periodogram(.),
      #pval = fisher.g.test(.))) %>% 
  #bind_rows( .id = "variable") 
#}

#FUNCTION 3 extract the centroids for the number of subsamples the user inputs 
#USAGE: data.frame (usr input), how many subsets of 100 contig time-series user want to use to estimate multiple centroids (usr input), distance type (tsclust options) and centroid type (tsclust options)

centroid_generator<-function(data,number_subsets,usr_dist_type,usr_centroid_type,size) {
  centroid_list = vector("list",length(number_subsets))
  available_data <- data #initializie available data
  data.sub<-data[FALSE,] #initialize data.sub with the same structure that "data"
   for (i in 1:number_subsets) {
    available_data <- anti_join(available_data,data.sub) #first cycle is full dataset compared to an empty one, the next cycle available_data will be a new df without the first cycle extract data nd like thut succesively. 
    data.sub<- sample_n(available_data,size) #THIS IS BEING SAMPLE WITHOUT REPLECEMENT IN THIS VERSION of the script, similar as generate_reprs_contigs_wo_replacement.R
    data_rpkms100.pam <- tsclust(data.sub , type="partitional", k=3L:(size-1), distance=usr_dist_type, centroid=usr_centroid_type)
    data_clust_stats100<-sapply(data_rpkms100.pam , cvi, type = "internal")
    maxSil<-max(data_clust_stats100[1,])
    estK<-(which(data_clust_stats100==maxSil,arr.ind=TRUE)[,2])+2
    bestdefined_clustering.pam <- tsclust(data.sub, type="partitional", k=estK, distance=usr_dist_type, clustering=usr_centroid_type)
    bestdefclust.centroids<-attr(bestdefined_clustering.pam@centroids,"series_id")
    centroid_list[[i]]<-c(rownames(data.sub[bestdefclust.centroids,]), paste(rownames(data.sub),"cluster",bestdefined_clustering.pam@cluster,sep = "_") )
  }
  centroid_list
}

#FUNCTION 4: This function was created for the specific dataset analyzed here
#Users are welcome to modify it depending of their own requierements. 

function_plotsTS <- function(plots,column_name,DataFrame){
cluster<- DataFrame[DataFrame[[column_name]]==plots,]$obs
#Here subset the zscores by rownames (this is different than the other)
cluster_df_zscore<-monthly_all_rpkms.normz.forplot.melt[monthly_all_rpkms.normz.forplot.melt$contig %in% cluster, ]
#PLOT
ggplot(cluster_df_zscore,aes(x = variable, y = value, group=contig)) +
geom_line(size=.2,aes( x=variable,y=value)) +
scale_x_discrete(expand = c(0.01,0))+
geom_vline(xintercept = c(2,14,26), color = "black")+theme_bw()+
guides(color=FALSE)+theme(axis.text.x = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),axis.ticks.y=element_blank(),axis.ticks.x=element_blank())
  }

#FUNCTION 5: This function estimates de euclidean distances
euclidean_distance_func <- function(x, y) {
  sqrt(sum((x - y) ^ 2))
}

#FUNCTION 6 Within_MeanEuclDist:
function_Within_MeanEuclDist <- function(plots,column_name,DataFrame){
cluster<- DataFrame[DataFrame[[column_name]]==plots,]$obs
#Here subset the zscores 
cluster_df_zscore<-monthly_all_rpkms.normz.forplot.melt[monthly_all_rpkms.normz.forplot.melt$contig %in% cluster, ]
#Unmelt
unmelted_zscore<-t(data.frame(cluster_df_zscore  %>% pivot_wider(names_from = contig, values_from = value)))
#remove the first non numeric line
colnames(unmelted_zscore) <- unmelted_zscore[1,]
# Remove the first row from the data frame (optional)
unmelted_zscore <- data.frame(unmelted_zscore[-1, ])
#turn the numbers into numeric not factor
unmelted_zscore <-unmelted_zscore %>% mutate_if(is.character, as.numeric)
pairwise_distances <- dist(as.matrix(unmelted_zscore), method = euclidean_distance_func)
mean_distance <- mean(pairwise_distances)
}

#Set the directory to the output inputted by the user
setwd(output_path) #I did this to output the reporting curves and heatmaps into the output dir - comment if this is not what you want.

countwect<-t(countwec) #transpose Time-series data frame
VIR_cl_env$sample<-rownames(VIR_cl_env)
VIR_cl_envdate<-VIR_cl_env[ , c("sample","Date_wavelet")]
VIR_cl_envdate$sample<-paste0("X",1:nrow(VIR_cl_envdate),".",VIR_cl_envdate$sample)
VIR_cl_envdate$Date_wavelet<-as.Date(VIR_cl_envdate$Date_wavelet,"%d/%m/%Y")

monthly_RPKM<-merge(VIR_cl_envdate,countwect,by=0, all = TRUE) ######I need to interpolate this
monthly_RPKM$Date_wavelet<-as.Date(monthly_RPKM$Date_wavelet,"%d/%m/%Y")
monthly_RPKM_sorted<-(monthly_RPKM[order(monthly_RPKM$Date_wavelet),]) #dim 32 26854

#Interpolation of the numerical columns. Im skipping 1:3 becuase these are other alphanumerical strings (Row.names-sample-Date_wavelet)
monthly_RPKM_sort_trim_int<- data.frame(monthly_RPKM_sorted[1:3], na.approx(monthly_RPKM_sorted[4:length(monthly_RPKM_sorted)])) #na.approx fn -> interpolation
monthly_RPKM_sort_trim_int_cl<-monthly_RPKM_sort_trim_int[,c(3:length(monthly_RPKM_sorted))]

################################################
#############SEASONALITY CHECK##################

#monthly_RPKM_sort_trim_int_cl_melt<-melt(monthly_RPKM_sort_trim_int_cl, id.vars=c("Date_wavelet")) #this is for the seasonal check

#Run the seasonal check
#seasonal.check.RPKM <-monthly_RPKM_sort_trim_int_cl_melt %>% 
  #split(.$variable,drop = T) %>% 
  #check_seasonality_fisher()

#select unique rows
#seasonal.check.RPKM_uniq<-seasonal.check.RPKM[!duplicated(seasonal.check.RPKM$variable),]

#split significant seasonal and non-seasonal (p-vale can be adjusted, we use 0.05)
#highly_seasonaldf<-seasonal.check.RPKM_uniq[ which(seasonal.check.RPKM_uniq$pval < 0.05), ][,1]
#non_seasonaldf<-seasonal.check.RPKM_uniq[ which(seasonal.check.RPKM_uniq$pval >= 0.05), ][,1]

#Subset these list from monthly_RPKM_sort_trim_int_cl_melt into two different data frames 
#Seasonal_significant_DF_melt <- monthly_RPKM_sort_trim_int_cl_melt[monthly_RPKM_sort_trim_int_cl_melt$variable %in% highly_seasonaldf, ]
#not_significant_DF_melt <- monthly_RPKM_sort_trim_int_cl_melt[monthly_RPKM_sort_trim_int_cl_melt$variable %in% non_seasonaldf, ]

#reshape the data frames
#Seasonal_significant_DF<-dcast(data = Seasonal_significant_DF_melt,formula = Date_wavelet~variable,value.var = "value") 
#not_significant_DF<-dcast(data = not_significant_DF_melt,formula = Date_wavelet~variable,value.var = "value") 

#Remove the first column (date wavelet)
#nums <- Seasonal_significant_DF[,2:length(Seasonal_significant_DF)] #Remove the wavelet column (1)
#nums_nots<-not_significant_DF[,2:length(not_significant_DF)] #Remove the wavelet column (1)

#transpose both - These are the objects that can also be inputted if user has already checked the seasonality of its data
#SignSeasonal_rpkms_subset<-t(nums) #contigs are rows #write.table(tnums, "/PATH/ctgs_interpolatedRPKMS.seasonal22k.tab", quote=FALSE, sep = "\t")
#notsign_rpkms_subset<-t(nums_nots) #write.table(tnums_nots, "/PATH/ctgs_interpolatedRPKMS.notseasonal4k.tab", quote=FALSE, sep = "\t")

#Rename the column headers -> automatize this
#colnames(SignSeasonal_rpkms_subset) <-VIR_cl_envdate$sample 
#colnames(notsign_rpkms_subset)<- VIR_cl_envdate$sample

#Sanity check missing values (this may be removed)
#summary(is.na(SignSeasonal_rpkms_subset)) (this may be removed)
#summary(is.na(notsign_rpkms_subset)) (this may be removed)

#z-score normalisation (this may be removed)
#SignSeasonal_rpkms.normz <- data.frame(zscore(SignSeasonal_rpkms_subset)) (this may be removed)
#notsign_rpkms.normz<-data.frame(zscore(notsign_rpkms_subset)) (this may be removed)

###############END OF SEASONALITY CHECK#######################

##work in progress#
#In the new version the seasonality check will be implemented. This means the user can indicate in the input paremeter 
#whether the time-series should be evaluated separately (seasonal vs non seasonal)
#In the current version, it is assumed the user has run XX code first to output the seasonal and not seasonal rows of the time-series
#then ChronoClustR.R is run with each of this output files. (User can used the same env file)

#####Reformat the monthly RPKM table so we can start from there (non splitting seasonal / non seasonal)

monthly_RPKM_sort_trim_int_cl$Date_wavelet<-as.Date(monthly_RPKM_sort_trim_int_cl$Date_wavelet,"%d%m/%Y")
monthly_RPKM_sort_trim_int_cl_merged<-merge(monthly_RPKM_sort_trim_int_cl,VIR_cl_envdate, by='Date_wavelet', all=TRUE)
rownames(monthly_RPKM_sort_trim_int_cl_merged)<-monthly_RPKM_sort_trim_int_cl_merged$sample
monthly_RPKM_sort_trim_int_cl_merged_numbers <- monthly_RPKM_sort_trim_int_cl_merged[-c(1, ncol(monthly_RPKM_sort_trim_int_cl_merged))]
monthly_RPKM_sort_trim_int_cl_tr<-t(monthly_RPKM_sort_trim_int_cl_merged_numbers)

############# STEP 1 z-score normalization ##########
monthly_all_rpkms.normz<-data.frame(zscore(monthly_RPKM_sort_trim_int_cl_tr))
monthly_all_rpkms.normz.forplot<-monthly_all_rpkms.normz
monthly_all_rpkms.normz.forplot$contig<-rownames(monthly_all_rpkms.normz.forplot)
monthly_all_rpkms.normz.forplot.melt<-melt(monthly_all_rpkms.normz.forplot, id.vars=c("contig"))


############# STEP 2 Run the function to extract the centroids of multiple 100 contigs subsets##########
#I created a function that subsets 100 contig time-series, analyze them, cluster them, and output the centroids of the unsupervised clustering as a list of names. 
#In a further step we will cluster them to see whether we can reach a saturation of centroids representation.
#I will use this function independently for each dataset (seasonal - not seasonal), subseting 100 each (for now), and bootstrapping to see the frequency of centroids we retrieve. 
#Unsupervised clustering will be sensitive to the data composition. That's why we need to iterate many times to retrieve a representative collection of patterns.
#Otherwise, unsupervised clustering will try to fit potentially different patterns into 1 to fulffil the parameters imposed (minimize distance within cluster and maximize it between clusters). 

#DEFAULT for:
#b) size_subsample = 74
datasets100 = as.integer(nrow(monthly_all_rpkms.normz)/size_subsample)

#declare an output list
FINAL_list <- list()

#num_iterations <- 10 #this can be inputted by the user

#####################
#Bootsrapping starts#
#####################

for (i in 1:num_iterations) {

#In the centroid generator function read the variable from std input
seasonal_centroids_10reps<-centroid_generator(monthly_all_rpkms.normz,datasets100,"Euclidean","pam",size_subsample) #After getting this, ask how do these centroids relate to each other?

#add name SubsetX to each element of the list to keep track of the 
NAME <- paste0("subset", 1:length(seasonal_centroids_10reps))
names(seasonal_centroids_10reps) <- NAME

#Extend the list saving its clustering representation. 
df_seasonal_centroids_10reps_tmp_tobesplited<-data.frame(ID = rep(names(seasonal_centroids_10reps), sapply(seasonal_centroids_10reps, length)), Obs = unlist(seasonal_centroids_10reps))

#Split the above data frame into representative centroids and full subsetted n groups of 100
df_seasonal_centroids_10reps<-subset(df_seasonal_centroids_10reps_tmp_tobesplited, !grepl('_cluster_', Obs)) #This set might not have pam.centroids. I need to unique here so the program don't crash

#select unique rows of df_seasonal_centroids_10reps
df_seasonal_centroids_10reps_uniq<-df_seasonal_centroids_10reps[!duplicated(df_seasonal_centroids_10reps$Obs),]

#######The object below will be used generate a potential training dataset: to be used in a following section [CHECK down below continuation]
df_seasonal_subsets_10reps<-subset(df_seasonal_centroids_10reps_tmp_tobesplited, grepl('_cluster_', Obs)) #this will includ the representative centroids with the subset number in one column and its cluster belongin embedded within the name in this form _cluster_X

#Subset, add the column of the subset number, melt and plot
df_seasonal_centroids_10reps_subset<-monthly_all_rpkms.normz[rownames(monthly_all_rpkms.normz) %in% df_seasonal_centroids_10reps_uniq$Obs, ] #centroids found total 124, but only 77 unique.
df_seasonal_centroids_10reps_subset$Obs<-rownames(df_seasonal_centroids_10reps_subset)

merged10reps<-merge(df_seasonal_centroids_10reps_subset, df_seasonal_centroids_10reps_uniq, by='Obs', all=TRUE)

##PENDING: EVALUATE WHY WE GET SO MANY REPETITIVE CENTROIDS FROM 10 SUBSETS.

merged10reps_melt<-melt(merged10reps, id.vars=c("Obs","ID")) #Table was collapsed from here, check this stuff in the next round of tests

#Change the names of merged10reps_melt
merged10reps_melt$Obs<- gsub(".viral.fraction__viral.spades.short__as.yet.unknown__", "", merged10reps_melt$Obs)
merged10reps_melt$Obs<- gsub("..full", "", merged10reps_melt$Obs)

first_clustrds<-ggplot(merged10reps_melt,aes(x = variable, y = value, group=Obs)) +
geom_line(size=1,aes( x=variable,y=value, color=ID)) +
scale_x_discrete(expand = c(0.01,0))+
#stat_smooth(size=1,aes( x=variable,y=value, color=ID), method = lm, formula = y ~ poly(x, 6), se = FALSE)+
geom_vline(xintercept = c(3,15,27), color = "black")+theme_bw()+
facet_wrap(ID~Obs, ncol = 10, strip.position="left") +
labs(x ="date", y = "z-score normalized RPKM")+
guides(color=FALSE)+theme(axis.text.x = element_blank())

############# STEP 3 cluster the centroids ##########
#Evaluate 3 methods of clustering the centroids#
#based on the first derivative and the approximation to lower the within distances of the clusters, generate an "optimal" Number of centroids representing the random subsetted 

#Cluster the centroids and evaluate
#merged10reps_melt<-readRDS("/Users/lb808/Documents/MyDrafts/InProgress/WEC_Viromes_ShortReads/ModelDesign/merged10reps_melt.rds")
#reshape if reading from RDS
#merged10reps = dcast(merged10reps_melt,  Obs~variable)

####Cluster the centroids of the multiple subsamplings

#generate a data frame to be inputted to an unsupervised clustering
rownames(merged10reps)<-merged10reps$Obs
merged10reps_toinput = merged10reps[,!(names(merged10reps)== c("Obs","ID"))] #drop off the contig name and the ID columns

#Evaluate stats generated from clustering 3 to 32 #I need to atomatize this from 3 to X (how to determine X). I mean, it does not really matter if we are going to iterate over what subsets we do this
merged10reps.pam <- tsclust(merged10reps_toinput , type="partitional", k=3L:(nrow(merged10reps_toinput)-1), distance="Euclidean", centroid="pam") #si cambias este rango (k:) se afecta

#evaluate the clusters depending on their number. cvi for each cluster
merged10reps.pam_clust_stats100<-sapply(merged10reps.pam, cvi, type = "internal")

#Get the highest Silhoute score
maxSil_merged10reps<-max(merged10reps.pam_clust_stats100[1,])

#Exctract the k kluster that maximizes the Silhoute
K_centroids<-(which(merged10reps.pam_clust_stats100==maxSil_merged10reps,arr.ind=TRUE)[,2])+2 #exctract the col number where the max value is located this equals at the cluster-2 (we started the eval at 3)

#cluster the non seasonal into K_seasonl value 
#This gets way better if we use the hierarchical / elbow methods down below. I don't really care to use the best point in the elbow. But to really minimize distance within the cluster
clust.seasonal.pam <- tsclust(merged10reps_toinput, type="partitional", k=K_centroids, distance="Euclidean", clustering="pam") #Pass the variable of the number of clusters that maximize Silhoute (above)

#PLOT clustering of centroids (12 clusters)
centroids.ctg_clusters<-cbind(merged10reps_toinput[,0], cluster = clust.seasonal.pam@cluster)
merged10reps_cls<-merged10reps_toinput
merged10reps_cls$cluster<-centroids.ctg_clusters[rownames(merged10reps_toinput),]
merged10reps_cls$cluster<- gsub(" ","",c(paste("cluster",merged10reps_cls$cluster)))
merged10reps_cls$contig<-rownames(merged10reps_cls)
merged10reps_cls_melt<-melt(merged10reps_cls, id.vars=c("contig","cluster"))


centroids10_clustered<- merged10reps_cls_melt %>% ggplot(aes(x = variable, y = value,group=contig, label=contig)) +
geom_line(linewidth=.6,aes( x=variable,y=value, color=cluster)) +
scale_x_discrete(expand = c(0.01,0))+
#stat_smooth(size=1,aes( x=variable,y=value, color=ID), method = lm, formula = y ~ poly(x, 6), se = FALSE)+
geom_vline(xintercept = c(3,15,27), color = "black")+theme_bw()+
# Labels defined here
#geom_label_repel(
#    data = merged10reps_cls_melt %>% group_by(contig) %>% filter(value == max(value)), 
#    aes(label = contig),color = "black",box.padding = 0.5, max.overlaps = Inf, size=1.9)+
facet_wrap(~cluster, ncol = 6, strip.position="left") +
labs(x ="date", y = "z-score normalized RPKM")+
guides(color=FALSE)+theme(axis.text.x = element_blank())

#svg("/Users/lb808/Documents/MyDrafts/InProgress/WEC_Viromes_ShortReads/ModelDesign/merged10reps_centroids_cls.svg", width=19,height=14)
#centroids10_clustered
#dev.off()

##### STEP 3.1 HeatMap and dendrogrom of the centroids #####

#This clustering is not bad at all, but I want it to identify slight "monthly" variation dynamics. 
#In this case I don't care about maximizing the distance between group distance, but make sure that representative clusters of centroids are tightly clustered. 
#Therefore I need to identify a threshold that minimizes the distance as much as I can between clusters

dist_ts <- TSclust::diss(SERIES = as.matrix(merged10reps_toinput), METHOD = "EUCL") #convert to matrix so it considers the centroids not the "samples/dates"
#print(dist_ts)

euc_dist <- function(x){TSclust::diss(SERIES = (x), METHOD = "EUCL")}

#This is not saving the heatmap image
heatmap_centroids<-merged10reps_toinput %>%
  t() %>% as.matrix() %>%
  gplots::heatmap.2 (
    # dendrogram control
    distfun = euc_dist,
    hclustfun = hclust,
    dendrogram = "column",
    Rowv = FALSE,
    tracecol=NA,
    labRow = FALSE
  )

#Use the distance that maximize the identityt between clusters -> Agglomerative clustering
#choose 3 main parameters: 
#a) Distance measure (quantify dissimilarity), 
#b) Prototype (summarizes characteristics of all series in a cluster), 
#c) Cluster algorithm (most common partitional or hierarchical). 
#To finish, evaluate results: cluster validity indices (CVI)

#Corroborate this by two methods 

#######Dendrogram#######
# the ward.D method tries to minimize the within cluster variance
dendrogram = hclust(d = dist(merged10reps_toinput, method = 'euclidean'), method = 'ward.D')

#Create a function that cuts the tree in different distances and then claculate the avarage distance within the cluster select the ona that maximizes 
#k-means, with several starting values, and the gap statistic (mirror) to determine the number of clusters that minimize the within-SS (Sum of squares)

coph_dists_dendro<-cophenetic(dendrogram) #extract cophenetic distances
max_coph_dists_dendro<-max(coph_dists_dendro) #in this case 66.10528, cutting here will output 1 cluster
min_coph_dists_dendro<-min(coph_dists_dendro) #in this case 0.6750531, cutting here will output n-1 clusters. 

#Below at some point: For loop to find the distance in the dendrogram that splits it in the closest number of clusters revelead by the derivative analysis

##### STEP 3.2 Kmeans of the centroids #####
set.seed(123)

#function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(merged10reps_toinput, k, nstart = 10, iter.max = 30)$tot.withinss #I want to minimize ALL the withinss not the total. But first we need to compute the total so we can remove the clusters that optimize what i don;t want
}

# Compute and plot wss for k = 1 to k = n-1
k.values <- 1:(nrow(merged10reps_toinput)-1)

# extract wss k.values
wss_values <- map_dbl(k.values, wss)

#svg("/Users/lb808/Documents/MyDrafts/InProgress/WEC_Viromes_ShortReads/ModelDesign/kvalues_wssvalues.svg", width=12)
kvalues_wssvalues<-plot(k.values, wss_values,
       type="b", pch = 19, frame = FALSE, 
       xlab="Number of clusters K",
       ylab="Total within-clusters sum of squares") #I don't want the total, I'd like to minimize the within cluster distance of all clusters. -> evaluate distributions of clusters after the elbow change
#dev.off()

####1st derivative (slope of the tangent line) cutoff: "This method finds the point along the curve where the slope is a given fraction of the maximum"

#1)PLOT and get the smooth curve
#geom_smooth() predicts for a internal sequence generated for the entire range of x, while loess in the example predicts only for the unique x values existing in the dataset (they have different results)
  
#Generate the smooth curve

#lw1 <- loess(wss_values ~ k.values) #is not a lineal curve
#plot(k.values, wss_values,type="l") CONTROL ONLY
lw2 <- loess(wss_values~ log(k.values)) #the k-means clusters vs distances almost all time follow a log curve, better use log to fit the model 
xl <- seq(min(k.values),max(k.values), (max(k.values) - min(k.values))/1000) #generate 1k datapoints to feed the formula
out = predict(lw2,log(xl)) #generate the wss_values that smooths the line under the loess formula! 
#lines(xl, out, col='red', lwd=2)

#retrieve the point at which the slope of tangent line (1st derivative) is 1/10 of its maximum value. Good way to get the elbow!
cutoff.point <- findCutoff(xl, out, method="first", 0.1)

#retrieve the closest integer on the x axis which is the number of clusters that "maximize distance between clusters and minimize distance within". 
#I want to evaluate the distance distribution of clusters below this point where the sum of the distances don't change a lot. Basically find the cluster configuration that minimizes the distancs of its clusters, regardles of the between cluster distance. 

rounded.cutoff.point<-round(cutoff.point$x) #round to the closest integer in X (number of clusters variable)

#Plot our cutoff point in the smooth version. CONTROL ONLY
#plot(xl, out, pch=20, col="gray")
#points(cutoff.point, col="red", cex=3, pch=20)

###withinss values from everything above rounded.cutoff.point
wss_all <- function(k) {
  kmeans(merged10reps_toinput, k, nstart = 10,iter.max = 30 )$withinss #I want to minimize ALL the withinss not the total. But first we need to compute the total so we can remove the clusters that optimize what i don;t want
}

k.values.aboverounded<-rounded.cutoff.point:(nrow(merged10reps_toinput)-1)

withinss.list<-sapply(k.values.aboverounded, wss_all) #we lost the minimim value here and it reseets, because we removed everything above the first derivative.

#Function to generate a dataframe with all the withinss values from the k-means clusters greater than the rounded cutoff (1/10th of the max first derivate)
summarystat<- function(x) {
  z1 <- mean(x)
  z2 <- median(x)
  z3 <- sd(x)
  return(list(mean=z1, median=z2, sd=z3))
}

list.stats<-lapply(withinss.list,summarystat)

#FIRST FILTER: median can not be 0 #We don't want a bunch of singletons.
list.stats.medianfilt<-Filter(function(x) x$median !=0, list.stats)

#SECOND FILTER, check the knee of the curve (1st derivative of the means) but get the one with the lowest standard deviation.
df.stats.medianfilt<-data.frame(matrix(unlist(list.stats.medianfilt), nrow=length(list.stats.medianfilt), byrow=TRUE))
colnames(df.stats.medianfilt) <- c("mean", "median", "sd")

means_dists_kclusters<-df.stats.medianfilt$mean
kclusters_number<-c(1:length(df.stats.medianfilt$mean))

#plot(kclusters_number,df.stats.medianfilt$mean,  pch=20, col="black") CONTROL ONLY
lw_means <- loess(df.stats.medianfilt$mean~ (kclusters_number^2)) #the k-means clusters vs distances almost all time follow a log curve, better use log to fit the model 
xl_klusters.means <- seq(min(kclusters_number),max(kclusters_number), (max(kclusters_number) - min(kclusters_number))/1000) #generate 1k datapoints to feed the formula
out_means = predict(lw_means,xl_klusters.means) #generate the wss_values that smooths the line under the loess formula! 
#lines(xl_klusters.means, out_means, col='red', lwd=2) CONTROL ONLY

cutoff.point_means <- findCutoff(xl_klusters.means, out_means, method="first", 0.25) #this is cuadratic not log, therefore I used 1/4 of the max and not 1:10 of the max. 
plot(xl_klusters.means, out_means, pch=20, col="gray")
points(cutoff.point_means, col="red", cex=3, pch=20)

best_number_of_clusters<-round(cutoff.point_means$x) #13, therofore 21 clusters

#THIRD FILTER: Once we have a rounded integer indicating the cluster which better fits the knee of the means curve. Explore 1/4 of the total of clusters looking for the smallest sd

#Get the integer representing 1/4 of the total clusters 
#limit_for_minSD<-round(round(length(kclusters_number)/4)/2) #I divide between 2 so I can have an upper and inferior limit to explore.

#min_pos_subset_df<-rounded.cutoff.point.k.means-limit_for_minSD
#max_pos_subset_df<-rounded.cutoff.point.k.means+limit_for_minSD

#df.stats.medianfilt.subset<-df.stats.medianfilt[min_pos_subset_df:max_pos_subset_df,]

#best_number_of_clusters<-(as.integer(rownames(df.stats.medianfilt.subset[which.min(df.stats.medianfilt.subset$sd),])))+rounded.cutoff.point #we hhave been evaluating number of clusters w/o the first fraction above the ctoff point so we need to add this value to select the exacto No of kluster

#I need to benchmark this by repeating multiple times this function to get a relationship of number of subsets of 100 contigs vs best number of clusters

##### STEP 3.3 put together the classification results using the optimal solution we found in best clusters number #####
#Once we determined the best number of clusters, sse Kmeans, tsclust and hierarchical clustering and compare the clustering. This is a sanity check step
#dataframe should be available to the user, this is not something I'll evaluate exhastively. 

##(A)##
#Kmeans with best_number_of_clusters
kmeans_best_numberk<-kmeans(merged10reps_toinput, best_number_of_clusters, nstart = 10 ,iter.max = 30) #caluclate kmeans with 26 clusters
kmeans_best_number_contigcluster<-kmeans_best_numberk$cluster #extract the relationship of contigId to cluster membership
kmeans_best_number_contigcluster_df<-data.frame(kmeans_best_number_contigcluster) #as df
colnames(kmeans_best_number_contigcluster_df)<-"kmeans" #change the name of the column with the cluster membership

#get the pam centroids, z-scores and cluster membership and write it in a table. This will be the training dataset
kmeans_best_number_contigcluster_df$kmeans<-paste0("cluster_",kmeans_best_number_contigcluster_df$kmeans)
training_dataset<-rownames(kmeans_best_number_contigcluster_df)
training_dataset_zscore<-merged10reps_toinput[rownames(merged10reps_toinput) %in% training_dataset, ]

training_dataset_final<-merge(kmeans_best_number_contigcluster_df, training_dataset_zscore, by='row.names', all=TRUE)

##(B)##
#Dendrogram
####find the distance needed to cut the closest integer to best_number_of_clusters in the dendrogram

#Generate a sequence of integers from the min distance to the max distance
dists_to_evaluate<-seq(min(min_coph_dists_dendro),max(min_coph_dists_dendro:max_coph_dists_dendro),.01)

#for loop to find the distance that is closer to generate the best number of clusters variable
#for (val in dists_to_evaluate) {
#  klusters_tree<-length(unique(cutree(dendrogram, h=val)))
#  if(round(klusters_tree) == best_number_of_clusters) klusters_tree_list<- cutree(dendrogram, h=val)} #we will take the last "decimal value" that fulffills the integer match
#print(klusters_tree_list)

klusters_tree_list<-cutree(dendrogram, k=best_number_of_clusters)
dendro_best_number_contigcluster_df<-data.frame(klusters_tree_list)

##(C)##
#TSCLUST unsupervised PAM - EUCL
clust.best_number_contig.pam <- tsclust(merged10reps_toinput, type="partitional", k=best_number_of_clusters, distance="Euclidean", clustering="pam") 
unsup_best_number_contigcluster_df<-cbind(merged10reps_toinput[,0], unsup_cluster =clust.best_number_contig.pam@cluster)
unsup_best_number_contigcluster_df$Row.names<-rownames(unsup_best_number_contigcluster_df)

#############
#Merge A,B,C#
#############
tmp_merged_klusters<-merge(kmeans_best_number_contigcluster_df, dendro_best_number_contigcluster_df, by=0, all=TRUE) 
rownames(tmp_merged_klusters)<-tmp_merged_klusters$Row.names
tmp2_merged_klusters<-merge(tmp_merged_klusters,unsup_best_number_contigcluster_df, by="Row.names", all=TRUE)

#variable merged_klusters contains the contigs organized in clusters (optimal solution). 
#merged_klusters<-tmp2_merged_klusters[,2:length(tmp2_merged_klusters)] #remove the duplicated rownames column

#Final clustering step: Extract the K-means centroids from the kmeans solution for the best number of clusters
#This can be printed at the end so the usr can check how many of the contigs were assigned to the same cluster using this different methods.
#All methods will give a different solution (due to inherent features), however these should be comparable and closely related centroids should be placed together reardless of the method

kmeans_best_numberk_centers<-data.frame(kmeans_best_numberk$centers)
#plot: visualize the optimal solution = kmeans centroids
kmeans_best_numberk_centers$centers<- gsub(" ","",c(paste("center_cluster",rownames(kmeans_best_numberk_centers))))
kmeans_best_numberk_centers_melt<-melt(kmeans_best_numberk_centers, id.vars=c("centers"))

kmeans_best_numberk_centers_melt_plot<- kmeans_best_numberk_centers_melt %>% ggplot(aes(x = variable, y = value,group=centers, label=centers)) +
geom_line(linewidth=.6,aes( x=variable,y=value, color=centers)) +
scale_x_discrete(expand = c(0.01,0))+
#stat_smooth(size=1,aes( x=variable,y=value, color=ID), method = lm, formula = y ~ poly(x, 6), se = FALSE)+
geom_vline(xintercept = c(3,15,27), color = "black")+theme_bw()+
# Labels defined here
#geom_label_repel(
#    data = kmeans_best_numberk_centers_melt %>% group_by(centers) %>% filter(value == max(value)), 
#    aes(label = centers),color = "black",box.padding = 0.5, max.overlaps = Inf, size=1.9)+
facet_wrap(~centers, ncol = 6, strip.position="left") +
labs(x ="date", y = "z-score normalized RPKM")+
guides(color=FALSE)+theme(axis.text.x = element_blank())

#Report for benchamrking (one line that will be concatenated with the other runs)
df_report_benchmarck <- setNames(data.frame(t(c(datasets100, nrow(kmeans_best_numberk_centers)))),c("datasets_sampled", "optimal_no_clusters"))

now <- Sys.time() #time stamp used as row name 
rownames(df_report_benchmarck) <- now
#df_report_benchmarck

#first unsupervised clustering of centroids
#svg("unsup_centroids_cls.svg", width=19,height=14)
#centroids10_clustered
#dev.off()

#Final collection of k-means centroids
#svg("merged_reps_centroids_cls.svg", width=19,height=14)
#kmeans_best_numberk_centers_melt_plot
#dev.off()

###CONTINUATION full subsets to link with the centroid clustering
df_seasonal_subsets_10reps$cluster<- sub(".*cluster_", "", df_seasonal_subsets_10reps$Obs) #cluster memebership in a new column ###NOTA tengo que corregir esto no puede ser full
df_seasonal_subsets_10reps$Obs<- gsub("_cluster_([0-9]*)", "", df_seasonal_subsets_10reps$Obs) #remove the information of cluster from the name of the contig
df_seasonal_subsets_10reps$Gral_ID<-paste(df_seasonal_subsets_10reps$ID,df_seasonal_subsets_10reps$cluster, sep="_") #Create a new column with a general ID
#Finally we need to create a centroid column linking the contig (a.k.a observation) to its correspending centroid in the selected subset
centroids_final<-tmp2_merged_klusters$Row.names

df_seasonal_subsets_10reps$centroid <- ifelse(df_seasonal_subsets_10reps$Obs %in% centroids_final, as.character(df_seasonal_subsets_10reps$Obs), "")

df_ctg_uniq<-df_seasonal_subsets_10reps[,c(4,5)] #just a link between centroid and uniqu cluster identifier
#remove columns without a centroid, this will leave only a dataframe with the centroid value and its unique id
df_ctg_uniq_NA<-mutate_all(df_ctg_uniq, list(~na_if(.,"")))
df_ctg_uniq_reduced<-na.omit(df_ctg_uniq_NA)

#######ARREGLAR ESTO Y YA HACER el run de prueba
#substitute centroid column with the corresponding centroid for the unique identifier #Quizas esto no esta funcionando
merged_df <- merge(df_seasonal_subsets_10reps, df_ctg_uniq_reduced, by = "Gral_ID", all.x = TRUE)

#remove the column ("centroid.x") that helped us to link the tables
final_repr_merged_df<-merged_df[ , -which(names(merged_df) %in% "centroid.x")]
#change the colname to match it in the next merge
colnames(tmp2_merged_klusters)[1]<-"centroid.y"

#Match the elemental contigs to the final centroid clustering determined by the different strategies 
Table_all_observations_clusters <- merge(final_repr_merged_df, tmp2_merged_klusters, by = "centroid.y", all.x = TRUE)

 FINAL_list[[i]] <- list(iteration = i, output_df = Table_all_observations_clusters)
}

##########Iteration ends

#Unfold the list -> necesitaria Table_all_observations_clusters mas una columna indicando la replica
bootstrap_df<-do.call(rbind, lapply(FINAL_list, function(x) {
  x$output_df$iteration <- x$iteration
  return(x$output_df)
}))

#Extract the clusters of each clustering method
data <- bootstrap_df[,c(4,7,9)] #klusters_tree_list (hierarchical)
colnames(data)<-c("ID", "Cluster", "Repetition")

data2 <- bootstrap_df[,c(4,6,9)] #(kmeans)
colnames(data2)<-c("ID", "Cluster", "Repetition") 

data3 <- bootstrap_df[,c(4,8,9)] #(TSclust)
colnames(data3)<-c("ID", "Cluster", "Repetition") 

unique_ids <- unique(data$ID)
unique_ids2 <- unique(data2$ID)
unique_ids3 <- unique(data3$ID)

#generate a coocurrence table in which the values are the total counts after running the machine learning step multiple times
#This value represents the number of times that 2 contigs cooccur in the same cluster after replicating the subsampling X times (in this case 50)
co_occurrence_matrix <- data %>%
  # Self-join to get all pairwise combinations of IDs
  full_join(data, by = "Cluster") %>%
  filter(ID.x != ID.y) %>% # Remove self-matches
  group_by(ID.x, ID.y) %>%
  # Count co-occurrences based on Cluster and Repetition
  summarise(co_occurrence_count = sum(Repetition.x == Repetition.y)) %>%
  # Convert the result to a matrix form
  pivot_wider(names_from = ID.y, values_from = co_occurrence_count, values_fill = 0) %>%
  column_to_rownames(var = "ID.x") %>%
  as.matrix()

#The data2 is giving an error: vector memory exhausted (limit reached?) I don't know why. 
co_occurrence_matrix2 <- data2 %>%
  # Self-join to get all pairwise combinations of IDs
  full_join(data2, by = "Cluster") %>%
  filter(ID.x != ID.y) %>% # Remove self-matches
  group_by(ID.x, ID.y) %>%
  # Count co-occurrences based on Cluster and Repetition
  summarise(co_occurrence_count = sum(Repetition.x == Repetition.y)) %>%
  # Convert the result to a matrix form
  pivot_wider(names_from = ID.y, values_from = co_occurrence_count, values_fill = 0) %>%
  column_to_rownames(var = "ID.x") %>%
  as.matrix()

co_occurrence_matrix3 <- data3 %>%
  # Self-join to get all pairwise combinations of IDs
  full_join(data3, by = "Cluster") %>%
  filter(ID.x != ID.y) %>% # Remove self-matches
  group_by(ID.x, ID.y) %>%
  # Count co-occurrences based on Cluster and Repetition
  summarise(co_occurrence_count = sum(Repetition.x == Repetition.y)) %>%
  # Convert the result to a matrix form
  pivot_wider(names_from = ID.y, values_from = co_occurrence_count, values_fill = 0) %>%
  column_to_rownames(var = "ID.x") %>%
  as.matrix()

#Express the value in terms of the number of iterations run 
#Get the last number of the vector: unique(data$Repetition)

total_iterations<-tail(sort(unique(data$Repetition)),1)
co_occurrence_matrix_rel<-co_occurrence_matrix/total_iterations
co_occurrence_matrix_rel2<-co_occurrence_matrix2/total_iterations
co_occurrence_matrix_rel3<-co_occurrence_matrix3/total_iterations

dendrogram = hclust(d = dist(co_occurrence_matrix_rel, method = 'euclidean'), method = 'ward.D')
dendrogram2 = hclust(d = dist(co_occurrence_matrix_rel2, method = 'euclidean'), method = 'ward.D')
dendrogram3 = hclust(d = dist(co_occurrence_matrix_rel3, method = 'euclidean'), method = 'ward.D')

#Cut the dendrogram using the highest number of optimal clusters estimated from all the iterations -> 57
highest_cl<-tail(sort(unique(data$Cluster)),1)

klusters_tree_list<-cutree(dendrogram, k=highest_cl)
klusters_tree_list2<-cutree(dendrogram2, k=highest_cl)
klusters_tree_list3<-cutree(dendrogram3, k=highest_cl)

dendro_best_number_contigcluster_df<-data.frame(klusters_tree_list)
dendro_best_number_contigcluster_df2<-data.frame(klusters_tree_list2)
dendro_best_number_contigcluster_df3<-data.frame(klusters_tree_list3)

dendro_best_number_contigcluster_df$obs<-row.names(dendro_best_number_contigcluster_df)
dendro_best_number_contigcluster_df2$obs<-row.names(dendro_best_number_contigcluster_df2)
dendro_best_number_contigcluster_df3$obs<-row.names(dendro_best_number_contigcluster_df3)

#Join by observation
dendro_best_number_contigcluster_df_merged<-merge(dendro_best_number_contigcluster_df, dendro_best_number_contigcluster_df2, by="obs")
dendro_best_number_contigcluster_df_merged<-merge(dendro_best_number_contigcluster_df_merged, dendro_best_number_contigcluster_df3, by="obs")

#Write the lists of the membership of the clusters for following steps
output_file <- paste0(output_path, "ChronoClustR_output_table.txt")
write.table(dendro_best_number_contigcluster_df_merged, file=output_file, quote=FALSE, sep = "\t")

#Check how many elements does each cluster have, this is just a sanity check step / uncomment if needed
#as.data.frame(table(dendro_best_number_contigcluster_df$klusters_tree_list))

#######Generate the BIG PLOT #########
######with panels=No of clusters######

#declare a list 
plot_list <- list()
plot_list2 <- list()
plot_list3 <- list()

#loop to generate all the plots of each cluster using the function function_plotsTS, this function extracts the subset of contigs for a specified cluster 
for (i in 1:highest_cl){
plot<-function_plotsTS(i,"klusters_tree_list",dendro_best_number_contigcluster_df) #HC
plot_list[[i]] <- plot}

for (i in 1:highest_cl){
plot2<-function_plotsTS(i,"klusters_tree_list2",dendro_best_number_contigcluster_df2) #kmeans
plot_list2[[i]] <- plot2}

for (i in 1:highest_cl){
plot3<-function_plotsTS(i,"klusters_tree_list3",dendro_best_number_contigcluster_df3) #TSclust
plot_list3[[i]] <- plot3}

#Print all, user can also print one by one by modifying the above for loop
output_plot_HC <- paste0(output_path, "ChronoClustR_treelist_HC.pdf")
plot_grob <- arrangeGrob(grobs=plot_list)
pdf(output_plot_HC)
grid.arrange(plot_grob)
dev.off()

output_plot_kmeans <- paste0(output_path, "ChronoClustR_treelist_kmeans.pdf")
plot_grob2 <- arrangeGrob(grobs=plot_list2)
pdf(output_plot_kmeans)
grid.arrange(plot_grob2)
dev.off()

output_plot_TSclust <- paste0(output_path, "ChronoClustR_treelist_TSclust.pdf")
plot_grob3 <- arrangeGrob(grobs=plot_list3)
pdf(output_plot_TSclust )
grid.arrange(plot_grob3)
dev.off()

#Uncomment to add the clustering information based on the clustering of the centroids

#Write the table comparing the clustering with the optimal number of clusters using the 3 different methods (centroids-representatives)
#write.table(tmp2_merged_klusters,"clusters_3methods.txt", quote = F, sep = "\t")

#Write the table of the elements and its corrsponding cluster -> this may be used in the future to generate a training dataset
#write.table(training_dataset_final,"pam_centroids_cluster.txt", quote = F, sep = "\t")

#Write a table with the final report number of datasets vs final number of distinct profile centroids for it. 
#write.table(df_report_benchmarck,"report_benchmarck.txt", quote = F, sep = "\t")

#Write a table with the final report with all the observations clustered by the three differetn methods, and the centroid used to classify them.
#write.table(Table_all_observations_clusters,"observations_clusters_3methods.txt", quote = F, sep = "\t")
