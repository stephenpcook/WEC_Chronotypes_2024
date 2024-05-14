#This script was used to generate the figure 03: Seasonal and non-seasonal chronotypes derived from an unsupervised machine learning clustering method

#README: This script represents steps 12 and 13 from the ChronoClustR analysis represented in Figure 1. The code below generates a co-occurence matrix based 
#on the replicates generated from steps 6-10 (11 n bootstrap replicates).

###DEFINE FUNCTIONS###

#NOW IN Z-SCORE as it is in the overall table
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

#FUNCTION 4 Euc_dist:
euclidean_distance_func <- function(x, y) {
  sqrt(sum((x - y) ^ 2))
}

# Step 4: Calculate pairwise distances using 'dist' and the custom function
# We use 'as.matrix' to convert the data frame to a matrix before computing the distances.
function_Within_MeanEuclDist <- function(plots,column_name,DataFrame){
cluster<- DataFrame[DataFrame[[column_name]]==plots,]$obs
#Here subset the zscores 
cluster_df_zscore<-NONSignSeasonal_rpkms.normz_melt[NONSignSeasonal_rpkms.normz_melt$contig %in% cluster, ]
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

#SEASONALITY TEST#
countwec <- read.table("rpkms_populations.txt", header=T, row.names=1, check.names=F)
countwect<-t(countwec)

VIR_cl_env<-read.table("virome_sr_environm_data_EQ_time.txt", header=T, row.names=1, check.names=F, sep ="\t")
VIR_cl_env$sample<-rownames(VIR_cl_env)
VIR_cl_envdate<-VIR_cl_env[ , c("sample","Date_wavelet")]
VIR_cl_envdate$sample<-paste0("X",1:nrow(VIR_cl_envdate),".",VIR_cl_envdate$sample)

monthly_RPKM<-merge(VIR_cl_envdate,countwect,by=0, all = TRUE) ######I need to interpolate this
monthly_RPKM$Date_wavelet<-as.Date(monthly_RPKM$Date_wavelet,"%d/%m/%Y")
monthly_RPKM_sorted<-(monthly_RPKM[order(monthly_RPKM$Date_wavelet),]) #dim 32 26854

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

#Subset these list from monthly_RPKM_sort_trim_int_cl_melt into two different data frames 
Seasonal_significant_DF_melt <- monthly_RPKM_sort_trim_int_cl_melt[monthly_RPKM_sort_trim_int_cl_melt$variable %in% highly_seasonaldf, ]
not_significant_DF_melt <- monthly_RPKM_sort_trim_int_cl_melt[monthly_RPKM_sort_trim_int_cl_melt$variable %in% non_seasonaldf, ]

#reshape the data frames
Seasonal_significant_DF<-dcast(data = Seasonal_significant_DF_melt,formula = Date_wavelet~variable,value.var = "value") #dim Seasonal_significant_DF 32 22697
not_significant_DF<-dcast(data = not_significant_DF_melt,formula = Date_wavelet~variable,value.var = "value") 

#Remove the first column (date wavelet)
nums <- Seasonal_significant_DF[,2:length(Seasonal_significant_DF)] #Remove the wavelet column (1)
nums_nots<-not_significant_DF[,2:length(not_significant_DF)] #Remove the wavelet column (1)

#transpose both - These are the objects that can also be inputted if user has already checked the seasonality of its data
SignSeasonal_rpkms_subset<-t(nums) #contigs are rows #write.table(tnums, "/PATH/ctgs_interpolatedRPKMS.seasonal22k.tab", quote=FALSE, sep = "\t")
notsign_rpkms_subset<-t(nums_nots) #write.table(tnums_nots, "/PATH/ctgs_interpolatedRPKMS.notseasonal4k.tab", quote=FALSE, sep = "\t")

#Rename the column headers -> automatize this
colnames(SignSeasonal_rpkms_subset) <-VIR_cl_envdate$sample 
colnames(notsign_rpkms_subset)<- VIR_cl_envdate$sample 

#z-score normalisation
SignSeasonal_rpkms.normz <- data.frame(zscore(SignSeasonal_rpkms_subset))
SignSeasonal_rpkms.normz$contig<-row.names(SignSeasonal_rpkms.normz)

NONSignSeasonal_rpkms.normz <- data.frame(zscore(notsign_rpkms_subset))
NONSignSeasonal_rpkms.normz$contig<-row.names(NONSignSeasonal_rpkms.normz)

SignSeasonal_rpkms.normz_melt<-melt(SignSeasonal_rpkms.normz)
NONSignSeasonal_rpkms.normz_melt<-melt(NONSignSeasonal_rpkms.normz)

#write.table(NONSignSeasonal_rpkms.normz_melt, "NONSignSeasonal_rpkms.normz_melt.txt", row.names = FALSE, col.names = TRUE, sep = "\t")

#Now we have the necessary data to extract from the arrange of clusters

#END of SEASONALITY TEST#

######MATRIX co-occurence######

#Read the file generated by merging the clusters of the 50 replicates (done individually in previous steps / looped bash script)
data <- read.table("repetition_datasets3k_window100_nonseasonal_100reps/klusters_tree_list_rep100.txt", header = FALSE)
colnames(data)<-c("ID", "Cluster", "Repetition")

data2<-read.table("repetition_datasets3k_window100_nonseasonal_100reps/kmeans_list_rep100.txt", header = FALSE)
colnames(data2)<-c("ID", "Cluster", "Repetition") 

data3<-read.table("repetition_datasets3k_window100_nonseasonal_100reps/ktsclust_list_rep100.txt", header = FALSE)
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

#The data2 is giving an error: vector memory exhausted (limit reached?) I don't know why. Same with non-seasonal 
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

#Express the value in terms of the number of iterations run (50)
co_occurrence_matrix_rel<-co_occurrence_matrix/100
co_occurrence_matrix_rel2<-co_occurrence_matrix2/1000
co_occurrence_matrix_rel3<-co_occurrence_matrix3/100

#sanity check - plot a heatmap
#heatmap_centroids<-co_occurrence_matrix_rel %>%
#  gplots::heatmap.2 (
    # dendrogram control
#    hclustfun = hclust,
#    dendrogram = "both",
#    tracecol=NA,
#    labRow = FALSE, labCol=FALSE
#  )

#If you see in the plot a diagonal of similiarities (colored squares through the diagonal), it means that we have meaningful co-occurences
#Now, basically we are going to generate the same dendrogram as the ones on the sides of the heatmap
dendrogram = hclust(d = dist(co_occurrence_matrix_rel, method = 'euclidean'), method = 'ward.D')
dendrogram2 = hclust(d = dist(co_occurrence_matrix_rel2, method = 'euclidean'), method = 'ward.D')
dendrogram3 = hclust(d = dist(co_occurrence_matrix_rel3, method = 'euclidean'), method = 'ward.D')

#Cut the dendrogram using the highest number of optimal clusters estimated from all the iterations
highest_cl<-46 #this is the same for all the clustering methods #Here we are just analyzing non-seasonal  (left panels of figure 03)

klusters_tree_list<-cutree(dendrogram, k=highest_cl)
klusters_tree_list2<-cutree(dendrogram2, k=highest_cl)
klusters_tree_list3<-cutree(dendrogram3, k=highest_cl)

dendro_best_number_contigcluster_df<-data.frame(klusters_tree_list)
dendro_best_number_contigcluster_df2<-data.frame(klusters_tree_list2)
dendro_best_number_contigcluster_df3<-data.frame(klusters_tree_list3)

#We turn the row.names into an "observations" column
dendro_best_number_contigcluster_df$obs<-row.names(dendro_best_number_contigcluster_df)
dendro_best_number_contigcluster_df2$obs<-row.names(dendro_best_number_contigcluster_df2)
dendro_best_number_contigcluster_df3$obs<-row.names(dendro_best_number_contigcluster_df3)

#Join by observation

dendro_best_number_contigcluster_df_merged<-merge(dendro_best_number_contigcluster_df, dendro_best_number_contigcluster_df3, by="obs")

#Write the lists of the membership of the clusters for following steps
#write.table(dendro_best_number_contigcluster_df_merged, "NONseas_HC_TsClust_contig_memb.txt", quote=FALSE, sep = "\t")

#Check how many elements does each cluster have, this is just a sanity check step / uncomment it if needed
#as.data.frame(table(dendro_best_number_contigcluster_df$klusters_tree_list))

#######Generate the multipanel plot showing the Z-score for each chronotype#########
######in the original we did not plot the number of chronotype, here this is included to track which panel represent what chronotype######

#FUNCTION 3 plots: This function needs a number (that will reflect which cluster)
#I dont want to do plot grid, just to have more versatility to explore specific clusters for now. 

#plots == a number
#column_name == this name changes depending on the method used to generate the cluster in the first script for example klusters_tree_list, kmeans, unsup_cluster
#If running for kmeans, remember to remove the string "cluster_" from the values, so the final input has only numbers. This function only runs on numbers

#Set the limits of the y coordinates: Based on the max and min value of the z score in SignSeasonal_rpkms.normz_melt$value
maxzscore<-ceiling(max(SignSeasonal_rpkms.normz_melt$value)) #set max 6
minzscore<-floor(min(SignSeasonal_rpkms.normz_melt$value)) #set min -4

function_plotsTS <- function(plots,column_name,DataFrame){
cluster<- DataFrame[DataFrame[[column_name]]==plots,]$obs
#Here subset the zscores 
cluster_df_zscore<-NONSignSeasonal_rpkms.normz_melt[NONSignSeasonal_rpkms.normz_melt$contig %in% cluster, ]
#PLOT
ggplot(cluster_df_zscore,aes(x = variable, y = value, group=contig)) +
geom_line(size=.2,aes( x=variable,y=value)) +
scale_x_discrete(expand = c(0.01,0))+
scale_y_continuous(limits = c(minzscore, maxzscore))+
geom_vline(xintercept = c(3,15,27), color = "black")+theme_bw()+
guides(color=FALSE)+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.text.x = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),axis.ticks.x=element_blank())
  }

#declare a list 
plot_list <- list()
#plot_list2 <- list()
plot_list3 <- list()

#loop to generate all the plots of each cluster using the function function_plotsTS, this function extracts the subset of contigs for a specified cluster 
for (i in 1:highest_cl){
plot<-function_plotsTS(i,"klusters_tree_list",dendro_best_number_contigcluster_df)
plot_list[[i]] <- plot}

#for (i in 1:highest_cl){
#plot2<-function_plotsTS(i,"klusters_tree_list2",dendro_best_number_contigcluster_df2)
#plot_list2[[i]] <- plot2}

for (i in 1:highest_cl){
plot3<-function_plotsTS(i,"klusters_tree_list3",dendro_best_number_contigcluster_df3)
plot_list3[[i]] <- plot3}

#Print all, user can also print one by one by modifying the above for loop
plot_grob <- arrangeGrob(grobs=plot_list)
#pdf("klusters_tree_list_rep100.pdf")
#svg("klusters_tree_list_rep100.svg")
grid.arrange(plot_grob)
#dev.off()

plot_grob2 <- arrangeGrob(grobs=plot_list2)
#pdf("kmeans_list_rep100.pdf")
grid.arrange(plot_grob2)
#dev.off()

plot_grob3 <- arrangeGrob(grobs=plot_list3)
#pdf("ktsclust_list_rep100.pdf")
#svg("ktsclust_list_rep100.svg")
grid.arrange(plot_grob3)
#dev.off()


######################################
#Estimate the eucladian distances within cluster and between clusters (mean) for the denoised HierClust tree 50 reps -> 107 points (distribution) stats between them

library(TSdist)

#declare a list
estimated_mean_euc_dist <- list()
estimated_mean_euc_dist2 <- list()
estimated_mean_euc_dist3 <- list()

#loop to generate all the mean euc distances of all vs all within a cluster. 
for (i in 1:highest_cl){
eucdist<-function_Within_MeanEuclDist(i,"klusters_tree_list",dendro_best_number_contigcluster_df)
estimated_mean_euc_dist[[i]] <- eucdist}

estimated_mean_euc_dist_vector<-unlist(estimated_mean_euc_dist, use.names=FALSE) #this vector has the 107 mean distances of the 107 clusters

#loop to generate all the mean euc distances of all vs all within a cluster. 
for (i in 1:highest_cl){
eucdist2<-function_Within_MeanEuclDist(i,"klusters_tree_list2",dendro_best_number_contigcluster_df2)
estimated_mean_euc_dist3[[i]] <- eucdist3}

estimated_mean_euc_dist_vector2<-unlist(estimated_mean_euc_dist2, use.names=FALSE)

#loop to generate all the mean euc distances of all vs all within a cluster. 
for (i in 1:highest_cl){
eucdist3<-function_Within_MeanEuclDist(i,"klusters_tree_list3",dendro_best_number_contigcluster_df3)
estimated_mean_euc_dist3[[i]] <- eucdist3}

estimated_mean_euc_dist_vector3<-unlist(estimated_mean_euc_dist3, use.names=FALSE)

#write out both vectors to 
#write(estimated_mean_euc_dist_vector, "klusters_tree_list_rep100.dists")
#write(estimated_mean_euc_dist_vector2, "kmeans_list_rep100.dists")
#write(estimated_mean_euc_dist_vector3, "ktsclust_list_rep100.dists")

#ALLL of the above was done to get the avg distances of the denoised. Do a stat test with Ho= all these distances derived from the same population. Meaning all are the same, there are no significant differences
estimated_mean_euc_dist_vector <- scan("klusters_tree_list_rep100.dists", what=numeric(), sep=" ")
estimated_mean_euc_dist_vector3<-scan("ktsclust_list_rep100.dists", what=numeric(), sep=" ")

kluster_means_df <- data.frame(HC = estimated_mean_euc_dist_vector , TsClust = estimated_mean_euc_dist_vector3)
kluster_means_df_melted<-melt(kluster_means_df)

#t.test(value ~ variable, data = kluster_means_df_melted)
#t = 0.94414, df = 89.971, p-value = 0.3476

#mean(kluster_means_df$HC) 2.543202
#mean(kluster_means_df$TsClust) 2.294201

boxplots_means<-ggplot(kluster_means_df_melted, aes(x = variable, y = value))+geom_boxplot(size = 0.2, color = "black",position = position_dodge2(0.75))+
geom_point() + ylab("within cluster mean distance")+ xlab("clustering method") + ylim(0,9)+
theme(panel.background = element_blank(),legend.position = "none",strip.background = element_blank(),axis.text.y=element_text(colour = 'black',size=14), axis.text.x=element_text(colour = 'black',size=14),axis.line = element_line(colour = "black"),axis.title.x=element_text(colour = 'black',size=14),axis.title.y=element_text(colour = 'black',size=12) )

#svg("boxplots_means_dists.svg", width=7, height=7)
boxplots_means
#dev.off()

#We will select TSclust clusters

################################################
#Function to retrieve de mediod of each cluster#
################################################

#We want to get this to use it as a representative - Use dendro_best_number_contigcluster_df3 first

#We are using two medoids of the cluster as representatives of if the cluster is larger than 2 (as a control, if the cluster is hiding other temporal patterns, if not, both will used for the corraltions. 

centroid_list = vector("list",length(highest_cl)) #Define a vector 

for (i in 1:highest_cl) {
  kluster_subset <- dendro_best_number_contigcluster_df3[dendro_best_number_contigcluster_df3$klusters_tree_list3 == i, ][,2] #subset the first ddataframe with the obs (names) and kluster membership
  if (length(kluster_subset)>2){
  cluster_df_zscore<-NONSignSeasonal_rpkms.normz_melt[NONSignSeasonal_rpkms.normz_melt$contig %in% kluster_subset, ] #subset the zscores based on the names from the above step
  unmelted_zscore<-t(data.frame(cluster_df_zscore  %>% pivot_wider(names_from = contig, values_from = value)))
  #remove the first non numeric line
  colnames(unmelted_zscore) <- unmelted_zscore[1,]
  unmelted_zscore <- data.frame(unmelted_zscore[-1, ])
  #turn the numbers into numeric not factor
  unmelted_zscore <-unmelted_zscore %>% mutate_if(is.character, as.numeric)
  #generate a 2 K clustering of the already defined cluster. Get Both as two representatives. #The underlying idea, is that if the cluster is tight either should work, if not, we use two to cover more of the dynamics. 
  bestdefined_clustering.pam <- tsclust(unmelted_zscore , type="partitional", k=2, distance="Euclidean", centroid="pam") #generate a mediod of this specific dataset
  bestdefclust.centroids<-attr(bestdefined_clustering.pam@centroids,"series_id")
  centroid_list[[i]]<-c(rownames(unmelted_zscore[bestdefclust.centroids,]))}
  else{
  centroid_list[[i]]<- kluster_subset}
  }

#Add the "kluster name" into the centroid_list
NAME <- paste0(1:length(centroid_list),"ns")
names(centroid_list) <- NAME

#Extend the list saving its clustering representation. 
nonseasonal_repr_centroids<-data.frame(kluster = rep(names(centroid_list), sapply(centroid_list, length)), Obs = unlist(centroid_list))

#write.table(nonseasonal_repr_centroids, "tsclust_reprs_NS.txt", quote=FALSE, sep = "\t")
