# ChronoClustR.R

Beta version - May 2024

## Instructions

(1) Create a specific conda environment

ChronoClustR was benchmarked in R version 4.3.0 (2023-04-21)

(2) Install R and the required R packages

(3) User can decide to run a seasonality test using the script `SeasonalityTest.R`

- User should specify the time-series data frame, environmental data (formatted as in the example `virome_sr_environm_data_EQ_time.txt`)
- Time-series should be equidistant with no missing values, therefore if a time-series have missing values, these should have been interpolated before running the seasonality test.

Command:

```sh
Rscript SeasonalityTest.R [Time-series_table.txt] [metadata.txt] [output.dir/] [start year] [start month] [end year] [end month] [frequency]
```

Example:

```sh
Rscript SeasonalityTest.R rpkms_populations_600.txt virome_sr_environm_data_EQ_time.txt Test_ChronoClustR/ 2018 11 2021 6 32
```

This script creates two files in the `Test_ChronoClustR/` dir:

- `Non_seasonal_timeseries.txt`
- `Seasonal_timeseries.txt`

(4) If seasonality was tested, each of both tables (seasonal and non-seasonal) can be run separately in `ChronoClustR.R`.

Command:

```sh
Rscript ChronoClusteR.R [Time-series_table.txt] [metadata.txt] [Bootstrap replicates] [Subsample size] [output.dir/]
```

Example:

```sh
Rscript ChronoClusteR.R Test_ChronoClustR/Non_seasonal_timeseries.txt virome_sr_environm_data_EQ_time.txt 20 21 Test_ChronoClustR/Non_seasonal_Clust/
Rscript ChronoClusteR.R Test_ChronoClustR/Seasonal_timeseries.txt virome_sr_environm_data_EQ_time.txt 20 55 Test_ChronoClustR/Seasonal_Clust/
```

> [!Note]
> User has the option to run clustering without a prior seasonality test just by running the full dataset in ChronoClustR and skipping step 3

(5) 5 files will be generated (the output files of the run are provided in the `Test_ChronoClustR` directory):
`ChronoClustR_output_table.txt`: This file contains the time-series identifier and the cluster membership estimated by the three methods (`HC`, `Kmeans`, `TSclust`)

`Column klusters_tree_list = Hierarchical Clustering`
`Column klusters_tree_list2 = Kmeans`
`Column klusters_tree_list3 = TSclust`

Only `TSclust` has been tested on our dataset.

Example:

```none
obs	klusters_tree_list	klusters_tree_list2	klusters_tree_list3
1	X2018.Dec.viral.fraction__viral.spades.short__as.yet.unknown__18062bp__001309..full	1	1	1
2	X2018.Dec.viral.fraction__viral.spades.short__as.yet.unknown__33569bp__000363..full	2	1	2
3	X2018.Dec.viral.fraction__viral.spades.short__as.yet.unknown__84964bp__000032..full	2	1	2
```

`Rplots.pdf`: contains plots of the clusters-distance within cluster for each bootstrap iteration, heatmaps of the centroids dendrogram used to generate the final clustering, and the maximum number of K used in the collapsed co-occurrence frequency can be derived from analyzing each of the K for each iteration in this report

Plots of individual time-series clustered in different panels depending their membership estimated by the three different methods described above

- `ChronoClustR_treelist_HC.pdf`
- `ChronoClustR_treelist_kmeans.pdf`
- `ChronoClustR_treelist_TSclust.pdf`

User should modify the `function_plotsTS` [Function 4] in order to generate the plots according to the time values of the input (this will be fixed in a following version)
