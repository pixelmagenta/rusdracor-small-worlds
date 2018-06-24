library("jsonlite")
library("curl")
library("igraph")

plays_to_remove <- list("blok-balaganchik", "blok-korol-na-ploschadi", "blok-neznakomka", "gogol-teatralnyi-razezd")

## Downloading plays
list_of_names <- fromJSON("https://dracor.org/api/corpus/rus")
sorted_ids <- list_of_names$dramas$id[sort.list(list_of_names$dramas$id)]
sorted_ids <- sorted_ids[sorted_ids != plays_to_remove]  ## removing of plays which do not represent social interactions
plays <- lapply(sorted_ids, function(x) read.csv(paste0("https://dracor.org/api/corpus/rus/play/", x, "/networkdata/csv"), stringsAsFactors = F))

## Remove 'Type' and 'Weight' variables
del_vars <- function(play){
  play$Type <- NULL
  play$Weight <- NULL
  return (play)
}
plays <- lapply(plays, del_vars)

metadata <- read.csv("https://dracor.org/api/corpus/rus/metadata.csv", stringsAsFactors = F)

metadata <- metadata[metadata$name != plays_to_remove,] ## removing of plays which do not represent social interactions

## Creating graphs of plays
graphs_of_plays <- lapply(plays, function(x) graph_from_data_frame(x, directed = F))

CC <- sapply(graphs_of_plays, transitivity)
APL <- sapply(graphs_of_plays, function(x) mean_distance(x, directed=F))

df <- subset(metadata, select = c(name, year, numOfSpeakers) )
df$Clustering_coefficient = CC
df$Average_path_length = APL

## Function for creating random graphs and calculating metrics for them
randomize_graph <- function(graph){
  random_graphs=list(1000)
  random_graphs <- lapply(random_graphs, function(x) x <- sample_gnm(gorder(graph), gsize(graph), directed = FALSE, loops = FALSE))
  mile_CC <- sapply(random_graphs, transitivity)
  mile_APL <- sapply(random_graphs, function(x) mean_distance(x, directed=F))
  results <- list()
  results$CC_rand <- mean(mile_CC)
  results$APL_rand <- mean(mile_APL)
  return(results)
}

metrics_for_rand_graphs <- lapply(graphs_of_plays, randomize_graph)
matrix_metrics <- t(sapply(metrics_for_rand_graphs, unlist))

df$Clustering_coefficient_rand=matrix_metrics[,1]
df$Average_path_length_rand=matrix_metrics[,2]

df <- transform(df, CC_dev = Clustering_coefficient / Clustering_coefficient_rand )
df <- transform(df, APL_dev = Average_path_length / Average_path_length_rand )

is.na(df) <- do.call(cbind,lapply(df, is.infinite))

## Calculating border values
CC_border <- mean(df$CC_dev, na.rm = TRUE)+2*sd(df$CC_dev, na.rm = TRUE)
APL_border_min <- mean(df$APL_dev, na.rm = TRUE)-2*sd(df$APL_dev, na.rm = TRUE)
APL_border_max <- mean(df$APL_dev, na.rm = TRUE)+2*sd(df$APL_dev, na.rm = TRUE)

## Applying criteria
df$crit_1 <- ifelse(df$CC_dev>=CC_border, TRUE, FALSE)
df$crit_2 <- ifelse((df$crit_1 == TRUE) & (df$APL_dev >= APL_border_min) & (df$APL_dev<=APL_border_max), TRUE, FALSE)

## The dataframe with small-world networks
small_worlds <- na.omit(df[df$crit_2 == TRUE,])