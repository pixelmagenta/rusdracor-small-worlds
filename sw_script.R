library("jsonlite")
library("curl")
library("igraph")
library("ggplot2")
library("parallel")

corpora <- "rus"
## plays_to_remove <- list("blok-balaganchik", "blok-korol-na-ploschadi", "blok-neznakomka", "gogol-teatralnyi-razezd")

## Downloading plays
list_of_names <- fromJSON(paste0("https://dracor.org/api/corpora/", corpora))
sorted_ids <- list_of_names$dramas$id[sort.list(list_of_names$dramas$id)]
## sorted_ids <- sorted_ids[sorted_ids != plays_to_remove]  ## removing of plays which do not represent social interactions
plays <- mclapply(sorted_ids, function(x) read.csv(paste0("https://dracor.org/api/corpora/", corpora, "/play/", x, "/networkdata/csv"), stringsAsFactors = F))

## Remove 'Type' and 'Weight' variables
del_vars <- function(play){
  play$Type <- NULL
  play$Weight <- NULL
  return (play)
}
plays <- mclapply(plays, del_vars)

metadata <- read.csv(paste0("https://dracor.org/api/corpora/", corpora, "/metadata.csv"), stringsAsFactors = F)
metadata <- metadata[order(metadata$name),]
## metadata <- metadata[metadata$name != plays_to_remove,] ## removing of plays which do not represent social interactions

## Creating graphs of plays
graphs_of_plays <- mclapply(plays, function(x) graph_from_data_frame(x, directed = F))

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


## PLOTS ##

theme_set(theme_gray(base_size = 15)) ## size of axis font
plot1 <- ggplot(na.omit(df), aes(x = year, y = CC_dev, 
                                 fill = crit_1, shape = crit_2, color = crit_1))+
  geom_point(size = 6)+
  scale_shape_manual(values=c(21, 22), name = "Criterion 2")+
  scale_color_manual(values=c("slateblue3", "indianred3"), name = "Criterion 1")+
  scale_fill_manual(values=c("slateblue1", "indianred1"), name = "Criterion 1")+
  labs(x="Year of creation", y = "Clustering coefficient deviation")
plot1


theme_set(theme_gray(base_size = 15)) ## size of axis font
plot2 <- ggplot(na.omit(df), aes(x = numOfSpeakers, y = CC_dev, 
                                   fill = crit_1, shape = crit_2, color = crit_1))+
  geom_point(size = 6)+
  scale_shape_manual(values=c(21, 22), name = "Criterion 2")+
  scale_color_manual(values=c("slateblue3", "indianred3"), name = "Criterion 1")+
  scale_fill_manual(values=c("slateblue1", "indianred1"), name = "Criterion 1")+
  labs(x="Number of characters", y = "Clustering coefficient deviation")+
  coord_cartesian(xlim = c(0, 105))+
  geom_text(aes(label=ifelse(numOfSpeakers>47,as.character(name),'')),hjust=0.5,vjust=1.5, color="black", size=4)
plot2


## Comparison of R-squared for power law ##

## preprocessing

degree_dist <- lapply(graphs_of_plays, function(x) degree(x, v = V(x), loops = FALSE, normalized = FALSE))
degree_dist_v <- lapply(degree_dist, as.vector)
number_of_nodes <- lapply(degree_dist_v, table)
distribution <- lapply(number_of_nodes, as.data.frame)

## regressions
fit <- lapply(distribution, function(x) lm(log(as.numeric(x$Freq)) ~ log(as.numeric(x$Var1))))
df$Rsqrt <- sapply (fit, function(x) summary(x)$r.squared)

## additional plot for play with PLAY_NUMBER from df
##ggplot(data = distribution[[PLAY_NUMBER]], aes(x = as.numeric(Var1), y = as.numeric(Freq)))+ 
##  geom_point(size = 3) + 
##  geom_smooth(method = "lm", formula=log(y)~log(x))+
##  labs(x="Node degree", y = "Number of nodes")

