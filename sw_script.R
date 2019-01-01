library("jsonlite")
library("curl")
library("igraph")
library("ggplot2")
library("parallel")

corpusname <- "rus"
##plays_to_remove <- list("blok-balaganchik", "blok-korol-na-ploschadi", "blok-neznakomka", "gogol-teatralnyi-razezd")

## Downloading plays
list_of_names <- fromJSON(paste0("https://dracor.org/api/corpora/", corpusname))
sorted_ids <- list_of_names$dramas$id[sort.list(list_of_names$dramas$id)]
##sorted_ids <- sorted_ids[sorted_ids != plays_to_remove]  ## removing of plays which do not represent social interactions

downloader <- function(playname){
  if (file.exists(paste0("csv/", playname, ".csv"))) {
  
  } else {
  download.file(paste0("https://dracor.org/api/corpora/", corpusname, "/play/", playname, "/networkdata/csv"), paste0("csv/", playname, ".csv"))
  }
  read.csv(paste0("csv/", playname, ".csv"), stringsAsFactors = F)
}

plays <- lapply(sorted_ids, downloader)

#plays <- mclapply(sorted_ids, function(x) read.csv(paste0("https://dracor.org/api/corpora/", corpusname, "/play/", x, "/networkdata/csv"), stringsAsFactors = F))

## Remove 'Type' and 'Weight' variables
del_vars <- function(play){
  play$Type <- NULL
  play$Weight <- NULL
  return (play)
}
plays <- mclapply(plays, del_vars)

metadata <- read.csv(paste0("https://dracor.org/api/corpora/", corpusname, "/metadata.csv"), stringsAsFactors = F)
metadata <- metadata[order(metadata$name),]
##metadata <- metadata[metadata$name != plays_to_remove,] ## removing of plays which do not represent social interactions

## Creating graphs of plays
graphs_of_plays <- lapply(plays, function(x) graph_from_data_frame(x, directed = F))

## Calculations of clustering coefficient and average path length
CC <- sapply(graphs_of_plays, transitivity)
APL <- sapply(graphs_of_plays, function(x) mean_distance(x, directed=F))

df <- subset(metadata, select = c(name, year, numOfSpeakers) )
df$CC = CC
df$APL = APL

## Function generates random graphs and calculates metrics (CC and APL) for them
set.seed(42)
randomize_graph <- function(graph){
  random_graphs=list(1000)
  for (i in 1:1000){
    random_graphs[[i]] <- sample_gnm(gorder(graph), gsize(graph), directed = F, loops = F)
  }
  ## random_graphs <- lapply(random_graphs, function(x) x <- sample_gnm(gorder(graph), gsize(graph), directed = FALSE, loops = FALSE))
  kilo_CC <- sapply(random_graphs, transitivity)
  kilo_APL <- sapply(random_graphs, function(x) mean_distance(x, directed=F))
  results <- list()
  results$CC_rand <- mean(kilo_CC)
  results$APL_rand <- mean(kilo_APL)
  return(results)
}

metrics_for_rand_graphs <- lapply(graphs_of_plays, randomize_graph)
matrix_metrics <- t(sapply(metrics_for_rand_graphs, unlist))

df$CC_rand=matrix_metrics[,1]
df$APL_rand=matrix_metrics[,2]

df <- transform(df, CC_dev = CC / CC_rand )
df <- transform(df, APL_dev = APL / APL_rand )

is.na(df) <- do.call(cbind,lapply(df, is.infinite))

## Calculating border values
CC_border <- mean(df$CC_dev, na.rm = TRUE)+2*sd(df$CC_dev, na.rm = TRUE)
APL_border_min <- mean(df$APL_dev, na.rm = TRUE)-2*sd(df$APL_dev, na.rm = TRUE)
APL_border_max <- mean(df$APL_dev, na.rm = TRUE)+2*sd(df$APL_dev, na.rm = TRUE)

## Applying criteria
df$crit_1 <- ifelse(df$CC_dev>=CC_border, TRUE, FALSE)
df$crit_2 <- ifelse((df$crit_1 == TRUE) & (df$APL_dev >= APL_border_min) & (df$APL_dev<=APL_border_max), TRUE, FALSE)

## Comparison of R-squared for power law ##

## preprocessing

degree_dist <- lapply(graphs_of_plays, function(x) degree(x, v = V(x), loops = FALSE, normalized = FALSE))
degree_dist_v <- lapply(degree_dist, as.vector)
number_of_nodes <- lapply(degree_dist_v, table)
distribution <- lapply(number_of_nodes, as.data.frame)

num_type <- function(x){
  x$Var1 <- as.numeric(as.character(x$Var1))
  x$Freq <- as.numeric(x$Freq)
  names(x) <- c("Node_degree", "Num_of_nodes")
  x
}

distribution <- lapply(distribution, num_type)

## regressions
fit <- lapply(distribution, function(x) lm(log(x$Num_of_nodes) ~ log(x$Node_degree)))
df$Rsqrd <- sapply (fit, function(x) summary(x)$r.squared)

## The dataframe with small-world networks
small_worlds <- na.omit(df[df$crit_2 == TRUE,])

## The list of small-world graphs
## small_worlds_g <- lapply(as.numeric(rownames(small_worlds)), function(x) graphs_of_plays[[x]])


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
  coord_cartesian(xlim = c(0, max(df$numOfSpeakers)+12))+
  geom_text(aes(label=ifelse(numOfSpeakers>75,as.character(name),'')),hjust=0.5,vjust=1.5, color="black", size=4)
plot2


power_law_plot <- function(name_of_play){
  id <- which(name_of_play == df$name)[[1]]
  loglog <- lm(log(distribution[[id]]$Num_of_nodes) ~ log(distribution[[id]]$Node_degree))
  loglog_df <- data.frame(x = distribution[[id]]$Node_degree, y = exp(fitted(loglog)))
  
  ggplot(data = distribution[[id]], aes(x = Node_degree, y = Num_of_nodes))+ 
    geom_point(size = 3) + 
    geom_line(data = loglog_df, aes(x = x, y = y), linetype = 2, color="blue", size = 1)+
    labs(x="Node degree", y = "Number of nodes")
}

power_law_plot("goethe-goetz-von-berlichingen-mit-der-eisernen-hand")
