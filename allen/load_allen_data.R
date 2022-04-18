library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(rdist)

df <- read.csv("../data/mouse/soma_subgraph_synapses_spines_v185.csv",  colClasses = "character")
loc <- read.csv("./data/mouse/soma_valence_v185.csv", colClasses = "character")
df$post_root_id %>% unique %>% length
df$pre_root_id %>% unique %>% length

W_counts = matrix(0, df$post_root_id %>% unique %>% length, df$post_root_id %>% unique %>% length)
W_weights = matrix(0, df$post_root_id %>% unique %>% length, df$post_root_id %>% unique %>% length)
rownames(W_counts) <- colnames(W_counts) <- (df$post_root_id %>% unique)
rownames(W_weights) <- colnames(W_weights) <- (df$post_root_id %>% unique)

for (i in 1:nrow(df)) {
  pre = df[[i, "pre_root_id"]]
  post = df[[i, "post_root_id"]]
  w = df[[i, "spine_vol_um3"]] %>% as.double
  W_counts[pre,post] = W_counts[pre,post]+1
  W_weights[pre,post] = W_weights[pre,post] + w
}



W_bin <- W_counts>0


loc <- loc %>% filter(pt_root_id %in% df$post_root_id)
loc <- loc[order(match((loc$pt_root_id), colnames(W_counts))),]
loc <- loc %>% mutate(clean_pos = str_sub(pt_position, 2, -2)) %>%  
  separate(clean_pos, c('x','y','z'), sep="\\s+") 

loc <- loc %>% mutate(x=as.double(x), y=as.double(y), z=as.double(z))

joined_on_pre <- left_join(df, loc %>% mutate(pre_root_id=pt_root_id), by="pre_root_id")
joined_on_post <- left_join(df, loc %>% mutate(post_root_id=pt_root_id), by="post_root_id")

# hist(loc$y, breaks=100)
# hist(joined_on_pre$y, breaks=100, add=T)

dm <- pdist(loc %>% select(x,y,z))
dm %>% dim
fit <- hclust(as.dist(dm))

dm <- dm[fit$order, fit$order]
W_bin <- W_bin[fit$order, fit$order]
W_counts <- W_counts[fit$order, fit$order]
W_weights <- W_weights[fit$order, fit$order]
loc <- loc[fit$order, ]

loc <- loc %>% mutate(indegs = colSums(W_bin),
                      outdegs = rowSums(W_bin))

# build g
g <- network(W_bin, directed = TRUE, vertex.attr = loc)

ggpairs(loc %>% select(x,y,z, indegs, outdegs))

