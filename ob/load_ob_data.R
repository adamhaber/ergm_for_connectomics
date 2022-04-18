# get data
library(R.matlab)
library(rmatio)
library(dplyr)
library(ggplot2)
library(reshape2)

# OB data used with permission from Andrian Wanner and Friedrich Rainer
newdata <- readMat(...)
somas <- readMat(...)
somas <- as.data.frame(t(rbind(somas$NeuronIDs,somas$SomaCenter)))
colnames(somas) <- c("ID", "x", "y", "z")

Q <- newdata$Q
names(Q) <- names(Q[,,1])

cellinfo <- rbind(Q$CellinfoMC, Q$CellinfoINall)[,1:7]
colnames(cellinfo) <- c("ID", "Original index", "cell type", "microcluster", "GlomID", "Length", "Row_W")
cellinfo <- merge(cellinfo, somas, by="ID")
cellinfo <- as_tibble(cellinfo)

# align neuron IDs with rows from the connectivity matrix
in_rows <- c()
in_ids <- c()
for (i in 1:length(Q$NeuronIDs1003)) {
  if (length(which(Q$IDs[Q$INinds]==Q$NeuronIDs1003[i]))>0) {
    in_rows <- c(in_rows,i)
    in_ids <- c(in_ids, Q$NeuronIDs1003[i])
  }
}

mc_rows <- c()
mc_ids <- c()
for (i in 1:length(Q$NeuronIDs1003)) {
  if (length(which(Q$IDs[Q$MCinds]==Q$NeuronIDs1003[i]))>0) {
    if (i %in% Q$CellinfoMC[,7]) {
      mc_rows <- c(mc_rows,i)  
      mc_ids <- c(mc_ids, Q$NeuronIDs1003[i])
    }
  }
}

mc_rows_ord <- mc_rows[order(match(mc_rows, Q$CellinfoMC[,7]))]
mc_ids_ord <- mc_ids[order(match(mc_rows, Q$CellinfoMC[,7]))]

all_rows <- c(mc_rows_ord, in_rows)
all_ids <- c(mc_ids_ord, in_ids)

cellinfo <- cellinfo %>% filter(ID %in% all_ids) %>% filter(Length>0) %>% arrange(match(ID, all_ids))
idx <- cellinfo$Row_W
cellinfo$colsums <- colSums(Q$W[idx,idx])
cellinfo <- cellinfo %>% filter(colsums > 0)

cellinfo <- cellinfo %>% head(446)
idx <- cellinfo$Row_W
n_mc = sum(cellinfo$`cell type`==1)
ecov <- cellinfo %>% select(x,y,z) %>% rdist %>% squareform
ecov <- ecov/mean(ecov)
binarized_connections <- Q$W[idx,idx]>0
connectome <- Q$W[idx,idx]

cellinfo <- cellinfo %>% mutate(len = Length, cell_type = `cell type`) %>% head(446)
# build g
g <- network(binarized_connections[1:446, 1:446], directed = TRUE, vertex.attr = cellinfo[1:446,])