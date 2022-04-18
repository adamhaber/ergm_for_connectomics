library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(rdist)

connectome <- read_xlsx("../data/worm/SI 5 Connectome adjacency matrices, corrected July 2020.xlsx",
                        sheet = "connectome") %>% as.matrix()
rownames(connectome) <- connectome[,1]
connectome <- connectome[,2:ncol(connectome)]
post_synaptic <- colnames(connectome)
pre_synaptic <- rownames(connectome)
connectome[is.na(connectome)] <- 0
connectome <- matrix(as.numeric(connectome),    # Convert to numeric matrix
                  ncol = ncol(connectome))
(post_synaptic %in% pre_synaptic) %>% length
connectome <- connectome[,(post_synaptic %in% pre_synaptic)]
rownames(connectome) <- pre_synaptic
colnames(connectome) <- post_synaptic[(post_synaptic %in% pre_synaptic)]

cellinfo <- read_xlsx("../data/worm/SI 4 Cell lists.xlsx", sheet=2) 
names(cellinfo) <- c("name","type","subtype",
                  "notes, sensilla, modality", "other")
cellinfo <- cellinfo %>% filter(name %in% rownames(connectome)) 

# taken from https://www.wormatlas.org/neuronalwiring.html
pos <- read_xls("../data/worm/celegans_neuron_positions.xls")
pos <- pos %>% 
  mutate(name=Neuron, ganglion = `AY Ganglion Designation`,
         pos = `Soma Position`, region = `Soma Region`) %>% 
  select(name, pos, region, Span, ganglion) %>% 
  filter(name %in% cellinfo$name)
cellinfo <- left_join(cellinfo, pos, by = c("name" = "name"))

birth <- read_xls("../data/worm/time_of_birth.xls") %>% 
  filter(Neuron %in% cellinfo$name) %>% 
  mutate(time = `Birth time (min.)`) %>% select(-`Birth time (min.)`)
cellinfo <- left_join(cellinfo, birth, by = c("name" = "Neuron")) 

rich_club_neurons <- c("AVAR", "AVBR", "AVDR", "AVER", "PVCR", "AVAL", "AVBL", "AVDL", "AVEL", "PVCL","DVA")
cellinfo <- cellinfo %>% mutate(rich_club = name %in% rich_club_neurons) 

cellinfo <- cellinfo %>% mutate(first_two = str_sub(name, 1, 2)) %>% 
  mutate(subtype_ext = if_else((first_two %in% c("VD", "VA", "VB", "DD", "DA", "DB", "AS"))&(type=="motorneuron"), first_two, subtype))

connectome <- connectome[cellinfo$name, cellinfo$name]

bin_connectome <- (connectome > 0)
g <- network(bin_connectome, vertex.attr = cellinfo)



