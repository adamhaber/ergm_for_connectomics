setwd("./celegans/")
library(stringr)
library(ROCR)
library(dplyr)
library(ggplot2)
filenames <- Sys.glob("./celegans/results/10_3_no_subtype_ext/res*.rds")
current_fnames <- filenames[1:length(filenames)]
res <- lapply(current_fnames, tryCatch(readRDS, error = function(e) "error"))
names(res) <- sapply(current_fnames, function(x) (substr(basename(x), 5, nchar(x))))

bics <- sapply(res, function(x) x$bic[1])
aics <- sapply(res, function(x) x$aic[1])
lls <- sapply(res, function(x) x$ll)
had_nans <- unlist(sapply(res, function(x) tryCatch(x$had_nans, error = function(e) -1)), recursive = T)
had_nans = c()
for (r in res) {
  if (is.null(r$had_nans)) {
    had_nans <- c(had_nans, -1) 
  } else {
    had_nans <- c(had_nans, r$had_nans) 
  }
}
# auc_from_g <- sapply(res, function(x) x$aucs$from_g)
auc_from_er <- sapply(res, function(x) x$aucs$from_er)
auc_from_cond <- sapply(res, function(x) x$aucs$conditional)
# pr_from_g <- sapply(res, function(x) x$prs$from_g)
pr_from_er <- sapply(res, function(x) x$prs$from_er)
pr_from_cond <- sapply(res, function(x) x$prs$conditional)
# ham_from_g <- sapply(res, function(x) mean(x$hammings$from_g))
ham_from_er <- sapply(res, function(x) mean(x$hammings$from_er))
name <- sapply(res, function(x) x$name)
job_id <- sapply(res, function(x) x$job_id)
n_features <- sapply(res, function(x) x$n_features)
indep <- sapply(res, function(x) x$indep)
n_params <- sapply(res, function(x) x$model$coefficients %>% length)
terminations <- sapply(res, function(x) x$termination)
ESS <- sapply(res, function(x) x$effectiveSampleSize)

nanmean <- function(x) mean(na.omit(x))
nanmed <- function(x) median(na.omit(x))

mean_train <- sapply(res, function(r) r$CV[,1] %>% nanmean)
mean_test <- sapply(res, function(r) r$CV[,2] %>% nanmean)
mean_test_w_train_mle <- sapply(res, function(r) r$CV[,3] %>% nanmean)

median_train <- sapply(res, function(r) r$CV[,1] %>% nanmed)
median_test <- sapply(res, function(r) r$CV[,2] %>% nanmed)
median_test_w_train_mle <- sapply(res, function(r) r$CV[,3] %>% nanmed)

mean_train_perm <- sapply(res, function(r) r$CV_perm[,1] %>% nanmean)
mean_test_perm <- sapply(res, function(r) r$CV_perm[,2] %>% nanmean)
mean_test_w_train_mle_perm <- sapply(res, function(r) r$CV_perm[,3] %>% nanmean)

median_train_perm <- sapply(res, function(r) r$CV_perm[,1] %>% nanmed)
median_test_perm <- sapply(res, function(r) r$CV_perm[,2] %>% nanmed)
median_test_w_train_mle_perm <- sapply(res, function(r) r$CV_perm[,3] %>% nanmed)

stderr <- function(x, na.rm=TRUE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

se_train <- sapply(res, function(r) r$CV[,1] %>% stderr)
se_test <- sapply(res, function(r) r$CV[,2] %>% stderr)
se_test_w_train_mle <- sapply(res, function(r) r$CV[,3] %>% stderr)

se_train_perm <- sapply(res, function(r) r$CV_perm[,1] %>% stderr)
se_test_perm <- sapply(res, function(r) r$CV_perm[,2] %>% stderr)
se_test_w_train_mle_perm <- sapply(res, function(r) r$CV_perm[,3] %>% stderr)

cv_auc <- sapply(res, function(r) r$CV_perm[,4] %>% nanmean)
cv_auc_err <- sapply(res, function(r) r$CV_perm[,4] %>% stderr)



# activity_mean <- sapply(res, function(x) mean(x$activity[2,]))
# activity_sd <- sapply(res, function(x) sd(x$activity[2,]))
# activity_rnn_rmse <- sapply(res, function(x) mean(x$RNN_activity_rmse))
# activity_rnn_corr <- sapply(res, function(x) mean(x$RNN_activity_cor, na.rm=T))
g <- res[[1]]$g
gs_motifs <- summary(g~triadcensus)
gs_ideg <- colSums(as.matrix(g))
gs_odeg <- rowSums(as.matrix(g))
iterations <- sapply(res, function(x) x$model$iterations)


motifs_accuracy <- sapply(res, function(x) rowMeans(x$motifs/gs_motifs)-1) %>% na.omit %>% abs %>% colMeans
ideg_accuracy <- sapply(res, function(x) tryCatch(cor(gs_ideg, rowMeans(x$idegs)), error = function(e) -1)) 
odeg_accuracy <- sapply(res, function(x) tryCatch(cor(gs_odeg, rowMeans(x$odegs)), error = function(e) -1))

source("../utils/utils.R")
source("./celegans/load_worm_data.R")
VAs <- cellinfo %>% filter(str_detect(name, "^VA")) %>% select(name) %>% unlist
VBs <- cellinfo %>% filter(str_detect(name, "^VB")) %>% select(name) %>% unlist
DAs <- cellinfo %>% filter(str_detect(name, "^DA")) %>% select(name) %>% unlist
DBs <- cellinfo %>% filter(str_detect(name, "^DB")) %>% select(name) %>% unlist
gals_neurons <- c("ASHR", "AVAR", "AVBR", "AVDR", "PVCR", "PHBR", "PHAR","ASHL", "AVAL", "AVBL", "AVDL", "PVCL", "PHBL", "PHAL",
                  VAs, VBs, DAs, DBs)
gals_neurons <- sort(gals_neurons)
gals_idx = which(cellinfo$name %in% gals_neurons)

gals_neurons_no_motor <- c("ASHR", "AVAR", "AVBR", "AVDR", "PVCR", "PHBR", "PHAR","ASHL", "AVAL", "AVBL", "AVDL", "PVCL", "PHBL", "PHAL")
gals_neurons_no_motor <- sort(gals_neurons_no_motor)
gals_idx_no_motor = which(cellinfo$name %in% gals_neurons_no_motor)

gals_data <- non_diag((g %>% as.matrix)[gals_idx, gals_idx])
auc_on_gals_circuit <- function(m) {
  predmat <- m$prediction_matrices$from_er
  rownames(predmat) <- colnames(predmat) <- cellinfo$name
  pred <- non_diag(predmat[gals_idx, gals_idx])
  performance(prediction(pred, gals_data), "auc")@y.values[[1]]
}
# gal_auc <- sapply(res, auc_on_gals_circuit)

df <- data.frame(filenames=current_fnames,
                 bics = bics %>% unlist(recursive = T),
                 aics=aics,
                 lls=lls,
                 # iter = iterations,
                 n_params = n_params,
                 # n_features=n_features,
                 auc_from_cond=auc_from_cond, auc_from_er=auc_from_er,
                 pr_from_cond=pr_from_cond, pr_from_er=pr_from_er,
                 ham_from_er=ham_from_er, indep=indep,
                 name=name, job_id=job_id,
                 had_nans=had_nans,
                 # ess=ESS,
                 # terminations=terminations)
                 mean_train = mean_train,
                 mean_test = mean_test,
                 mean_test_w_train_mle = mean_test_w_train_mle,
                 # gal_auc = gal_auc,
                 mean_train_perm = mean_train_perm,
                 mean_test_perm = mean_test_perm,
                 mean_test_w_train_mle_perm = mean_test_w_train_mle_perm,
                 
                 median_train_perm = median_train_perm,
                 median_test_perm = median_test_perm,
                 median_test_w_train_mle_perm = median_test_w_train_mle_perm,
                 
                 se_train = se_train,
                 se_test = se_test,
                 se_test_w_train_mle = se_test_w_train_mle,
                 
                 se_train_perm = se_train_perm,
                 se_test_perm = se_test_perm,
                 se_test_w_train_mle_perm = se_test_w_train_mle_perm,
                 cv_auc = cv_auc,
                 
                 # activity_mean=activity_mean,
                 # activity_sd=activity_sd,
                 
                 is_mutual=sapply(current_fnames, function(x) str_detect(x, "mutual")),
                 is_ext=sapply(current_fnames, function(x) str_detect(x, "ext")),
                 is_absdiff= sapply(current_fnames, function(x) str_detect(x, "absdiff")),
                  is_rich=sapply(current_fnames, function(x) str_detect(x, "rich")),
                 #  is_esp=sapply(current_fnames, function(x) str_detect(x, "esp")),
                 # is_type=sapply(current_fnames, function(x) str_detect(x, "type")),
                 # activity_rnn_corr=activity_rnn_corr,
                 # activity_rnn_rmse=activity_rnn_rmse,
                 motifs_accuracy=motifs_accuracy,
                 ideg_accuracy=ideg_accuracy, odeg_accuracy=odeg_accuracy
)


##########

df <- df %>% arrange(job_id)


forms <- list(is_type = c(0,1,1),
                 is_ganglion = c(0,1,1),
                 is_region = c(0,1),
                 is_span = c(0,1),
                 is_rich_club = c(0,1),
                 is_time = c(0,1),
                 is_pos = c(0,1),
                 is_distance = c(0,1),
                 is_reciprocity = c(0,1))


df <- df %>% filter(!str_detect(name, "sender"))
terms_df <- do.call(expand.grid, forms)
terms_df <- terms_df %>% mutate(form_id = 1:n())
terms_df <- terms_df[terms_df$form_id %in% df$job_id,] %>% select(-form_id)
n_features <- 1+rowSums(terms_df)
df <- cbind(df, terms_df, n_features)



#########

ggplot(df, aes(ymin=mean_test_w_train_mle_perm-se_test_w_train_mle_perm, 
               y=mean_test_w_train_mle_perm,
               ymax=mean_test_w_train_mle_perm+se_test_w_train_mle_perm, 
               x=n_features)) +
  geom_pointrange(position = position_jitter(w = 0.15, h = 0), size=0.2) +
  theme_classic() + xlab("Number of biological features") +
  ylab("Cross-validated LL") +
  scale_x_continuous(breaks = 0:10) + theme(legend.position='none') -> p_n_features_ll
# geom_label_repel(data = subset(df %>% filter(n_features>0) %>% filter(job_id==412)), aes(label="Model"),
# nudge_x = -2,
# nudge_y = 200)
# geom_point(data = tmp %>% filter(job_id==412), aes(x=n_features, y=mean_test_w_train_mle_perm), col="black", size=3) -> p
p_n_features_ll

library(ggallin)
ggplot(df %>% filter(lls>quantile(df$lls, .05)), aes(x=lls,
                                       # xmin=mean_test_w_train_mle_perm-se_test_w_train_mle_perm,
                                       # xmax=mean_test_w_train_mle_perm+se_test_w_train_mle_perm,
                                       y=auc_from_er, 
                                       # ymin=cv_auc-cv_auc_err,
                                       # ymax=cv_auc+cv_auc_err, 
                                       color=as.factor(n_features))) +
  geom_point() +
  # geom_errorbar() +
  # geom_errorbarh() +
  # geom_pointrangeh(position = position_jitter(w = 0, h = 0), size=0.3) +
  theme_classic() + 
  scale_x_continuous(trans = pseudolog10_trans) +
  xlab("Cross-validated LL") +
  # xlim(c(-4.5,-4)) +
  # scale_x_log10() +
  ylab("Cross-validated AUROC") +
  scale_color_brewer(palette = "Blues", direction = -1) +
  ylab("Cross-validated LL") + 
  theme(legend.position = 'none') -> p_ll_auc_n_features_xerr
p_ll_auc_n_features_xerr
# ggsave("./individual_figures/worm/pdfs/p_ll_auc_n_features_xerr.pdf", plot=p_ll_auc_n_features_xerr)
# ggsave("./individual_figures/worm/pngs/p_ll_auc_n_features_xerr.png", plot=p_ll_auc_n_features_xerr)

# ggplot(df, aes(x=mean_test_w_train_mle_perm,
#                                        y=cv_auc, 
#                                        ymin=cv_auc-cv_auc_err,
#                                        ymax=cv_auc+cv_auc_err, 
#                                        color=as.factor(n_features))) +
#   geom_point() +
#   geom_errorbar() +
#   # geom_errorbarh() +
#   # geom_pointrangeh(position = position_jitter(w = 0, h = 0), size=0.3) +
#   theme_classic() + 
#   xlab("Cross-validated LL") +
#   ylab("Cross-validated AUROC") +
#   scale_color_brewer(palette = "Blues", direction = -1) +
#   ylab("Cross-validated LL") -> p_ll_auc_n_features
# p_ll_auc_n_features
# ggsave("./individual_figures/worm/p_ll_auc_n_features.pdf", plot=p_ll_auc_n_features)


tmp = df %>% filter(is_mutual==T, is_type==T)
plot(tmp$auc_from_er, tmp$lls)
q = 0.9
top_ll <- quantile(tmp$lls,  q)
top_auc <- quantile(tmp$auc_from_er,  q)
best_m_names <- tmp %>% filter(lls>top_ll) %>% 
  filter(auc_from_er>top_auc) %>% 
  arrange(n_features, n_params) %>% 
  select(n_features, n_params, name)
best_m_names %>% select(name) %>% unlist %>% first

best_m <- res$`939_edges+nmix(~subtype)+nicov(~poly(time,2))+nocov(~poly(time,2))+absdiff(~pos)+mutual.rds`
samples <- readRDS("./celegans/results/10_3_no_subtype_ext/samples_939_edges+nmix(~subtype)+nicov(~poly(time,2))+nocov(~poly(time,2))+absdiff(~pos)+mutual.rds")
er_model <- res$`1_edges.rds`
er_samples <- readRDS("./celegans/results/10_3_no_subtype_ext/samples_1_edges.rds")

df %>% filter(str_detect(name, "nmix\\(~subtype")) %>%
  filter(str_detect(name, "nocov\\(~poly\\(time,2\\)")) %>% 
  filter(str_detect(name, "mutual")) %>% 
  filter(str_detect(name, "absdiff\\(~pos")) %>% 
  arrange(-n_features) %>% 
  select(job_id) %>% head %>% unlist

full_model <- res$`1113_edges+nmix(~subtype)+nifactor(~ganglion)+nofactor(~ganglion)+nmix(~region)+nifactor(~Span)+nofactor(~Span)+nicov(~poly(time,2))+nocov(~poly(time,2))+nicov(~poly(pos,2))+nocov(~poly(pos,2))+absdiff(~pos)+mutual.rds`

animal = "worm"
animal_path = "worm"
cmap = "darkblue"
ncuts = 3

# highest_rocs_per_n_features = list()
# roc_objs <- list()
# library(rlist)
# 
# for (n in 1:max(df$n_features)) {
#   tmp = df %>% filter(n_features==n) %>% arrange(auc_from_er) %>% tail(1) %>% select(job_id) %>% unlist
#   highest_rocs_per_n_features[[n]] = list.filter(res, job_id==tmp)
#   roc_objs[[n]] <- roc(non_diag(highest_rocs_per_n_features[[n]][[1]]$g %>% as.matrix), 
#                        non_diag(highest_rocs_per_n_features[[n]][[1]]$prediction_matrices$from_er %>% as.matrix))
# }
# roc_objs[[(1+length(roc_objs))]] <- roc(non_diag(best_m$g %>% as.matrix), non_diag(best_m$prediction_matrices$from_er %>% as.matrix))
# 
# library(RColorBrewer)
# 
# pal <- brewer.pal(10, "Blues")
# pal[[length(roc_objs)]] <- "black"
# p_rocs_all <- ggroc(roc_objs) + theme_classic()  + coord_fixed() + scale_color_manual(values = pal)
# p_rocs_all

weight_breaks = c(-0.1, 1:4, Inf)
weight_breaks
Ng = network.size(g)
degree_dist_upper_bound = 60
source("../utils/utils.R")
source("./celegans/load_worm_data.R")

save.image(file = "./worm_workspace_with_res.RData")

# df_orig <- df

res = NULL
conf_samples = NULL
# save.image(file = "./worm_workspace.RData")

source("./gen_figures_for_best_model.R")
save.image(file = "./worm_figures.RData")










stop()
############### old

ggplot(df, aes(x=auc_from_er, y=cv_auc, color=as.factor(is_mutual))) + geom_point(position = position_jitter(w = 0.2, h = 0)) +
  geom_abline(intercept=0, slope=1)+
  theme_classic()

# mutuality is critical for motifs
ggplot(df, aes(x=motifs_accuracy, color=is_mutual)) + geom_density() + xlim(0,4) +
  theme_classic()

ggplot(df %>% filter(is_mutual==T), aes(x=n_features, y=cv_auc, color=motifs_accuracy)) + geom_point(position = position_jitter(w = 0.2, h = 0)) +
  theme_classic()

ggplot(df %>% filter(is_mutual==T), aes(x=n_params, y=cv_auc, color=motifs_accuracy)) + geom_point(position = position_jitter(w = 0.2, h = 0)) +
  theme_classic()

# df %>% filter(is_mutual==T) %>% View

ggplot(df, aes(x=n_params, y=cv_auc, color=as.factor(is_mutual))) + geom_point(position = position_jitter(w = 0.2, h = 0)) +
  theme_classic()

# library(plotly)
# p <- ggplot(df %>% filter(is_ext==T) %>% filter(is_mutual==T) %>% filter(is_rich==T), aes(x=n_features, y=cv_auc, color=motifs_accuracy,
#                     text = paste("names:",job_id))) + geom_point(position = position_jitter(w = 0.2, h = 0)) +
#   theme_classic()
# 
# fig <- ggplotly(p)
# 
# fig

# df %>% filter(is_mutual==F) %>% View

# this one really misses the degree distribution
# best_m <- res$`1156_edges+nmix(~subtype_ext)+absdiff(~pos)+mutual.rds`
# samples <- readRDS("./celegans/results/18_11_w_samples/samples_1156_edges+nmix(~subtype_ext)+absdiff(~pos)+mutual.rds")

ggplot(df %>% filter(is_ext==T), aes(x=gal_auc, fill=is_mutual)) + geom_density()
# View(df %>% filter(is_ext==T) %>% filter(is_mutual==F) %>% filter(is_rich==T) %>% arrange(gal_auc))

#best model
# best_m <- res$`403_edges+nmix(~subtype)+nifactor(~ganglion)+nofactor(~ganglion)+nmix(~region)+absdiff(~pos).rds`
# samples <- readRDS("./celegans/results/18_11_w_samples/samples_403_edges+nmix(~subtype)+nifactor(~ganglion)+nofactor(~ganglion)+nmix(~region)+absdiff(~pos).rds")
# model <- readRDS("./celegans/results/18_11_w_samples/model_403_edges+nmix(~subtype)+nifactor(~ganglion)+nofactor(~ganglion)+nmix(~region)+absdiff(~pos).rds")

# best_m <- res$`448_edges+nmix(~subtype_ext)+nmix(~region)+nmix(~rich_club)+absdiff(~pos).rds`
# samples <- readRDS("./celegans/results/18_11_w_samples/samples_448_edges+nmix(~subtype_ext)+nmix(~region)+nmix(~rich_club)+absdiff(~pos).rds")
# model <- readRDS("./celegans/results/18_11_w_samples/model_448_edges+nmix(~subtype_ext)+nmix(~region)+nmix(~rich_club)+absdiff(~pos).rds")

# best_m <- res$`1503_edges+nmix(~subtype)+nmix(~region)+nmix(~rich_club)+nicov(~poly(time,2))+nocov(~poly(time,2))+nicov(~poly(pos,2))+nocov(~poly(pos,2))+absdiff(~pos)+mutual.rds`
# samples <- readRDS("./celegans/results/18_11_w_samples/samples_1503_edges+nmix(~subtype)+nmix(~region)+nmix(~rich_club)+nicov(~poly(time,2))+nocov(~poly(time,2))+nicov(~poly(pos,2))+nocov(~poly(pos,2))+absdiff(~pos)+mutual.rds")

animal = "worm"
animal_path = "worm"
cmap = "darkblue"
ncuts = 3

# connectome %>% as.vector %>% hist(breaks=100)
# connectome[connectome>0] %>% as.vector %>% table
weight_breaks = c(-0.1, 1:4, Inf)
weight_breaks
Ng = network.size(g)
degree_dist_upper_bound = 60
source("../utils/utils.R")
source("./celegans/load_worm_data.R")
er_model <- res$`1_edges.rds`
er_samples <- readRDS("./celegans/results/18_11_w_samples/samples_1_edges.rds")

#gal's circuit
create_heatmap(g[gals_idx_no_motor, gals_idx_no_motor], cmap=cmap) -> p_gal_data
create_heatmap(best_m$prediction_matrices$from_er[gals_idx_no_motor, gals_idx_no_motor], cmap=cmap) -> p_gal_pred
create_heatmap(best_m$prediction_matrices$conditional[gals_idx_no_motor, gals_idx_no_motor], cmap=cmap) -> p_gal_cond

library(GGally)
ggnet2(g[gals_idx, gals_idx], size = 3, color = "black", edge.size = 0.2, edge.color = "grey", arrow.size = 2, arrow.gap = 0.015,
       mode = "circle", label=g[gals_idx, gals_idx] %>% rownames, label.size=3, vjust = -1) -> p_gals_graph
# ggsave("./final figures/worm/gals_neurons_graph.pdf", plot=p_gals_graph)

gal_cmap = "mediumvioletred"
create_heatmap(g[gals_idx, gals_idx], cmap=gal_cmap) + ggtitle("Data") -> p_gal_data
create_heatmap(best_m$prediction_matrices$from_er[gals_idx, gals_idx], cmap=gal_cmap) + ggtitle("Mean") -> p_gal_pred
create_heatmap(best_m$closest$from_er[gals_idx, gals_idx]*1, cmap=gal_cmap) + ggtitle("Sample") -> p_gal_sample
grid.arrange(p_gal_data, p_gal_pred, p_gal_sample, nrow=1) -> p_gal_heatmaps

# ggsave("./final figures/worm/gals_neurons.pdf", plot=p_gal_heatmaps)
# ggsave("./final figures/worm/gals_neurons.png", plot=p_gal_heatmaps)


#### study mismatches
neuron_aucs_out <- sapply(1:Ng, function(i) {
  p <- best_m$prediction_matrices$from_er[i,]
  return(tryCatch(performance(prediction(p, g[i,]), "auc")@y.values[[1]],
                  error = function(e) NA))})

neuron_aucs_in <- sapply(1:Ng, function(i) {
  p <- best_m$prediction_matrices$from_er[,i]
  return(tryCatch(performance(prediction(p, g[,i]), "auc")@y.values[[1]],
                  error = function(e) NA))})

neuron_ents_out <- sapply(1:Ng, function(i) {
  p <- best_m$prediction_matrices$from_er[i,]
  return(-sum(p*log(p)+(1-p)*log(1-p),na.rm = T)/(length(p)-1))})

neuron_ents_in <- sapply(1:Ng, function(i) {
  p <- best_m$prediction_matrices$from_er[,i]
  return(-sum(p*log(p)+(1-p)*log(1-p),na.rm = T)/(length(p)-1))})

cellinfo_aug <- cellinfo %>% 
  mutate(Output_entropy = neuron_ents_out,
         Input_entropy = neuron_ents_in,
         Output_AUC = neuron_aucs_out,
         Input_AUC = neuron_aucs_in, 
         gen_auc = (neuron_aucs_out+neuron_aucs_in)/2,
         Indegree = colSums(g %>% as.matrix),
         Outdegree = rowSums(g %>% as.matrix))
cellinfo_aug %>% select(Indegree, Outdegree, Input_AUC, Output_AUC, Input_entropy, Output_entropy) %>% na.omit %>% cor

library(GGally)
my_fn <- function(data, mapping, method="loess", ...){
  p <- ggplot(data = data, mapping = mapping) + 
    geom_point() + 
    geom_smooth(method=method, ...)
  p
}

# Default loess curve    
# cellinfo_aug %>%
#   select(Indegree, Outdegree, Input_AUC, Output_AUC, Input_entropy, Output_entropy) %>%
#   na.omit %>%
#   ggpairs(lower = list(continuous = my_fn)) +
#   theme_classic() -> p
# p
# ggsave("./final figures/worm/worm_entropies_aucs_degrees_pairplot.pdf", plot = p)
# ggsave("./final figures/worm/worm_entropies_aucs_degrees_pairplot.png", plot = p)

# p1 <- ggplot(cellinfo_aug, aes(y=gen_auc, x=time)) + geom_point() + theme_classic()
# p1
# p2 <- ggplot(cellinfo_aug, aes(y=Output_entropy, color=time)) + geom_boxplot() + 

# ggplot(cellinfo_aug, aes(y=gen_auc, x=type)) + geom_boxplot() + theme_classic(base_size=15)+ ylim(0,1)
# ggsave("./final figures/worm/worm_input_AUC_per_type.png")
# ggplot(cellinfo_aug, aes(y=Input_AUC, x=ganglion)) + geom_boxplot() + theme_classic()+ ylim(0,1)
# ggsave("./final figures/worm/worm_input_AUC_per_ganglion.png")

ggplot(cellinfo_aug, aes(y=gen_auc, x=pos)) + 
  geom_point() + 
  geom_smooth(data=cellinfo_aug, aes(y=Input_AUC, x=pos)) +
  theme_classic(base_size = 15) + ylim(0,1) + xlab("A-P position") +
  ylab("AUC") -> p_auc_position
# ggsave("./final figures/worm/worm_input_AUC_per_pos.png", plot=p_auc_position)
# ggsave("./final figures/worm/worm_input_AUC_per_pos.pdf", plot=p_auc_position)

p_auc_types <- ggplot(cellinfo_aug %>% mutate(Type=type), aes(y=Input_AUC, x=Output_AUC, col=Type)) + 
  geom_point() + 
  theme_classic(base_size = 15) + xlim(0.5, 1) +
  ylim(0.5, 1) +
  geom_abline(intercept=0, slope=1) +
  coord_fixed() 
# ggsave("./final figures/worm/worm_AUC_per_type.png", plot=p_auc_types)
# ggsave("./final figures/worm/worm_AUC_per_type.pdf", plot=p_auc_types)

res = NULL
conf_samples = NULL
save.image(file = "./worm_workspace.RData")

source("./gen_figures_for_best_model.R")
save.image(file = "./worm_figures.RData")
stop()




library(plotly)
p_auc_types<-ggplot(cellinfo_aug, aes(y=Input_AUC, x=Output_AUC, col=type, text=name)) + 
  geom_point() + 
  theme_classic(base_size = 15) + xlim(0.5, 1) +
  ylim(0.5, 1) +
  geom_abline(intercept=0, slope=1) +
  coord_fixed() 

fig <- ggplotly(p)

fig



cellinfo_aug %>% names

# sapply(best_m$odegs, function(x) unlist(x)[1:271]) %>% rowMeans %>% hist(breaks=100, xlim=c(0,50))
# g %>% as.matrix %>% rowSums %>% hist(breaks=100)
# 
# p <- do.call("grid.arrange", c(list(p1,p2,p3,p4), ncol=2))
# ggsave("./final figures/worm/entropies_aucs_boxplot.pdf", plot=p)
# ggsave("./final figures/worm/entropies_aucs_boxplot.png", plot=p)

res = NULL
conf_samples = NULL
save.image(file = "./worm_workspace.RData")
stop()








# ###### old stuff
# ggplot(df, aes(x=n_params, y=lls, color=as.factor(is_mutual))) + geom_point(position = position_jitter(w = 0, h = 0)) +
#   theme_class
# 
# ggplot(df, aes(x=n_params, y=motifs_accuracy, color=is_mutual)) + geom_point(position = position_jitter(w = 0, h = 0)) +
#   theme_classic()
# 
# ggplot(df, aes(x=n_params, y=auc_from_er, color=is_mutual)) + geom_point(position = position_jitter(w = 0, h = 0)) +
#   theme_classic()
# 
# ggplot(df, aes(x=n_params, y=auc_from_er, color=is_mutual)) + geom_point(position = position_jitter(w = 0, h = 0)) +
#   theme_classic()
# 
# # we want subtype_ext for high AUCs
# ggplot(df, aes(x=n_params, y=gal_auc, color=is_mutual)) + geom_point(position = position_jitter(w = 0, h = 0)) +
#   theme_classic()
# 
# #mutuality is important for the motifs accuracy
# ggplot(df, aes(x=motifs_accuracy, color=is_mutual)) + geom_density() +
#   theme_classic()
# 
# #mutuality is important for the motifs accuracy
# ggplot(df, aes(x=auc_from_er, y=gal_auc, color=n_params)) + geom_point(position = position_jitter(w = 0, h = 0)) +
#   theme_classic()
# 
# 
# ggplot(df %>% filter(is.finite(median_test_w_train_mle_perm)) %>% 
#          filter(median_test_w_train_mle_perm > -1e20), aes(x=n_params, y=gal_auc, color=as.factor(is_absdiff))) + geom_point(position = position_jitter(w = 0., h = 0)) +
#   theme_classic()
# 
# 
# library(plotly)
# p<-ggplot(df %>% filter(n_params<500) %>% filter(motifs_accuracy<2), aes(x=n_params, y=ideg_accuracy, color=motifs_accuracy, text=job_id)) + geom_point(position = position_jitter(w = 0., h = 0)) +
#   theme_classic()
# 
# fig <- ggplotly(p)
# 
# fig
# 
# df %>% filter(n_params<250) %>% View
# 
# ggplot(df %>% filter(lls>-5e4), aes(x=n_params, y=odeg_accuracy, color=as.factor(is_mutual))) + geom_point(position = position_jitter(w = 0.2, h = 0)) +
#   theme_classic()
# 
# ggplot(df %>% filter(lls>-5e4), aes(x=n_params, y=ideg_accuracy, color=as.factor(is_mutual))) + geom_point(position = position_jitter(w = 0.2, h = 0)) +
#   theme_classic()
# 
# ggplot(df %>% filter(lls>-5e4), aes(x=n_params, y=motifs_accuracy, color=as.factor(is_mutual))) + geom_point(position = position_jitter(w = 0.2, h = 0)) +
#   theme_classic()
# 
# ggplot(df %>% filter(lls>-5e4), aes(x=n_params, y=lls, color=as.factor(is_mutual))) + geom_point(position = position_jitter(w = 0.2, h = 0)) +
#   theme_classic()
# 
# ggplot(df %>% filter(lls>-5e4), aes(x=n_params, y=aics, color=as.factor(is_mutual))) + geom_point(position = position_jitter(w = 0.2, h = 0)) +
#   theme_classic()
# 
# ggplot(df %>% filter(lls>-5e4), aes(x=n_params, y=ham_from_er, color=as.factor(is_mutual))) + geom_point(position = position_jitter(w = 0.2, h = 0)) +
#   theme_classic()
# 
# ggplot(df %>% filter(lls>-5e4), aes(x=n_params, y=auc_from_cond, color=as.factor(is_mutual))) + geom_point(position = position_jitter(w = 0.2, h = 0)) +
#   theme_classic()
# 
# ggsave("./final figures/worm/worm_indep_models_ll_vs_n_param.png")
# 
# ggplot(df %>% filter(lls>-5e4), aes(x=n_params, y=auc_from_er, color=as.factor(is_mutual))) + geom_point(position = position_jitter(w = 0.2, h = 0)) +
#   theme_classic()
# ggsave("./final figures/worm/worm_indep_models_ll_vs_n_param.png")
# 
# 
# 
# source("../utils/utils.R")
# par(mfrow=c(3,3))
# single_features = df %>% filter(n_features<=1) %>% arrange(lls) %>% rownames
# ps <- list()
# for (f in single_features) {
#   ps[[f]] <- create_heatmap(res[[f]]$prediction_matrices$from_er, cmap,
#                             title=strsplit(res[[f]]$name, "_")[[1]][2])
#   # image(res[[f]]$prediction_matrices$from_er, main=strsplit(res[[f]]$name, "_")[[1]][2])
# }
# p_all <- do.call("grid.arrange", c(ps[1:9], ncol=3))
# ggsave("./final figures/worm/worm_all_single_feature_models.png", p_all)
# ggsave("./final figures/worm/worm_all_single_feature_models.pdf", p_all)
# create_heatmap(g %>% as.matrix, cmap, title="Data")
# ggsave("./final figures/worm/worm_data.png")
# 
# ### best model
# # best_m <- res$`192_edges+nmix(~subtype)+nmix(~ganglion)+nmix(~region)+nifactor(~Span)+nofactor(~Span)+nmix(~rich_club)+nicov(~poly(time,2))+nocov(~poly(time,2))+absdiff(~pos).rds`
# # best_m <- res$`380_edges+nmix(~type)+nmix(~region)+nifactor(~Span)+nofactor(~Span)+nmix(~rich_club)+nicov(~poly(time,2))+nocov(~poly(time,2))+absdiff(~pos)+nicov(~poly(pos,2))+nocov(~poly(pos,2)).rds`
# # best_m <- res$`465_edges+nmix(~subtype)+nmix(~region)+nmix(~rich_club)+nicov(~poly(time,2))+nocov(~poly(time,2))+mutual.rds`
# # best_m <- res$`555_edges+nmix(~subtype)+nmix(~rich_club)+nicov(~poly(time,2))+nocov(~poly(time,2))+absdiff(~pos)+mutual.rds`
# # best_m <- res$`832_edges+nmix(~subtype_ext)+nmix(~region)+nmix(~rich_club)+mutual.rds`
# best_m <- res$`440_edges+nmix(~subtype_ext)+nifactor(~ganglion)+nofactor(~ganglion)+nmix(~rich_club)+absdiff(~pos).rds`
# # best_m <- res$`1348_edges+nmix(~subtype_ext)+nicov(~poly(pos,2))+nocov(~poly(pos,2))+absdiff(~pos)+mutual.rds`
# # best_m <- res$`100_edges+nmix(~ganglion)+absdiff(~pos).rds`
# animal = "worm"
# cmap = "darkblue"
# ncuts = 3
# 
# connectome %>% as.vector %>% unique
# weight_breaks = c(-0.1, 0.1, 1:4*3)
# weight_breaks
# Ng = network.size(g)
# degree_dist_upper_bound = 60
# 
# conf_m = res$`1537_sender+receiver.rds`
# source("./gen_figures_for_best_model.R")
# 
# 
# 
# 
# predmat <- best_m$prediction_matrices$conditional
# rownames(predmat) <- colnames(predmat) <- cellinfo$name
# create_heatmap(predmat[gals_neurons, gals_neurons], cmap)
# create_heatmap((g %>% as.matrix)[gals_neurons, gals_neurons], cmap, title="Data")
# 
# #conditional on the rest
# g_cond <- g
# idx = which(cellinfo$name %in% gals_neurons)
# g_cond[idx, idx] <- NA
# condsims <- simulate(best_m$model$formula %>% as.formula, coef = best_m$model$coefficients, basis=g_cond,
#                      constraints = ~observed, nsim = 100,
#                      control = control.simulate.formula(MCMC.effectiveSize = NULL, MCMC.interval = 1e6))
# condpreds <- lapply(condsims, function(nw) nw[idx,idx])
# condmean <- mean_connectivity(condsims)
# 
# create_heatmap(predmat[gals_neurons, gals_neurons], cmap)
# create_heatmap((g %>% as.matrix)[gals_neurons, gals_neurons], cmap)
# create_heatmap(predmat[gals_neurons,gals_neurons], cmap)
# create_heatmap(condmean[gals_neurons,gals_neurons], cmap)
# 
# 
# plot(predmat[idx,idx], condmean)
