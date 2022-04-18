setwd("./ob/")
library(stringr)
library(ROCR)
library(dplyr)
source("../utils/utils.R")
library(ggplot2)
base_path = "./ob/results/24_3_new_microglom/"
filenames <- Sys.glob("./ob/results/24_3_new_microglom/res*.rds")
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
# n_terms <- sapply(res, function(x) x$n_terms)
indep <- sapply(res, function(x) x$indep)
n_params <- sapply(res, function(x) x$model$coefficients %>% length)
terminations <- sapply(res, function(x) x$termination)
ESS <- sapply(res, function(x) x$effectiveSampleSize)

mean_train <- sapply(res, function(r) r$CV[,1] %>% mean)
mean_test <- sapply(res, function(r) r$CV[,2] %>% mean)
mean_test_w_train_mle <- sapply(res, function(r) r$CV[,3] %>% mean)

median_train <- sapply(res, function(r) r$CV[,1] %>% median)
median_test <- sapply(res, function(r) r$CV[,2] %>% median)
median_test_w_train_mle <- sapply(res, function(r) r$CV[,3] %>% median)

mean_train_perm <- sapply(res, function(r) r$CV_perm[,1] %>% mean)
mean_test_perm <- sapply(res, function(r) r$CV_perm[,2] %>% mean)
mean_test_w_train_mle_perm <- sapply(res, function(r) r$CV_perm[,3] %>% mean)

median_train_perm <- sapply(res, function(r) r$CV_perm[,1] %>% median)
median_test_perm <- sapply(res, function(r) r$CV_perm[,2] %>% median)
median_test_w_train_mle_perm <- sapply(res, function(r) r$CV_perm[,3] %>% median)

stderr <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}

se_train <- sapply(res, function(r) r$CV[,1] %>% stderr)
se_test <- sapply(res, function(r) r$CV[,2] %>% stderr)
se_test_w_train_mle <- sapply(res, function(r) r$CV[,3] %>% stderr)

se_train_perm <- sapply(res, function(r) r$CV_perm[,1] %>% stderr)
se_test_perm <- sapply(res, function(r) r$CV_perm[,2] %>% stderr)
se_test_w_train_mle_perm <- sapply(res, function(r) r$CV_perm[,3] %>% stderr)

cv_auc_perm <- sapply(res, function(r) r$CV_perm[,4] %>% mean)
cv_auc_perm_err <- sapply(res, function(r) r$CV_perm[,4] %>% stderr)
se_test_perm <- sapply(res, function(r) r$CV_perm[,2] %>% stderr)
se_test_w_train_mle_perm <- sapply(res, function(r) r$CV_perm[,3] %>% stderr)


activity_mean <- sapply(res, function(x) mean(x$activity[2,]))
activity_sd <- sapply(res, function(x) sd(x$activity[2,]))

# activity_rnn_rmse <- sapply(res, function(x) mean(x$RNN_activity_rmse))
# activity_rnn_corr <- sapply(res, function(x) mean(x$RNN_activity_cor, na.rm=T))
g <- res[[1]]$g
ob_g <- g
library(ergm)
gs_motifs <- summary(g ~ triadcensus)
gs_ideg <- colSums(as.matrix(g))
gs_odeg <- rowSums(as.matrix(g))
iterations <- sapply(res, function(x) x$model$iterations)

# motifs_accuracy <- sapply(res, function(x) sqrt(mean((log(gs_motifs)-log(rowMeans(x$motifs)))^2)))
motifs_accuracy <- sapply(res, function(x) rowMeans(x$motifs/gs_motifs)-1) %>% na.omit %>% abs %>% colMeans
ideg_accuracy <- sapply(res, function(x) tryCatch(cor(gs_ideg, rowMeans(x$idegs)), error = function(e) -1)) 
odeg_accuracy <- sapply(res, function(x) tryCatch(cor(gs_odeg, rowMeans(x$odegs)), error = function(e) -1))

df <- data.frame(filenames=current_fnames,
                 bics = bics %>% unlist(recursive = T),
                 aics=aics,
                 lls=lls,
                 iter = iterations,
                 n_params = n_params,
                 # n_terms=n_terms,
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
                 mean_train_perm = mean_train_perm,
                 mean_test_perm = mean_test_perm,
                 mean_test_w_train_mle_perm = mean_test_w_train_mle_perm,
                 
                 median_train_perm = median_train_perm,
                 median_test_perm = median_test_perm,
                 median_test_w_train_mle_perm = median_test_w_train_mle_perm,
                 
                 se_train = se_train,
                 se_test = se_test,
                 se_test_w_train_mle = se_test_w_train_mle,
                 cv_auc_perm = cv_auc_perm,
                 cv_auc_perm_err = cv_auc_perm_err,
                 
                 
                 se_train_perm = se_train_perm,
                 se_test_perm = se_test_perm,
                 se_test_w_train_mle_perm = se_test_w_train_mle_perm,
                 
                 activity_mean=activity_mean,
                 activity_sd=activity_sd,
                 
                 motifs_accuracy=motifs_accuracy,
                 ideg_accuracy=ideg_accuracy, odeg_accuracy=odeg_accuracy
                 )


##### 

df <- df %>% arrange(job_id)

forms <- list(is_mut = c(0, 1),
                 is_cell_type = c(0, 1),
                 is_distance = c(0, 1),
                 is_olen = c(0, 1),
                 is_ilen = c(0, 1),
                 is_o = c(0, 1),
                 is_i = c(0, 1),
                 is_glom = c(0, 1))

terms_df <- do.call(expand.grid, forms)
n_features <- 1+rowSums(terms_df)
df <- cbind(df, terms_df, n_features)
# df %>% select(name, n_features) %>% View
# df_full <- df %>% mutate(single_terms = (n_terms==1))
# df <- df_full %>% filter(is_sender==F)

#### anotate_single_feature_models 

# ordering to get the grouping of the channels-by-color right
single_feature_models <- (df %>% filter(n_features <= 2) %>% rownames)[c(1,2,5,7,6,8,4,3,9)]
single_feature_models
# df$single_feauture_models <- names(res) %in% single_feature_models
df$resnames <- rownames(df)
# df$colors <- "black"
library(RColorBrewer)
# df$colors[df$single_feauture_models==T] <- brewer.pal(name="Paired", n=9)

######################## n terms
ggplot(df, aes(ymin=mean_test_w_train_mle_perm-se_test_w_train_mle_perm, 
                                     y=mean_test_w_train_mle_perm,
                                     ymax=mean_test_w_train_mle_perm+se_test_w_train_mle_perm, 
                                     x=n_features)) +
  geom_pointrange(position = position_jitter(w = 0.15, h = 0), size=0.2) +
  theme_classic() + xlab("Number of biological features") +
  ylab("Cross-validated LL") +
  scale_x_continuous(breaks = 0:10) + theme(legend.position='none') -> p_n_features_ll
  # geom_label_repel(data = subset(df %>% filter(n_terms>0) %>% filter(job_id==412)), aes(label="Model"),
                   # nudge_x = -2,
                   # nudge_y = 200)
  # geom_point(data = tmp %>% filter(job_id==412), aes(x=n_terms, y=mean_test_w_train_mle_perm), col="black", size=3) -> p
p_n_features_ll
# # ggsave("./individual_figures/ob/p_cv_ll_features.pdf", plot=p_n_features_ll)
# # ggsave("./final figures/fig1e.png", plot=p_n_features_ll)
# # ggsave("./final figures/fig1e.pdf", plot=p_n_features_ll)

# auc vs ll vs n_terms
ggplot(df, aes(x=mean_test_w_train_mle_perm,
                                       xmin=mean_test_w_train_mle_perm-se_test_w_train_mle_perm,
                                       xmax=mean_test_w_train_mle_perm+se_test_w_train_mle_perm,
                                       y=cv_auc_perm, 
                                       ymin=cv_auc_perm-cv_auc_perm_err,
                                       ymax=cv_auc_perm+cv_auc_perm_err, 
                                       color=as.factor(n_features))) +
  geom_point() +
  geom_errorbar() +
  geom_errorbarh() +
  # geom_pointrangeh(position = position_jitter(w = 0, h = 0), size=0.3) +
  theme_classic() + 
  xlab("Cross-validated LL") +
  ylab("Cross-validated AUROC") +
  scale_color_brewer(palette = "Reds", direction = -1) -> p_ll_auc_n_features_xerr
p_ll_auc_n_features_xerr
# # ggsave("./individual_figures/ob/p_ll_auc_n_features_err.pdf", plot=p_ll_auc_n_features_xerr)

ggplot(df, aes(x=mean_test_w_train_mle_perm,
                                       y=cv_auc_perm, 
                                       ymin=cv_auc_perm-cv_auc_perm_err,
                                       ymax=cv_auc_perm+cv_auc_perm_err, 
                                       color=as.factor(n_features))) +
  geom_point() +
  geom_errorbar() +
  # geom_errorbarh() +
  # geom_pointrangeh(position = position_jitter(w = 0, h = 0), size=0.3) +
  theme_classic() + 
  xlab("Cross-validated LL") +
  ylab("Cross-validated AUROC") +
  scale_color_brewer(palette = "Reds", direction = -1) -> p_ll_auc_n_features
p_ll_auc_n_features
# # ggsave("./individual_figures/ob/p_ll_auc_n_features.pdf", plot=p_ll_auc_n_features)


# geom_label_repel(data = subset(df %>% filter(n_features>0) %>% filter(job_id==412)), aes(label="Model"),
# nudge_x = -2,
# nudge_y = 200)
# geom_point(data = tmp %>% filter(job_id==412), aes(x=n_features, y=mean_test_w_train_mle_perm), col="black", size=3) -> p




# library(ggrepel)
jitter = rnorm(nrow(df), 1, 0.02)
tmp = df %>% mutate(x_jit = n_features*jitter)
tmp2 <- tmp %>% filter(n_features<=2)

ggplot(tmp %>% filter(n_features>2),
       aes(y=cv_auc_perm,
           ymin=cv_auc_perm-cv_auc_perm_err,
           ymax=cv_auc_perm+cv_auc_perm_err,
           x=x_jit)) +
  geom_pointrange(color="black", size=0.3) +    
  geom_pointrange(data = tmp %>% filter(n_features<=2) %>% 
               arrange(match(resnames, single_feature_models)),
             color=brewer.pal(9, "Paired"), size=0.3) +

  theme_classic() + xlab("Number of biological features") +
  ylab("Cross-validated AUROC") +
  scale_x_continuous(breaks = 1:10) -> p_n_features_auc_with_colors
p_n_features_auc_with_colors
# # ggsave("./individual_figures/ob/pngs/p_n_features_auc_with_colors.png", plot=p_n_features_auc_with_colors)
# # ggsave("./individual_figures/ob/pdfs/p_n_features_auc_with_colors.pdf", plot=p_n_features_auc_with_colors) 

ggplot(tmp %>% filter(n_features>2),
       aes(y =mean_test_w_train_mle_perm,
           ymin=mean_test_w_train_mle_perm-se_test_w_train_mle_perm,
           ymax=mean_test_w_train_mle_perm+se_test_w_train_mle_perm,
           x=x_jit)) +
  geom_pointrange(color="black", size=0.3) +
  # geom_pointrange(data = tmp %>% filter(offset_n_features<=2) %>% filter(colors!="black") %>% 
  #              arrange(factor(resnames, levels = factor(single_feature_models))),
  #            color=brewer.pal(9, "Paired")) +
  geom_pointrange(data = tmp2 %>% filter(n_features<=2) %>% 
                    arrange(match(resnames, single_feature_models)),
                  color=brewer.pal(9, "Paired"), size=0.3) +
  
  # geom_errorbar() +
  theme_classic() + xlab("Number of biological features") +
  ylab("Cross-validated LL") +
  scale_x_continuous(breaks = 1:10) + theme(legend.position='none') -> p_n_features_ll_colored
p_n_features_ll_colored

q = .90
top_ll <- quantile(df$lls,  q)
top_auc <- quantile(df$auc_from_er,  q)
best_m_names <- df %>% filter(lls>top_ll) %>% 
  filter(auc_from_er>top_auc) %>% 
  arrange(n_features, n_params) %>% 
  select(n_features, n_params, name, job_id) #%>% View

### best model
source("./ob/load_ob_data.R")

best_m <- res$`160_edges+mutual+nmix(~cell_type)+edgecov(ecov)+nocov(~len)+nicov(~len)+nmatch(~GlomID)+nmatch(~microcluster,levels=c(2:134)).rds`
samples <- readRDS("./ob/results/24_3_new_microglom/samples_160_edges+mutual+nmix(~cell_type)+edgecov(ecov)+nocov(~len)+nicov(~len)+nmatch(~GlomID)+nmatch(~microcluster,levels=c(2:134)).rds")
full_model <- res$`256_edges+mutual+nmix(~cell_type)+edgecov(ecov)+nocov(~len)+nicov(~len)+nocov(~x)+nocov(~y)+nocov(~z)+nicov(~x)+nicov(~y)+nicov(~z)+nmatch(~GlomID)+nmatch(~microcluster,levels=c(2:134)).rds`

er_model <- res$`1_edges.rds`
er_samples <- readRDS("./ob/results/24_3_new_microglom/samples_1_edges.rds")

best_m$activity[2,] %>% mean
best_m$activity[2,] %>% sd

animal = "fish"
animal_path = "ob"
cmap = "darkred"
unique(connectome %>% as.vector)
ncuts = 3
degree_dist_upper_bound = 100
weight_breaks = c(-0.1, 0.1, 1:5*5)
Ng = network.size(g)
# source("../utils/utils.R")

df %>% arrange(activity_mean) %>% tail(1) %>% select(job_id, activity_mean)
worst_decor_model <- res$`57_edges+nocov(~len)+nicov(~len)+nocov(~x)+nocov(~y)+nocov(~z).rds`

df %>% arrange(activity_mean) %>% head(1) %>% select(job_id, activity_mean)
best_decor_model <- res$`168_edges+mutual+nmix(~cell_type)+edgecov(ecov)+nocov(~x)+nocov(~y)+nocov(~z)+nmatch(~GlomID)+nmatch(~microcluster,levels=c(2:134)).rds`

# single feature heatmaps - means
par(mfrow=c(3,3)) 
# single_features = df %>% filter(n_terms<=1) %>% arrange(lls) %>% rownames
# single_features

library(pROC)
ps_means <- list()
for (f in single_feature_models) {
  if (str_detect(f, "mut")) {
    mat <- res[[f]]$prediction_matrices$from_er
  } else {
    mat <- res[[f]]$prediction_matrices$conditional
  }
  ps_means[[f]] <- create_heatmap(mat, "darkred", 
                            title=strsplit(res[[f]]$name, "_")[[1]][2], lims = c(0, 0.25))
}

ps_means_colored <- list()
for (i in 1:length(single_feature_models)) {
  f = single_feature_models[[i]]
  if (str_detect(f, "mut")) {
    mat <- res[[f]]$prediction_matrices$from_er
  } else {
    mat <- res[[f]]$prediction_matrices$conditional
  }
  c = brewer.pal(9, "Paired")[[i]]
  ps_means_colored[[f]] <- create_heatmap(mat, c, 
                                  title=strsplit(res[[f]]$name, "_")[[1]][2], lims = c(0, 0.25), with_guide = "colourbar")
}



ps_single_sample <- list()
for (f in single_feature_models) {
  ps_single_sample[[f]] <- create_heatmap(res[[f]]$closest$from_er*1, "darkred", 
                            title=strsplit(res[[f]]$name, "_")[[1]][2])
}

ps_single_sample_colored <- list()
for (i in 1:length(single_feature_models)) {
  f = single_feature_models[[i]]
  c = brewer.pal(9, "Paired")[[i]]
  ps_single_sample_colored[[f]] <- create_heatmap(res[[f]]$closest$from_er*1, c, 
                                                  title=strsplit(res[[f]]$name, "_")[[1]][2], with_guide = "colourbar")
}


########### distance figure
coords <- cbind(res$`1_edges.rds`$g%v%'x', res$`1_edges.rds`$g%v%'y', res$`1_edges.rds`$g%v%'z')
dm <- rdist::pdist(coords)

distance_dep <- function(g, dm, interval) {
  dm_binned <- cut_width(dm, interval)
  ds <- c()
  ps <- c()
  for (l in levels(dm_binned)[-1]) {
    ds <- c(ds, mean(dm[dm_binned==l]))
    ps <- c(ps, mean(g[dm_binned==l]))
  }
  return(list(ps=ps, ds=ds))
}

g_mat <- res$`1_edges.rds`$g %>% as.matrix
ground_truth = distance_dep(g_mat, dm, 5000)
dm_binned <- cut_width(dm, 5000)
ds_ps_per_model = list(data.frame(distances = ground_truth$ds, type="Data", probabilities = ground_truth$ps))
roc_single_features <- list()
for (i in 1:length(single_feature_models)) {
  p <- single_feature_models[[i]]
  type = p
  tmp = distance_dep(res[[p]]$prediction_matrices$conditional, dm, 5000)
  ds_ps_per_model[[p]] <- data.frame(distances = tmp$ds, type=type, probabilities = tmp$ps)
  roc_single_features[[i]] <- roc(non_diag(res[[p]]$g %>% as.matrix), 
                       non_diag(res[[p]]$prediction_matrices$from_er %>% as.matrix))
  
} 

pal <- brewer.pal(length(roc_single_features), "Paired")
pal[[length(roc_single_features)+1]] <- "black"

names(roc_single_features) <- single_feature_models
roc_single_features[["full"]] <- roc(non_diag(full_model$g %>% as.matrix), 
                                     non_diag(full_model$prediction_matrices$from_er %>% as.matrix))
p_rocs_single_features <- ggroc(roc_single_features) + theme_classic()  + coord_fixed() + scale_color_manual(values = pal)
p_rocs_single_features

p_rocs_single_features_no_legend <- p_rocs_single_features + theme(legend.position = 'none')
p_rocs_single_features_no_legend

mses <- matrix(0, length(single_feature_models), length(single_feature_models))
corrs <- matrix(0, length(single_feature_models), length(single_feature_models))
for (i in (1:length(single_feature_models))) {
  if (single_feature_models[[i]]!="2_edges+mutual.rds") {
    mi = res[[single_feature_models[[i]]]]$prediction_matrices$conditional  
  } else {
    mi = res[[single_feature_models[[i]]]]$prediction_matrices$from_er  
  }
  for (j in (i:length(single_feature_models))) {
    if (single_feature_models[[j]]!="2_edges+mutual.rds") {
      mj = res[[single_feature_models[[j]]]]$prediction_matrices$conditional  
    } else {
      mj = res[[single_feature_models[[j]]]]$prediction_matrices$from_er  
    }
    mses[i,j] = mean((mi-mj)^2)
    mses[j,i] = mses[i,j]
    corrs[i,j] = cor(non_diag(mi), non_diag(mj))
    corrs[j,i] = corrs[i,j]
  }
}

melted_mat <- melt(mses)
s_mses = mses %>% as.matrix %>% dim %>% .[1]
p_mses <- ggplot(melted_mat , aes(x = Var1, y = Var2, fill=value)) +
  geom_raster(interpolate=FALSE) +
  scale_fill_gradient(low="white", high="black") +
  theme_classic() +
  # theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2)) +
        # axis.text.x=element_blank(),
        # axis.title.x=element_blank(),
        # axis.ticks.x=element_blank(),
        # axis.text.y=element_blank(),
        # axis.title.y=element_blank(),
        # axis.ticks.y=element_blank(),
        # plot.title = element_text(hjust = 0.5),
        # panel.background = element_blank()
  # ) + 
  labs(x = "", y = "") +
  # ggtitle(title) +
  coord_fixed() + 
  geom_rect(aes(xmin = .5, xmax = s_mses+.5, ymin = .5, ymax = s_mses+.5),
            fill = "transparent", color = "black", size = 0.1)
p_mses

# # ggsave("./individual_figures/ob/pngs/p_mses.png", plot=p_mses)
# # ggsave("./individual_figures/ob/pdfs/p_mses.pdf", plot=p_mses)


small_cors <- corrs[3:9,3:9]
melted_mat <- melt(small_cors)
s_corrs <- small_cors %>% as.matrix %>% dim %>% .[1]
p_corrs <- ggplot(melted_mat , aes(x = Var1, y = Var2, fill=value)) +
  geom_raster(interpolate=FALSE) +
  scale_fill_gradient(low="white", high="black") +
  theme_classic() +
  # theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2)) +
  # axis.text.x=element_blank(),
  # axis.title.x=element_blank(),
  # axis.ticks.x=element_blank(),
  # axis.text.y=element_blank(),
  # axis.title.y=element_blank(),
  # axis.ticks.y=element_blank(),
  # plot.title = element_text(hjust = 0.5),
  # panel.background = element_blank()
  # ) + 
  labs(x = "", y = "") +
  # ggtitle(title) +
  coord_fixed() + 
  geom_rect(aes(xmin = .5, xmax = s_corrs+.5, ymin = .5, ymax = s_corrs+.5),
            fill = "transparent", color = "black", size = 0.1)

tmp <- do.call(rbind, ds_ps_per_model) 

ggplot(tmp %>% filter(type!="Data"), aes(x=distances, y=probabilities, color=type)) + 
  geom_line() +
  scale_color_brewer(palette="Paired") +
  geom_line(data= tmp %>% filter(type=="Data"), aes(x=distances, y=probabilities), color="black", size=1) +
  theme_classic()  -> p_single_feature_distances
p_single_feature_distances

p_single_features_distances_no_legend <- p_single_feature_distances + theme(legend.position='none')
p_single_features_distances_no_legend


save.image(file = "./ob_workspace_with_res.RData")

res = NULL
conf_samples = NULL
# save.image(file = "./ob_workspace.RData")

source("./gen_figures_for_best_model.R")



########## functional

ex_fig_5 <- function(W_orig, n_mc) {
  all_rows
  if (any(class(W_orig)=="edgelist")) {
    W_orig <- as.network(W_orig) %>% as.matrix
  } else {
    W_orig <- as.matrix(W_orig)
  }
  
  backdivfactor=5;
  ifactor=3;
  INthresh=2.8;
  
  W_in_to_mc <- W_orig[(1+n_mc):dim(W_orig)[1], 1:n_mc]
  W_mc_to_in <- W_orig[1:n_mc, (1+n_mc):dim(W_orig)[1]]
  
  # Normalize them
  for (c in 1:dim(W_in_to_mc)[2]) {
    W_in_to_mc[,c] <- W_in_to_mc[,c]/mean(W_in_to_mc[,c])
  }
  W_in_to_mc[!is.finite(W_in_to_mc)] <- 0.01
  
  for (c in 1:dim(W_mc_to_in)[2]) {
    W_mc_to_in[,c] <- W_mc_to_in[,c]/mean(W_mc_to_in[,c])
  }
  W_mc_to_in[!is.finite(W_mc_to_in)] <- 0.01
  
  # Transpose them
  WINpost=t(W_mc_to_in)
  WMCpost=t(W_in_to_mc)
  
  # Odor indices
  BAodors <- 1:4
  AAodors <- 5:8
  
  # Extract activity in T1 - average over the 5 time steps
  DMref <- apply(Q$DMotMC[,Q$twinstretch1,],c(1,3),mean)
  DMref <- DMref/mean(DMref)
  
  # Check mean correlation between BA stimuli responses
  mean(cor(DMref[,BAodors])[upper.tri(cor(DMref[,BAodors]))])
  
  # Predict IN activations
  INact= WINpost %*% DMref
  
  # Normalize
  help=INact[,BAodors]
  INact[,BAodors] = INact[,BAodors]/sd(help)
  help=INact[,AAodors]
  INact[,AAodors] = INact[,AAodors]/sd(help)
  
  # Threshold
  mask=(INact>=INthresh)
  INact=INact*mask
  
  # Compute and normalize feedback
  MCback= WMCpost %*% INact;
  backmat=MCback;
  backmat=backmat/mean(backmat)
  
  # Subtractive inhibition
  DMnew=DMref-backmat*ifactor;
  mask=DMnew>-0;
  DMnew=DMnew*mask
  
  # Divisive inhibition
  backmatDiv=backmat/max(backmat)
  backmatDiv=backmatDiv*backdivfactor;
  backmatDiv=matrix(1, nrow(backmatDiv), ncol(backmatDiv))-backmatDiv;
  backmatDiv=backmatDiv*(backmatDiv>=0);
  DMnewDiv=DMref*backmatDiv;
  mask=DMnewDiv>-0;
  DMnewDiv=DMnewDiv*mask
  
  # Check correlation in both approaches
  list(vec=c(mean(cor(DMref[,BAodors])[upper.tri(cor(DMref[,BAodors]))]),
             mean(cor(DMnew[,BAodors])[upper.tri(cor(DMnew[,BAodors]))]),
             mean(cor(DMnewDiv[,BAodors])[upper.tri(cor(DMnewDiv[,BAodors]))])),
       new = cor(DMnew), newDiv = cor(DMnewDiv))
  
}




# ggplot(df, aes(x=activity_mean)) + #geom_density(alpha=0.3, show.legend = F) +
#   geom_histogram(data=subset(df,is_mut==T & is_glom==T),fill = "red", alpha = 0.5) +
#   geom_histogram(data=subset(df,is_mut==T & is_glom==F),fill = "blue", alpha = 0.5)
# 
# fit <- lm(df$activity_mean ~ df %>% select(starts_with("is_")) %>% as.matrix)



# # ggsave("./individual_figures/ob/ob_function_all.png", plot=p_func_models)
# # ggsave("./individual_figures/ob/ob_function_all.pdf", plot=p_func_models)

# ggplot(df, aes(x=activity_mean, fill=as.factor(is_mut))) + #geom_histogram(alpha=0.3, show.legend = F) +
#   geom_histogram(data=subset(df,is_mut == T),fill = "red", alpha = 0.5) +
#   geom_histogram(data=subset(df,is_mut == F),fill = "blue", alpha = 0.5) +
#   theme_classic() + xlim(0,.5) + xlab("Mean decorrelation at t2") + 
#   scale_y_continuous(expand = c(0,0)) -> p_func_mut
# p_func_mut
# # ggsave("./individual_figures/ob/ob_function_mutual.png", plot = p_func_mut)
# # ggsave("./individual_figures/ob/ob_function_mutual.pdf", plot = p_func_mut)

  # # ggsave("./individual_figures/ob/ob_function_mutual_glom.png", plot=p_func_both)
# # ggsave("./individual_figures/ob/ob_function_mutual_glom.pdf", plot=p_func_both)

# p_func <- grid.arrange(p_func_models + xlab(""), p_func_mut+ xlab(""), p_func_both)
# # ggsave("./individual_figures/ob/ob_function.png", plot=p_func)
# # ggsave("./individual_figures/ob/ob_function.pdf", plot=p_func)

# p_func_worst + scale_linetype_manual(values = c(1,1,2)) + geom_vline(xintercept = tr_av)
tr1 = Q$DMotMCtr2[,Q$twinstretch2 ,] %>% apply(c(1,3), mean) %>% cor %>% .[1:4,1:4] %>% upper_tri %>% mean()
tr2 = Q$DMotMCtr1[,Q$twinstretch2 ,] %>% apply(c(1,3), mean) %>% cor %>% .[1:4,1:4] %>% upper_tri %>% mean()
tr_av = Q$DMotMC[,Q$twinstretch2 ,] %>% apply(c(1,3), mean) %>% cor %>% .[1:4,1:4] %>% upper_tri %>% mean()

ggplot(data.frame(x=best_m$activity[2,]), aes(x=x)) + geom_histogram(alpha=1, show.legend = F, fill="#8F2727") +
  theme_classic() + xlim(0,.55) + ggtitle("Model") + xlab("") +
  geom_vline(xintercept = tr_av, colour = "black") -> p_func_samples
ggplot(data.frame(x=er_model$activity[2,]), aes(x=x)) + geom_histogram(alpha=1, show.legend = F, fill="#8F2727") +
  theme_classic() + xlim(0,.55)  + ggtitle("ER") + xlab("") +
  geom_vline(xintercept = tr_av, colour = "black") -> p_func_er
ggplot(data.frame(x=worst_decor_model$activity[2,]), aes(x=x)) + geom_histogram(alpha=1, show.legend = F, fill="#8F2727") +
  theme_classic() + xlim(0,.55)  + ggtitle("Smallest decorrelation") + xlab("") +
  geom_vline(xintercept = tr_av, colour = "black")-> p_func_worst
ggplot(data.frame(x=best_decor_model$activity[2,]), aes(x=x)) + geom_histogram(alpha=1, show.legend = F, fill="#8F2727") +
  theme_classic() + xlim(0,.55) + xlab("Stimulus orrelation at t2")  + ggtitle("Highest decorrelation") +
  geom_vline(xintercept = tr_av, colour = "black")-> p_func_best
p_func_samples_combined <- grid.arrange(p_func_worst, p_func_er, p_func_samples, p_func_best, nrow=4)

# # ggsave("./individual_figures/ob/ob_function_samples.png", plot=p_func_samples_combined)
# # ggsave("./individual_figures/ob/ob_function_samples.pdf", plot=p_func_samples_combined)

print("should get here")

#functional - mean
n_stim = 4
upper_bound = 0.7
get_melted_cormat <- function(m, n_stim) {
  diag(m) <- 0
  rownames(m) <- colnames(m) <- (Q$odorlist2 %>% unlist)
  melted_cormat <- melt(m[1:n_stim, 1:n_stim]) %>% mutate(max_val = if_else(value>0.7, 0.7, value))
  return(melted_cormat)
}
worst_decor = df %>% arrange(activity_mean) %>% tail(1) %>% select(activity_mean)
best_decor = df %>% arrange(activity_mean) %>% head(1) %>% select(activity_mean)

model_mean_cormat <- Reduce('+', lapply(best_m$activity_all, function(x) x$new))/length(best_m$activity_all)
worst_mean_cormat <- Reduce('+', lapply(worst_decor_model$activity_all, function(x) x$new))/length(best_m$activity_all)
best_mean_cormat <- Reduce('+', lapply(best_decor_model$activity_all, function(x) x$new))/length(best_m$activity_all)
er_mean_cormat <- Reduce('+', lapply(er_model$activity_all, function(x) x$new))/length(best_m$activity_all)
G_cormat <- ex_fig_5(g %>% as.matrix, n_mc)$new

ub = G_cormat %>% upper_tri %>% max
lb = G_cormat %>% upper_tri %>% min

worst_mean_cormat[worst_mean_cormat>ub] = ub
model_mean_cormat[model_mean_cormat>ub] = ub
best_mean_cormat[best_mean_cormat>ub] = ub

worst_mean_cormat[worst_mean_cormat<lb] = lb
model_mean_cormat[model_mean_cormat<lb] = lb
best_mean_cormat[best_mean_cormat<lb] = lb


#model
ggplot(data = get_melted_cormat(G_cormat, n_stim), aes(x=Var1, y=Var2, fill=max_val)) + 
  geom_tile() + scale_fill_continuous(limits=c(lb, ub)) + #scale_fill_viridis_c(limits = cmap_lim) +
  theme_classic() +
  xlab("Stimulus 1") +
  ylab("Stimulus 2") +
  ggtitle("Data") + 
  coord_equal() -> p_decor_base
p_decor_base

model_cor_rmse <- sqrt(mean((non_diag(G_cormat[1:4, 1:4])-non_diag(model_mean_cormat[1:4, 1:4]))^2))
best_cor_rmse <- sqrt(mean((non_diag(G_cormat[1:4, 1:4])-non_diag(best_mean_cormat[1:4, 1:4]))^2))
worst_cor_rmse <- sqrt(mean((non_diag(G_cormat[1:4, 1:4])-non_diag(worst_mean_cormat[1:4, 1:4]))^2))


(p_decor_data <- p_decor_base)
(p_decor_model <- p_decor_base %+% get_melted_cormat(model_mean_cormat, n_stim) + 
  ggtitle(paste0("Model (RMSE=", round(model_cor_rmse,2), ")")))
(p_decor_best <- p_decor_base %+% get_melted_cormat(best_mean_cormat, n_stim) +
    ggtitle(paste0("Best (RMSE=", round(best_cor_rmse,2), ")")))
(p_decor_worst <- p_decor_base %+% get_melted_cormat(worst_mean_cormat, n_stim) +
    ggtitle(paste0("Worst (RMSE=", round(worst_cor_rmse,2), ")")))


# functional plots - best model
ggplot(df, aes(x=activity_mean)) + geom_histogram(show.legend = F, fill="grey") +
  theme_classic() + xlim(0,.5) + xlab("Mean correlation at t2") + 
  scale_y_continuous(expand = c(0,0)) + 
  geom_vline(xintercept = tr_av, colour = "black") -> p_func_models
p_func_models

tmp_df = df %>% mutate(class = factor(is_mut + 2*is_glom))
ggplot(tmp_df, aes(x=activity_mean)) + #geom_density(alpha=0.3, show.legend = F) +
  geom_histogram(data=subset(tmp_df, class==0), aes(fill=class), alpha = 0.5) +
  geom_histogram(data=subset(tmp_df, class==1), aes(fill=class),  alpha = 0.5) +
  geom_histogram(data=subset(tmp_df, class==2), aes(fill=class), alpha = 0.5) +
  geom_histogram(data=subset(tmp_df, class==3), aes(fill=class), alpha = 0.5) +
  theme_classic() + xlim(0,.5) + xlab("Mean correlation at t2") + 
  geom_vline(xintercept = tr_av, colour = "black") +
  scale_fill_manual(values=c("darkorange","blue","green","red"), labels=c("1","2","3","4")) +
  scale_y_continuous(expand = c(0,0)) -> p_func_both
p_func_both

save.image(file = "./ob_figures.RData")
stop()
