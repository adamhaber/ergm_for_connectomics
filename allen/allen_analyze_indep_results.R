setwd("./allen/")
library(stringr)
library(ROCR)
library(dplyr)
library(ggplot2)
RDS_files <- # path to where all the saved runs are at
filenames <- Sys.glob(RDS_files)
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
n_terms <- sapply(res, function(x) x$n_terms)
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
median_test_w_train_ml_perme <- sapply(res, function(r) r$CV_perm[,3] %>% median)

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


cv_auc <- sapply(res, function(r) r$CV_perm[,4] %>% mean)
cv_auc_err <- sapply(res, function(r) r$CV_perm[,4] %>% stderr)

g <- res[[1]]$g
gs_motifs <- summary(g~triadcensus)
gs_ideg <- colSums(as.matrix(g))
gs_odeg <- rowSums(as.matrix(g))
iterations <- sapply(res, function(x) x$model$iterations)


motifs_accuracy <- sapply(res, function(x) rowMeans(x$motifs/gs_motifs)-1) %>% na.omit %>% abs %>% colMeans
ideg_accuracy <- sapply(res, function(x) tryCatch(cor(gs_ideg, rowMeans(x$idegs)), error = function(e) -1)) 
odeg_accuracy <- sapply(res, function(x) tryCatch(cor(gs_odeg, rowMeans(x$odegs)), error = function(e) -1))

df <- data.frame(filenames=current_fnames,
                 lls=lls,
                 iter = iterations,
                 n_params = n_params,
                 n_terms=n_terms,
                 auc_from_cond=auc_from_cond, auc_from_er=auc_from_er,
                 pr_from_cond=pr_from_cond, pr_from_er=pr_from_er,
                 ham_from_er=ham_from_er, indep=indep,
                 name=name, job_id=job_id,
                 had_nans=had_nans,
                 mean_train = mean_train,
                 mean_test = mean_test,
                 mean_test_w_train_mle = mean_test_w_train_mle,
                 mean_train_perm = mean_train_perm,
                 mean_test_perm = mean_test_perm,
                 mean_test_w_train_mle_perm = mean_test_w_train_mle_perm,
                 
                 se_train = se_train,
                 se_test = se_test,
                 se_test_w_train_mle = se_test_w_train_mle,
                 cv_auc = cv_auc,
                 cv_auc_err = cv_auc_err,
                 
                 se_train_perm = se_train_perm,
                 se_test_perm = se_test_perm,
                 se_test_w_train_mle_perm = se_test_w_train_mle_perm,
                 
                 is_mutual=sapply(current_fnames, function(x) str_detect(x, "mutual")),
                 is_dm=sapply(current_fnames, function(x) str_detect(x, "edgecov")),
                 is_twopath=sapply(current_fnames, function(x) str_detect(x, "twopath")),
                 motifs_accuracy=motifs_accuracy,
                 ideg_accuracy=ideg_accuracy, odeg_accuracy=odeg_accuracy
)

df$offset_n_terms <- df$n_terms + 1

df <- df %>% filter(!str_detect(name, "sender"))
######################## n terms
ggplot(df, aes(ymin=mean_test_w_train_mle_perm-se_test_w_train_mle_perm, 
               y=mean_test_w_train_mle_perm,
               ymax=mean_test_w_train_mle_perm+se_test_w_train_mle_perm, 
               x=offset_n_terms)) +
  geom_pointrange(position = position_jitter(w = 0.15, h = 0), size=0.2) +
  theme_classic() + xlab("Number of biological features") +
  ylab("Cross-validated LL") +
  scale_x_continuous(breaks = 0:10) + theme(legend.position='none') -> p_n_features_ll
# geom_label_repel(data = subset(df %>% filter(n_terms>0) %>% filter(job_id==412)), aes(label="Model"),
# nudge_x = -2,
# nudge_y = 200)
# geom_point(data = tmp %>% filter(job_id==412), aes(x=n_terms, y=mean_test_w_train_mle_perm), col="black", size=3) -> p
p_n_features_ll
# ggsave("./individual_figures/allen/p_cv_ll_features.pdf", plot=p_n_features_ll)
# ggsave("./final figures/fig1e.png", plot=p_n_features_ll)
# ggsave("./final figures/fig1e.pdf", plot=p_n_features_ll)

# auc vs ll vs n_terms
tmp <- df
tmp$n_terms[tmp$n_terms==0] <- 1
ggplot(tmp %>% filter(n_terms>-1), aes(x=mean_test_w_train_mle_perm,
                                       xmin=mean_test_w_train_mle_perm-se_test_w_train_mle_perm,
                                       xmax=mean_test_w_train_mle_perm+se_test_w_train_mle_perm,
                                       y=cv_auc, 
                                       ymin=cv_auc-cv_auc_err,
                                       ymax=cv_auc+cv_auc_err, 
                                       color=as.factor(n_terms))) +
  geom_point() +
  geom_errorbar() +
  geom_errorbarh() +
  # geom_pointrangeh(position = position_jitter(w = 0, h = 0), size=0.3) +
  theme_classic() + 
  xlab("Cross-validated LL") +
  ylab("Cross-validated AUROC") +
  scale_color_brewer(palette = "Greens", direction = -1) +
  ylab("Cross-validated LL") +
  theme(legend.position = 'none') -> p_ll_auc_n_terms_err
p_ll_auc_n_terms_err
# ggsave("./individual_figures/allen/pdfs/p_ll_auc_n_terms_err.pdf", plot=p_ll_auc_n_terms)
# ggsave("./individual_figures/allen/pngs/p_ll_auc_n_terms_err.png", plot=p_ll_auc_n_terms)


ggplot(tmp %>% filter(n_terms>-1), aes(x=mean_test_w_train_mle_perm,
                                       y=cv_auc, 
                                       ymin=cv_auc-cv_auc_err,
                                       ymax=cv_auc+cv_auc_err, 
                                       color=as.factor(n_terms))) +
  geom_point() +
  geom_errorbar() +
  # geom_errorbarh() +
  # geom_pointrangeh(position = position_jitter(w = 0, h = 0), size=0.3) +
  theme_classic() + 
  xlab("Cross-validated LL") +
  ylab("Cross-validated AUROC") +
  scale_color_brewer(palette = "Greens", direction = -1) +
  ylab("Cross-validated LL") -> p_ll_auc_n_terms
p_ll_auc_n_terms
# ggsave("./individual_figures/allen/p_ll_auc_n_terms.pdf", plot=p_ll_auc_n_terms)


# geom_label_repel(data = subset(df %>% filter(n_terms>0) %>% filter(job_id==412)), aes(label="Model"),
# nudge_x = -2,
# nudge_y = 200)
# geom_point(data = tmp %>% filter(job_id==412), aes(x=n_terms, y=mean_test_w_train_mle_perm), col="black", size=3) -> p




# library(ggrepel)
ggplot(df,
       aes(ymin=cv_auc-cv_auc_err,
           y=cv_auc,
           ymax=cv_auc+cv_auc_err,
           x=offset_n_terms)) +
  geom_pointrange(position = position_jitter(w = 0.2, h = 0), size=0.2) +
  theme_classic() + xlab("Number of biological features") +
  ylab("Cross-validated AUROC") +
  scale_x_continuous(breaks = 1:10) + theme(legend.position='none') -> p_n_features_auc
# geom_label_repel(data = subset(df %>% filter(n_terms>0) %>% filter(job_id==412)), aes(label="Model"),
#                  nudge_x = -2,
#                  nudge_y = 0.01) +
# geom_point(data = tmp %>% filter(job_id==412), aes(x=n_terms, y=cv_auc_perm), col="black", size=3) -> p
p_n_features_auc
# ggsave("./individual_figures/allen/p_cv_auc.pdf", plot=p_n_features_auc) 

q = .90
top_ll <- quantile(df$lls,  q)
top_auc <- quantile(df$auc_from_er,  q)
best_m_names <- df %>% filter(lls>top_ll) %>% 
  filter(auc_from_er>top_auc) %>% 
  arrange(n_terms, n_params) %>% 
  select(name)
best_m_names %>% unlist %>% first





ggplot(df, aes(x=n_terms, y=cv_auc, color=as.factor(is_twopath))) + geom_point(position = position_jitter(w = 0.2, h = 0)) +
  theme_classic()

best_m <- res$`50_edges+edgecov(dm)+nicov(~poly(x,2))+nicov(~poly(y,2))+nicov(~poly(z,2))+nocov(~poly(x,2))+nocov(~poly(y,2))+nocov(~poly(z,2)).rds`
samples <- readRDS("./allen/results/15_11_with_samples/samples_50_edges+edgecov(dm)+nicov(~poly(x,2))+nicov(~poly(y,2))+nicov(~poly(z,2))+nocov(~poly(x,2))+nocov(~poly(y,2))+nocov(~poly(z,2)).rds")
full_model <- res$`108_edges+edgecov(dm)+asymmetric+nicov(~poly(x,2))+nicov(~poly(y,2))+nicov(~poly(z,2))+nocov(~poly(x,2))+nocov(~poly(y,2))+nocov(~poly(z,2))+twopath.rds`

animal = "allen"
animal_path = "allen"
cmap = "darkgreen"
ncuts = 3

df_orig <- df
source("./allen/load_allen_data.R")
source("../utils/utils.R")
Ng = network.size(g)
degree_dist_upper_bound = 25
W_weights[W_weights>0] %>% as.vector %>% hist(breaks=100)
connectome <- W_counts
weight_breaks = -1:6
er_model <- res$`1_edges.rds`
er_samples <- readRDS("./allen/results/15_11_with_samples/samples_1_edges.rds")

res = NULL
conf_samples = NULL
save.image(file = "./allen_workspace.RData")
library(pROC)
source("./gen_figures_for_best_model.R")
save.image(file = "./allen_figures.RData")
