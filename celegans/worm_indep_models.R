rm(list=ls())
library(pryr)
setwd("./celegans/")
source("../utils/utils.R")

cmap <- "darkblue"
results_dir <- "./celegans/results/15_3_diff_inference/"

# get job_id
job_id <- as.integer(Sys.getenv(c("LSB_JOBINDEX")))

# an attempt to avoid the problem of all the scripts accessing the same file at the same time
Sys.sleep(runif(1, min=0, max=10))

if (job_id==1) {
  file.copy("./celegans/worm_indep_models.R", file.path(results_dir, "worm_indep_models.R"))
}


# get data
source("./celegans/load_worm_data.R")

Ng = network.size(g)
N_cv = 10
N_sim = 500

suffixes <- list(type = c("", "+ nodemix(~subtype)"),
                 ganglion = c("", "+ nodemix(~ganglion)"),
                 region = c("", "+ nodemix(~region)"),
                 rich_club = c("", "+ nodemix(~rich_club)"),
                 time = c("", "+ nodeicov(~poly(time,2)) + nodeocov(~poly(time,2))"),
                 pos = c("", "+ nodeicov(~poly(pos,2)) + nodeocov(~poly(pos,2))"),
                 distance = c("", "+ absdiff(~pos)"),
                 reciprocity = c("", "+ mutual"))
                 # first_char = c("", "+ nodemix(~first_char_low_res)"),
                 # last_char = c("", "+ nodemix(~last_char_low_res)"),
                 # trans = c("", "+ nodemix(~trans)")
                 # )
formulas <- gen_formulas(suffixes)

control_list_version = (job_id-1)%/%length(formulas)+1
job_id = (job_id-1)%%length(formulas)+1

# formulas[[+1]] <- "g ~ sender + receiver"

term_counts <- do.call(expand.grid, lapply(suffixes, function(l) sapply(l, function(x) as.double(nchar(x)>0))))

names <- get_names(formulas)

name <- names[[job_id]]
name <- paste0("algo_", control_list_version, "_", name)
print("name")
print(name)

form <- formulas[[job_id]]
if (job_id <= (length(formulas)-2)) {
  n_terms <- rowSums(as.matrix(term_counts))[[job_id]]
} else if (job_id == (length(formulas)-1)) {
  n_terms <- 2
} else {
  n_terms <- 3
}

if (control_list_version==1) {
  control_list = snctrl(MCMC.prop = ~strat(attr = ~subtype, empirical = TRUE) + sparse,
                        seed = 123, 
                        main.method = "Stochastic-Approximation")
} else if (control_list_version==2) {
  control_list = snctrl(MCMC.prop = ~strat(attr = ~subtype, empirical = TRUE) + sparse,
                        MCMLE.termination = "Hotelling",
                        seed = 123)
} else if (control_list_version==3) {
  control_list = snctrl(MCMC.prop = ~strat(attr = ~subtype, empirical = TRUE) + sparse,
                        MCMLE.termination = "confidence",
                        seed = 123)
} else if (control_list_version==4) {
  control_list = snctrl(MCMC.prop = ~strat(attr = ~subtype, empirical = TRUE) + sparse,
                        MCMLE.termination = "Hotelling",
                        seed = 123, 
                        MCMLE.effectiveSize = NULL,
                        MCMC.burnin = 5e5,
                        MCMC.interval = 5e4,
                        MCMC.samplesize = 7500)
                        
} else if (control_list_version==5) {
  control_list = snctrl(MCMC.prop = ~strat(attr = ~subtype, empirical = TRUE) + sparse,
                        MCMLE.termination = "confidence",
                        seed = 123, 
                        MCMLE.effectiveSize = NULL,
                        MCMC.burnin = 5e5,
                        MCMC.interval = 5e4,
                        MCMC.samplesize = 7500)
}


#####################
# same for random perm
#####################

# train_test <- function(rep) {
#   print(rep)
#   train = c()
#   test = c()
#   
#   # stratified by subtype
#   for (subtype in (cellinfo$subtype %>% unique)) {
#     idx = which(cellinfo$subtype==subtype)
#     idx = sample(idx)
#     l = length(idx)
#     if (l%%2 == 0) {
#       train = c(train, idx[1:(l%/%2)])
#       test = c(test, idx[((l%/%2)+1):l])
#     } else {
#       if (length(train)>length(test)) {
#         train = c(train, idx[1:(l%/%2)])
#         test = c(test, idx[((l%/%2)+1):l])
#       } else {
#         train = c(train, idx[1:(l%/%2+1)])
#         test = c(test, idx[(l%/%2+2):l])
#       }
#     }
#   }
#   
#   
#   # build g_train
#   g_train <- network(bin_connectome[train,train], vertex.attr = cellinfo[train,])
#   # build g_test
#   g_test <- network(bin_connectome[test,test], vertex.attr = cellinfo[test,], estimate = "MPLE")
#   
#   m <- ergm(as.formula(paste0("g_train ", substr(form, 3, nchar(form)))), control=control_list)
#   m2 <- ergm(as.formula(paste0("g_test ", substr(form, 3, nchar(form)))), control=control_list)
#   prev_ll <- logLik(m2)
#   
#   m2$estimate <- "MLE"
#   m2$MPLE_is_MLE <- F
#   m2$coefficients <- coef(m)
#   
#   x_perm <- simulate(m2, nsim = N_sim, control = control.simulate.ergm(MCMC.effectiveSize = NULL, MCMC.burnin = 1e6, MCMC.interval = 1e5), 
#                      verbose = T)
#   samples_from_er_perm <- x_perm[!unlist(lapply(x_perm, is.null))]
#   mean_samples_from_er_perm <- mean_connectivity(samples_from_er_perm)
#   test_auc = get_metrics(mean_samples_from_er_perm, g_test, "auc")
#   
#   return(list(ll_train = logLik(m), ll_test = prev_ll, ll_test_new = logLik(m2, force.reeval = T), test_auc = test_auc))
# }
# 
# res_CV_perm = lapply(1:N_cv, function(rep) tryCatch(train_test(rep), 
#                                                   error = function(e) list(ll_train = NA, ll_test = NA, ll_test_new = NA, test_auc = NA)))
# 
# tmp_perm <- do.call(rbind, res_CV_perm)
# tmp2_perm = matrix(tmp_perm %>% as.double, ncol=4)

# fit model
model <- ergm(as.formula(form), verbose = T, control=control_list)

if (is.dyad.independent(model)) {
  samples_from_er <- list()
  conditional <- predict(model, conditional=T, output="matrix")
  mean_samples_from_er <- matrix(0, Ng, Ng)
  for (i in 1:N_sim) {
    current <- matrix(runif(Ng^2), Ng, Ng)
    samples_from_er[[i]] <- as.matrix(current<conditional, Ng, Ng)
    mean_samples_from_er = mean_samples_from_er + as.integer(current<conditional)
  }
  mean_samples_from_er <- mean_samples_from_er/N_sim
} else {
  print("dep model")
  x <- simulate(model, nsim = N_sim, control = control.simulate.ergm(MCMC.effectiveSize = NULL, MCMC.interval = 1e6), 
                verbose = T)
  samples_from_er <- x[!unlist(lapply(x, is.null))]
  mean_samples_from_er <- mean_connectivity(samples_from_er)
}

# hamming distance of samples
hamming_from_er <- hdist(samples_from_er, g)
hammings <- list(from_er=hamming_from_er)

# closest single samples
closest_sample_from_er <- samples_from_er[[which.min(hamming_from_er)]]
closest <- list(from_er=as.matrix(closest_sample_from_er))

conditional <- predict(model, conditional=T, output="matrix")
prediction_matrices <- list(from_er=mean_samples_from_er, conditional=conditional)

aucs <- lapply(prediction_matrices, function(m) get_metrics(m, g, "auc"))
prs <- lapply(prediction_matrices, function(m) get_metrics(m, g, "aucpr"))

model_summary <- summary(model)

motifs <- sapply(samples_from_er, function(nw) summary(nw ~ triadcensus))
if (samples_from_er[[1]] %>% class == "matrix") {
  distances <- sapply(samples_from_er, function(nw) {tmp <- ergm.geodistdist(nw %>% as.network(directed=T)); return(c(tmp[length(tmp)], tmp[1:6]))})
  
} else {
  distances <- sapply(samples_from_er, function(nw) {tmp <- ergm.geodistdist(nw); return(c(tmp[length(tmp)], tmp[1:6]))})
  
}
idegs <- sapply(samples_from_er, function(nw) sna::degree(as.network(nw), cmode="indegree"))
odegs <- sapply(samples_from_er, function(nw) sna::degree(as.network(nw), cmode="outdegree"))

slim_model <- model[which(unlist(sapply(model, function(x) ifelse(object_size(x)>1e6,F,T))))]

results <- list(model = slim_model,
                hammings = hammings,
                closest = closest,
                prediction_matrices = prediction_matrices,
                aucs = aucs,
                prs = prs,
                summary = model_summary,
                aic = model_summary$aic,
                bic = model_summary$bic,
                ll = logLik(model),
                n_terms = n_terms,
                name = name,
                job_id = job_id, 
                # gofs = gofs,
                g = g,
                form=form, 
                motifs = motifs,
                idegs = idegs,
                odegs = odegs,
                indep = is.dyad.independent(model),
                distances = distances,
                # CV_perm = tmp2_perm,
                control_list_version = control_list_version

)

# save results                
saveRDS(results, file = file.path(results_dir, paste0("res_", name, ".rds")))
saveRDS(samples_from_er, file = file.path(results_dir, paste0("samples_", name, ".rds")))
saveRDS(model, file = file.path(results_dir, paste0("model_", name, ".rds")))

