rm(list=ls())
library(pryr)
setwd("./allen")
source("../utils/utils.R")

cmap <- "Oranges"
results_dir <- "./allen/results/15_11_with_samples/"

# get job_id
job_id <- as.integer(Sys.getenv(c("LSB_JOBINDEX")))

if (job_id==1) {
  file.copy("./allen/allen_indep_CV.R", file.path(results_dir, "allen_indep_CV.R"))
}

# get data
source("./allen/load_allen_data.R")
dm = dm/mean(dm)
orig_dm = dm
Ng = network.size(g)
N_cv = 10
N_sim = 500

suffixes <- list(distance = c("", "+ edgecov(dm)"),
                 reciprocity = c("", "+ mutual", "+ asymmetric"),
                 spatial_linear_i = c("", "+ nodeicov(~x) + nodeicov(~y) + nodeicov(~z)", "+ nodeicov(~poly(x,2)) + nodeicov(~poly(y,2)) + nodeicov(~poly(z,2))"),
                 spatial_linear_o = c("", "+ nodeocov(~x) + nodeocov(~y) + nodeocov(~z)", "+ nodeocov(~poly(x,2)) + nodeocov(~poly(y,2)) + nodeocov(~poly(z,2))"),
                 paths = c("", "+ twopath"))
formulas <- gen_formulas(suffixes)
formulas[[length(formulas)+1]] <- "g ~ sender + receiver"
formulas[[length(formulas)+1]] <- "g ~ sender + receiver + mutual"
term_counts <- do.call(expand.grid, lapply(suffixes, function(l) sapply(l, function(x) as.double(nchar(x)>0))))

names <- get_names(formulas)

name <- names[[job_id]]
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

############ cv part
set.seed(10)
res_CV_perm = list()
for (rep in 1:N_cv) {
  print(rep)
  perm = sample(1:334)
  train = perm[1:167]
  test = perm[168:334]

  # build g_train
  dm <- orig_dm[train, train]
  g_train <- network(W_bin[train, train], directed = TRUE, vertex.attr = loc[train,])
  m <- ergm(as.formula(paste0("g_train ", substr(form, 3, nchar(form)))))
  
  # build g_test
  dm <- orig_dm[test, test]
  g_test <- network(W_bin[test, test], directed = TRUE, vertex.attr = loc[test,])
  m2 <- ergm(as.formula(paste0("g_test ", substr(form, 3, nchar(form)))), estimate="MPLE")
  
  prev_ll <- logLik(m2)
  
  m2$estimate <- "MLE"
  m2$MPLE_is_MLE <- F
  m2$coefficients <- coef(m)
  
  x_perm <- simulate(m2, nsim = N_sim, control = control.simulate.ergm(MCMC.effectiveSize = NULL, MCMC.burnin = 1e6, MCMC.interval = 1e5), 
                     verbose = T)
  samples_from_er_perm <- x_perm[!unlist(lapply(x_perm, is.null))]
  mean_samples_from_er_perm <- mean_connectivity(samples_from_er_perm)
  test_auc = get_metrics(mean_samples_from_er_perm, g_test, "auc")
  
  res_CV_perm[[rep]] = list(ll_train = logLik(m), ll_test = prev_ll, ll_test_new = logLik(m2, force.reeval = T), test_auc = test_auc)
}

##############
tmp_perm <- do.call(rbind, res_CV_perm)
tmp2_perm = matrix(tmp_perm %>% as.double, ncol=4)

# fit model
dm <- orig_dm
model <- ergm(as.formula(form), verbose = T)

if (is.dyad.independent(model)) {
  samples_from_er <- list()
  conditional <- predict(model, conditional=T, output="matrix")
  mean_samples_from_er <- matrix(0, Ng, Ng)
  for (i in 1:N_sim) {
    current <- matrix(runif(Ng^2), Ng, Ng)
    samples_from_er[[i]] <- as.matrix(current<conditional, Ng, Ng)
    mean_samples_from_er = mean_samples_from_er + as.integer(current<conditional)/N_sim
  }
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
                # CV = tmp2,
                CV_perm = tmp2_perm
)

# save results                
saveRDS(results, file = file.path(results_dir, paste0("res_", name, ".rds")))
saveRDS(samples_from_er, file = file.path(results_dir, paste0("samples_", name, ".rds")))
saveRDS(model, file = file.path(results_dir, paste0("model_", name, ".rds")))



