rm(list=ls())
library(pryr)
setwd("./ob")
source("../utils/utils.R")

cmap <- "darkred"
results_dir <- "./ob/results/"

# get job_id
job_id <- as.integer(Sys.getenv(c("LSB_JOBINDEX")))

if (job_id==1) {
  file.copy("./ob/indep_models_CV.R", file.path(results_dir, "indep_models_CV.R"))
}

# get data
source("./ob/load_ob_data.R")
orig_ecov <- ecov
Ng = network.size(g)
N_cv = 10
N_sim = 500

# create formulas
suffixes <- list(mut = c("", "+ mutual"),
                 cell_type = c("","+ nodemix(~cell_type)"),
                 distance = c("", "+ edgecov(ecov)"),
                 olen = c("", " + nodeocov(~len)"),
                 ilen = c("", " + nodeicov(~len)"),
                 o = c("", " + nodeocov(~x) + nodeocov(~y) + nodeocov(~z)"),
                 i = c("", " + nodeicov(~x) + nodeicov(~y) + nodeicov(~z)"),
                 glom = c("", "+ nodematch(~GlomID) + nodematch(~microcluster, levels=c(2:134))"))

formulas <- gen_formulas(suffixes)
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

##############
set.seed(100)
res_CV_perm = list()
for (rep in 1:N_cv) {
  print(rep)
  perm = sample(1:Ng)
  train = perm[1:(Ng/2)]
  test = perm[(Ng/2+1):Ng]
  
  # build g_train
  ecov <- orig_ecov[train, train]
  g_train <- network(binarized_connections[train, train], directed = TRUE, vertex.attr = cellinfo[train,])
  m <- ergm(as.formula(paste0("g_train ", substr(form, 3, nchar(form)))))
  
  ecov <- orig_ecov[test, test]
  # build g_test
  g_test <- network(binarized_connections[test, test], directed = TRUE, vertex.attr = cellinfo[test,])
  m2 <- ergm(as.formula(paste0("g_test ", substr(form, 3, nchar(form)))), , estimate="MPLE")
  
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
ecov <- orig_ecov
model <- ergm(as.formula(form), verbose = T,
              control = snctrl(MCMC.prop = ~strat(attr = ~cell_type, empirical = TRUE) + sparse, seed = 123))

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

gof_from_er <- createGOF(lapply(samples_from_er, as.network) %>% network.list, list(g))
gofs <- list(from_er=gof_from_er)

slim_model <- model[which(unlist(sapply(model, function(x) ifelse(object_size(x)>1e6,F,T))))]

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

activity_all <- lapply(lapply(samples_from_er, function(nw) {as.matrix(as.network(nw))}), function(s) ex_fig_5(s, n_mc))
activity <- sapply(activity_all, function(l) l$vec)

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
                gofs = gofs,
                g = g,
                form=form, 
                motifs = motifs,
                idegs = idegs,
                odegs = odegs,
                indep = is.dyad.independent(model),
                distances = distances,
                # CV = tmp2,
                CV_perm = tmp2_perm,
                activity = activity,
                activity_all = activity_all
                )

# save results                
library(Matrix)
saveRDS(results, file = file.path(results_dir, paste0("res_", name, ".rds")))
saveRDS(samples_from_er, file = file.path(results_dir, paste0("samples_", name, ".rds")))
# saveRDS(lapply(samples_from_er, function(m) as(m, "sparseMatrix")), file = file.path(results_dir, paste0("sparse_samples_", name, ".rds")))
saveRDS(model, file = file.path(results_dir, paste0("model_", name, ".rds")))

