library(ergm)
library(network)
library(sna)
library(ggplot2)
library(GGally)
library(reshape2)
library(statnet)
library(Bergm)
library(gridExtra)
# library(abind)
library(rdist)
library(latentnet)
# library(btergm)
library(readxl)
library(Matrix)
library(stringr)
library(xergm.common)
library(ROCR)
library(Rfast)
library(readr)
library(dplyr)
library(igraph)
library(intergraph)
source("./gofplots.R")
source("./gofs.R")
source("./edgeprob_utils.R")



run_forumula <- function (formula) {
  model <- ergm(formula)
  print(summary(model))
  model
}

get_gof <- function (m) {
  gof(m, GOF = ~distance + espartners + idegree + odegree + model)
}

plot_gof <- function (g,i) {
  pdf(paste('model_', i,'_gof.pdf', sep=''))
  par(mfrow = c(3,3), mar=c(4,4,4,4))
  plot(g, cex.lab=1.6, cex.axis=1.6, plotlogodds = FALSE)
  dev.off()
}

# # terribly slow, not sure why
# gen_basis <- function(g, n) {
#   graphs_array <- rgraph(network.size(g), m=n, tprob=network.density(g))
#   graphs_list <- lapply(1:n, function(x) graphs_array[ x , , ])
#   lapply(graphs_list, function(x) update(network.copy(g), x, matrix.type="adjacency"))
# }

get_samples_from_er <- function(model, g, nsim=100, MCMC.interval = 1e6, MCMC.burnin = 1e6) {
  basis <- update(network.copy(g), rgraph(network.size(g), m=1, tprob=network.density(g)), matrix.type="adjacency")
  samples <- simulate(model, nsim = nsim, basis=basis, control = control.simulate.ergm(MCMC.effectiveSize = 100), verbose = T)
  return(samples)
}

get_samples_from_g <- function(model, nsim=100, MCMC.interval = 1e6, MCMC.burnin = 1e6) {
  samples <- simulate(model, nsim = nsim, control = control.simulate.ergm(MCMC.effectiveSize = 100), verbose = T)
  return(samples)
}

hamming_from_samples <- function(g, samples) {
  return(hdist(samples, g))
}


# sample_model <- function(m, g, n=1000) {
#   basis <- update(network.copy(g), rgraph(network.size(g), m=1, tprob=network.density(g)), matrix.type="adjacency")
#   m_list <- lapply(simulate(m, nsim=n, sequential = T, basis = basis), FUN=as.matrix)
# }

# sample_model_from_basis <- function(m, basis, n=1000, MCMC.burnin=1000, MCMC.interval=100, seq=T) {
#   m_list <- lapply(simulate(m, nsim=n, sequential = seq, basis = basis, control=control.simulate.ergm(
#     MCMC.burnin=MCMC.burnin,
#     MCMC.interval=MCMC.interval)
#     ),
#     FUN=as.matrix)
# }

mean_connectivity <- function(samples) {
  return(mean_sample(lapply(samples, as.matrix)))  
}

mean_connectivity_edgelist <- function(samples, N) {
  l <- length(samples)
  mat <- matrix(0, N, N)
  for (s in samples) {
    mat[as.matrix(s)] = mat[as.matrix(s)] + 1
  }
  return(mat/l)
}

mean_sample <- function(l) {
  Reduce("+", l) / length(l)
}

create_heatmap <- function(mat, cmap, title = "", lims=NULL, bounding_box=TRUE,
                           with_guide=F, transpose=T,
                           xlab = "Postsynaptic", ylab = "Presynaptic") {
  if (transpose==T) {
    mat <- t(mat)
  }
  rownames(mat) <- colnames(mat) <- NULL
  melted_mat <- melt(mat)
  p <- ggplot(melted_mat , aes(x = Var1, y = Var2, fill=value)) +
    geom_raster(interpolate=FALSE) +
    scale_fill_gradient(low="white", high=cmap, guide = with_guide, limits=lims) + 
    theme(legend.position = "right",
          axis.text.x=element_blank(),
          # axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          # axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.background = element_blank()
          ) + 
    labs(x = xlab, y = ylab) +
    ggtitle(title) +
    coord_fixed()
  s = mat %>% as.matrix %>% dim %>% .[1]
  if (bounding_box==T) {
    return(p + geom_rect(aes(xmin = .5, xmax = s+.5, ymin = .5, ymax = s+.5),
                             fill = "transparent", color = "black", size = 0.1))
  } else {
    return(p)
  }
}


create_heatmap_tile <- function(mat, cmap, title = "", lim=1) {
  melted_mat <- melt(mat)
  p <- ggplot(melted_mat , aes(x = Var1, y = Var2, fill=value)) +
    geom_tile(interpolate=FALSE) +
    scale_fill_gradient(low="white", high=cmap, guide = F) + 
    theme(legend.position = "right",
          axis.text.x=element_blank(),
          axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.background = element_blank()) + 
    ggtitle(title) +
    coord_fixed()
  
  p;
}

indep_lp <- function(m_orig, g) {
  m <- m_orig
  m$coef <- ifelse(is.infinite(m$coef),-10,m$coef)
  t <- edgeprob(m, verbose = T)
  tt <- matrix(0,dim(as.matrix(g))[[1]],dim(as.matrix(g))[[2]])
  counter <- 1 
  for(i in 1:dim(as.matrix(g))[[2]]) {
    for (j in 1:dim(as.matrix(g))[[2]]) {
      if (i == j) {
        next
      }
      else {
        tt[i,j] <- t[counter,dim(t)[2]]
        counter <- counter + 1
      }
    }
  }
  tt
}

indep_lp_sym <- function(m_orig, g) {
  m <- m_orig
  m$coef <- ifelse(is.infinite(m$coef),-10,m$coef)
  t <- edgeprob(m, verbose = T)
  n <- dim(as.matrix(g))[[1]]
  tt <- squareform(t[,dim(t)[2]])
  tt + t(tt)
}

non_diag <- function(m) {
  c(upper_tri(m),lower_tri(m))
}

gen_formulas <- function(suffixes) {
  all_models_suffixes <- do.call(expand.grid, suffixes)
  all_models_suffixes <- do.call(paste, all_models_suffixes)
  all_models_suffixes <- lapply(as.list(all_models_suffixes), trimws)
  all_formulas <- lapply(all_models_suffixes, function(x) paste("g ~ edges ", x))
  all_formulas
}

get_names <- function(formulas) {
  lapply(1:length(formulas), function(i) str_replace_all(str_remove_all(paste0(i, "_", substr(formulas[i],4,nchar(formulas[i]))), '"| |g~'),"node","n")) 
}

get_metrics <- function(m, g, metric) {
  pred <- prediction(non_diag(m),non_diag(as.matrix(g)>0))
  perf <- performance(pred, metric)@y.values[[1]]
  return(perf)
}



run_comp_table <- function(wd_path, g, W, formula, suffixes, control_list){
  
  gof <- createGOF(samples, list(g))
  
  marginals <- tryCatch(indep_lp_custom(model, g), error = function(e) {print("returning NULL marginals"); matrix(runif(279^2), ncol=279)})
  
  l <- list(g = g,
            s = summary(model), 
            m = model,
            gof = gof,
            job_id = job_id,
            mean_mat = mean_sample1,
            mean_mat_g = mean_sample_g,
            marginals = marginals,
            hamming_from_g = hamming_from_g,
            hamming_from_init = hamming_from_init,
            min_hamming_g = min_hamming_g,
            min_hamming_init = min_hamming_init,
            samples_from_g = samples_from_g,
            samples_from_er = samples,
            terms = form,
            n_terms = n_terms)
  
  # l2 <- list(g = l$g,
  #            motif_counts = sapply(l$samples_from_g, function(x) motifs(intergraph::asIgraph(x))),
  #            s = summary(l$m), 
  #            m = l$m,
  #            bic = summary(l$m)$bic,
  #            aic = summary(l$m)$aic,
  #            job_id = l$job_id,
  #            mean_mat = l$mean_mat,
  #            mean_mat_g = l$mean_mat_g,
  #            marginals = l$marginals,
  #            pred_marginals = prediction(non_diag(l$marginals),non_diag(as.matrix(l$g)>0)),
  #            pred_mean_mat = prediction(non_diag(l$mean_mat_g),non_diag(as.matrix(l$g)>0)),
  #            # infer_out = mean_sample2,
  #            # infer_in = mean_sample3,
  #            hamming_from_g = l$hamming_from_g,
  #            hamming_from_init = l$hamming_from_init,
  #            # min_hamming_g = l$min_hamming_g,
  #            # min_hamming_init = l$min_hamming_init,
  #            terms = l$terms,
  #            n_terms = n_terms)
  
  # save(l2, file=paste("Table summary",job_id,"with metadata.Rda"))
  
  
  pdf(paste("Model ",name,".pdf"))
  par(mfrow=c(2,2))
  image(as.matrix(g), main="Data")
  image(as.matrix(min_hamming_g), main="Closest (g)")
  image(as.matrix(min_hamming_init), main="Closest (ER)")
  image(marginals, main="Analytic")
  
  par(mfrow=c(2,2))
  # image(mean_sample2, main="Infer out")
  # image(mean_sample3, main="Infer in")
  image(mean_sample1,main="MCMC; initial = random ER")
  image(mean_sample_g, main="MCMC; initial = original graph")
  
  par(mfrow=c(2,2))
  lim = c(min(hamming_from_g, hamming_from_init), max(hamming_from_g, hamming_from_init))
  hist(hamming_from_g, xlim = lim, main="Hamming (g)")
  hist(hamming_from_init, xlim = lim, main="Hamming (ER)")
  hist(marginals[W>0], main="P(syn|syn exists)", xlim = c(0,1))
  hist(marginals[W==0], main="P(syn|syn doesn't exist)", xlim = c(0,1))
  
  par(mfrow=c(1,1))
  
  plot(gof)
  
  dev.off()
  l
}
