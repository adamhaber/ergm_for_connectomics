source("../utils/utils.R")
source("./celegans/load_worm_data.R")

library(optimization)
library(dfoptim)
library(pso)

job_idx <- as.integer(Sys.getenv(c("LSB_JOBINDEX")))

cellinfo <- cellinfo %>% mutate(M = if_else(region=="M", "M", "E"))
idx = cellinfo %>% filter(type=="motorneuron") %>% filter(M=="M") %>% select(name) %>% unlist
motor_connectome <- connectome[idx, idx]
bin_connectome <- (motor_connectome > 0)

full_g <- g
g <- network(bin_connectome>0, vertex.attr = cellinfo %>% filter(name %in% idx))

baseline_m <- ergm(g ~ absdiff(~pos) + nodemix(~first_two))
off <- baseline_m$coefficients[[1]]

best_m <- readRDS("./celegans/results/10_3_no_subtype_ext/res_1014_edges+nmix(~subtype)+nifactor(~ganglion)+nofactor(~ganglion)+nicov(~poly(pos,2))+nocov(~poly(pos,2))+absdiff(~pos)+mutual.rds")
tmp <- best_m$prediction_matrices$from_er
colnames(tmp) <- rownames(tmp) <- cellinfo$name
p_best_model_motor <- create_heatmap(best_m$prediction_matrices$from_er[idx,idx], "darkblue")
p_data_motor <- create_heatmap(as.matrix(g), "darkblue")

N_LL = 10000
k = job_idx+1

###### brute force
set.seed(1)

classes = sample(1:k, network.size(g), replace=T)
original_classes <- classes

g%v%"latent" <- classes
current_model <- ergm(g ~ nodemix("latent") + offset(absdiff("pos")), verbose = F, offset.coef = c(off))

#maximize LL
ll <- as.double(logLik(current_model))
all_lls <- rep(0, N_LL)
for (n in 1:N_LL) {
  print(n)
  all_lls[n] <- ll
  g%v%"latent" <- classes
  best_classes <- classes
  
  random_neuron <- sample(1:network.size(g),1)
  current_class <- classes[random_neuron]
  
  for (random_class in 1:k) {
    if (random_class!=current_class) {
      new_classes <- classes
      new_classes[random_neuron] <- random_class
      g%v%"latent" <- new_classes
      new_model <- ergm(g ~ nodemix("latent") + offset(absdiff("pos")), verbose = F, offset.coef = c(off))
      # new_ll <- get_grade(new_model)
      new_ll <- as.double(logLik(new_model))
      if (new_ll > ll) {
        ll <- new_ll
        best_classes <- new_classes
      }
    }
  }
  classes <- best_classes
}
ll_classes <- classes

save.image(file=paste0("./offset_worm_motor_k_", k, ".RData"))
stop()

# p_optim <- data.frame(x=1:length(all_lls), y=all_lls) %>% ggplot(aes(x=x, y=y)) + 
#   geom_line() + theme_classic() + geom_point()
# p_optim

# g%v%"classes" <- ll_classes
# table(ll_classes, g%v%"subtype_ext")

# #m_ll <- ergm(g ~ nodemix("classes") + absdiff("pos"))
# p_motor_mean <- create_heatmap(as.matrix(predict(m_ll, output="matrix", conditional=T)), "darkblue")
# p_motor_samples <- list()
# sims <- simulate(m_ll, nsim = 10)
# for (i in 1:10) {
#   p_motor_samples[[i]] <- create_heatmap(as.matrix(sims[[i]]), "darkblue")
# }

# save.image(file="./worm_motor_after_optim.RData")

# new_best_m_formula <- best_m$form
# new_types <- cellinfo$subtype
# new_types[which(cellinfo$name %in% idx)] <- as.character(ll_classes)
# full_g%v%"new_types" <- new_types
# new_form <- str_replace(str_replace(best_m$form, "g ~", "full_g ~"), "subtype", "new_types")
# control_list = snctrl(MCMC.prop = ~strat(attr = ~subtype, empirical = TRUE) + sparse, #MCMLE.termination = termination,
#                       # MCMC.burnin = 1e6, MCMC.interval = 1e4,
#                       seed = 123, 
#                       MPLE.type = "glm", main.method = "Stochastic-Approximation")
 
# stop()
