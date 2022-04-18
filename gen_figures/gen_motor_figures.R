source("../utils/utils.R")
source("./celegans/load_worm_data.R")

cellinfo <- cellinfo %>% mutate(M = if_else(region=="M", "M", "E"))
idx = cellinfo %>% filter(type=="motorneuron") %>% filter(M=="M") %>% select(name) %>% unlist
motor_connectome <- connectome[idx, idx]
bin_connectome <- (motor_connectome > 0)

full_g <- g
g <- network(bin_connectome>0, vertex.attr = cellinfo %>% filter(name %in% idx))

best_m <- readRDS("./celegans/results/10_3_no_subtype_ext/res_939_edges+nmix(~subtype)+nicov(~poly(time,2))+nocov(~poly(time,2))+absdiff(~pos)+mutual.rds")
tmp <- best_m$prediction_matrices$from_er
colnames(tmp) <- rownames(tmp) <- cellinfo$name
p_best_model_motor <- create_heatmap(best_m$prediction_matrices$from_er[idx,idx], "darkblue", with_guide = "colourbar")
p_data_motor <- create_heatmap(as.matrix(g), "darkblue")

aucs <- c()
lls <- c()
figs <- list()
for (k in 2:10) {
  l <- list()
  load(paste0("./worm_motor_k_",k,".RData"))  

  class_df = data.frame(Latent_type = as.factor(ll_classes),
                        # orig = as.factor(g%v%"subtype_ext"),
                        Neuron_class = as.factor(g%v%"first_two"),
                        loc = g%v%"pos")
  ggplot(class_df, aes(x=loc, y=Neuron_class)) + 
    facet_wrap(vars(Latent_type)) +
    geom_point(size=5) +
    # theme_classic() + 
    ylab("Neuron class") +
    xlab("A-P position") +
    theme_bw() +
    scale_color_brewer(palette = "Paired") -> p_latents
  p_latents
  l$p_latents <- p_latents
  
  p_optim <- data.frame(x=1:length(all_lls), y=all_lls) %>% ggplot(aes(x=x, y=y)) +
    geom_line() + theme_classic() + geom_point()
  p_optim
  
  g%v%"classes" <- ll_classes
  # m_ll <- ergm(g ~ nodemix("classes") + offset(absdiff("pos")), offset.coef = -13)
  m_ll <- ergm(g ~ nodemix("classes") + absdiff("pos"))
  mat <- predict(m_ll, conditional=T, output="matrix")
  auc <- get_metrics(mat, g, "auc")
  l$mat <- mat
  l$classes <- ll_classes
  ll <- logLik(m_ll)
  figs[[k]] <- l
  lls <- c(lls, ll)
  aucs <- c(aucs, auc)
  
}
data.frame(k = 2:10,
           auc = aucs,
           lls = lls) %>% 
  ggplot(aes(x=k, y=lls)) + geom_line() + theme_classic() -> p_motor_lls

data.frame(k = 2:10,
           auc = aucs,
           lls = lls) %>% 
  ggplot(aes(x=k, y=aucs)) + geom_line() + theme_classic() -> p_motor_aucs

p_model <- create_heatmap(figs[[7]]$mat, "darkblue", with_guide = "colourbar")
m2 <- ergm(g ~ nodemix(~first_two) + absdiff(~pos))
mat <- predict(m2, conditional=T, output="matrix")
auc <- get_metrics(mat, g, "auc")
tmp <- table(figs[[7]]$classes, g%v%"first_two")
fit_cols <- hclust(as.dist(1-cor(tmp)))
fit_rows <- hclust(as.dist(1-cor(t(tmp))))
mat <- tmp[fit_rows$order,fit_cols$order]
p_latent_table <- create_heatmap(mat, "black") + xlab("") + ylab("")

#p_optim
ggsave(paste0("./individual_figures/motor/pngs/", "p_best_model_motor", ".png"), plot=p_best_model_motor, height = 5, width = 5, units = "cm")
ggsave(paste0("./individual_figures/motor/pngs/", "p_data_motor", ".png"), plot=p_data_motor, height = 5, width = 5, units = "cm")
ggsave(paste0("./individual_figures/motor/pngs/", "p_model", ".png"), plot=p_model, height = 5, width = 5, units = "cm")
ggsave(paste0("./individual_figures/motor/pdfs/", "p_latents", ".pdf"), plot=p_latent_table, height = 5, width = 5, units = "cm")
