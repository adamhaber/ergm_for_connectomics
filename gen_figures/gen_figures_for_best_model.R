"""
This file creates the figures for the best model for each dataset.
Putting all the code in one file makes it easier to standardize the figures across datasets.
"""

roc_objs <- list()
roc_objs$full <- roc(non_diag(full_model$g %>% as.matrix), 
                     non_diag(full_model$prediction_matrices$from_er %>% as.matrix)) 
roc_objs$chosen <- roc(non_diag(best_m$g %>% as.matrix), 
                       non_diag(best_m$prediction_matrices$from_er %>% as.matrix)) 
roc_objs$er <- roc(non_diag(er_model$g %>% as.matrix), 
                       non_diag(er_model$prediction_matrices$from_er %>% as.matrix)) 
best_m_auc <- paste0("AUROC = ", round(best_m$aucs$from_er,2))
ggroc(roc_objs) + theme_classic()  + coord_fixed() + 
  scale_color_manual(values = c("black", cmap, "grey")) + 
  annotate("text", x=0.75, y=1, size=5, label=best_m_auc) -> p_rocs_full_chosen

if (cmap=="darkred") {
  brew = "Reds"
} else if (cmap=="darkblue") {
  brew = "Blues"
} else {
  brew = "Greens"
}


chosen_full_preds_df <- data.frame(full = non_diag(full_model$prediction_matrices$from_er %>% as.matrix), 
                       chosen = non_diag(best_m$prediction_matrices$from_er %>% as.matrix), 
                        Synapse = as.factor(non_diag(full_model$g %>% as.matrix))) %>% 
  mutate(Synapse = if_else(Synapse==0, "No", "Yes"))

chosen_full_preds_df_sample <- chosen_full_preds_df %>% 
  filter(chosen>0, full>0) %>% 
  group_by(Synapse) %>%
  sample_n(50) %>% 
  mutate(chosen_err = sqrt(chosen*(1-chosen)/500),
         full_err = sqrt(full*(1-full)/500))

ggplot(chosen_full_preds_df_sample, 
       aes(x=full, xmin = full-full_err, xmax = full+full_err,
           y=chosen, ymin = chosen-chosen_err, ymax = chosen+chosen_err,
           color=Synapse, fill=Synapse)) + 
  geom_errorbar(alpha=1) +
  geom_errorbarh(alpha=1) +
  geom_point(alpha=1, shape=21) +
  # scale_color_brewer(palette=brew) +
  scale_color_manual(values=c("black", cmap)) +
  scale_fill_manual(values=c("white", cmap)) +
  theme_classic() + 
  xlab("Full model") + 
  ylab("Chosen model") +
  coord_fixed() +
  geom_abline(intercept = 0, slope=1) -> p_full_vs_chosen_probs
p_full_vs_chosen_probs

p_full_vs_chosen_probs_log <- p_full_vs_chosen_probs +
  scale_x_log10(limits=c(5e-4, 0.5)) +
  scale_y_log10(limits=c(5e-4, 0.5)) 
p_full_vs_chosen_probs_log

p_model_predictions <- create_heatmap(best_m$prediction_matrices$from_er, cmap, title="Synaptic probabilities", with_guide = "colourbar", lims = heatmap_lims)
p_model_predictions
p_model_predictions_no_guide <- p_model_predictions + theme(legend.position = 'none')

p_model_sample <- create_heatmap(best_m$closest$from_er %>% as.network %>% as.matrix, cmap, title="Single sample")
p_model_sample

p_data <- create_heatmap(g %>% as.matrix, cmap, title="Data")
p_data

er_motifs = er_model$motifs

er_df <- data.frame(
  mean = apply(er_motifs,1,mean),
  gs = gs_motifs,
  q5 = apply(er_motifs,1,function(x) quantile(x, .05)),
  q95 = apply(er_motifs,1,function(x) quantile(x, .95)),
  type = "Erdos-Reyni",
  name = rownames(er_motifs)) 

model_motifs = (best_m$motifs)
model_df <- data.frame(
  gs = gs_motifs,
  mean = apply(model_motifs,1,mean),
  q5 = apply(model_motifs,1,function(x) quantile(x, .05)),
  q95 = apply(model_motifs,1,function(x) quantile(x, .95)),
  type = "ERGM",
  name = rownames(model_motifs)) 

motifs_df <- rbind(model_df, er_df)

# motifs with ER
ggplot(data = motifs_df %>% mutate(mean = mean, q5 = q5, q95 = q95) %>% filter(gs>0),
       aes(x=gs, y=mean, ymax=q95, ymin=q5,colour=type, fill=type)) +
  geom_point(size=3) + 
  geom_errorbar() +
  scale_colour_manual("Model", values = c(str_replace(cmap,"dark",""), "grey60")) +
  scale_fill_manual("Model", values = c(str_replace(cmap,"dark",""), NA)) +
  scale_alpha(c(1,.1))+ 
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_classic() + geom_abline(intercept = 0, slope=1) +
  xlab("Motifs counts [data]") + 
  ylab("Motifs counts [model]") + 
  coord_fixed(
    xlim=c(0.1*min(gs_motifs[gs_motifs>0]), 2*max(gs_motifs[gs_motifs>0])),
    ylim=c(0.1*min(gs_motifs[gs_motifs>0]), 2*max(gs_motifs[gs_motifs>0]))
    )-> p_motifs_with_er
p_motifs_with_er


# motifs without ER
ggplot(data = model_df %>% mutate(mean = mean, q5 = q5, q95 = q95) %>% filter(gs>0),
       aes(x=gs, y=mean, ymax=q95, ymin=q5,colour=type, fill=type)) +
  geom_point(size=3) + 
  geom_errorbar() +
  scale_colour_manual("Model", values = c(str_replace(cmap,"dark",""), "grey60")) +
  scale_fill_manual("Model", values = c(str_replace(cmap,"dark",""), NA)) +
  scale_alpha(c(1,.1))+ 
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  theme_classic() + geom_abline(intercept = 0, slope=1) +
  xlab("Motifs counts [data]") + 
  ylab("Motifs counts [model]") + 
  coord_fixed(
    xlim=c(0.1*min(gs_motifs[gs_motifs>0]), 2*max(gs_motifs[gs_motifs>0])),
    ylim=c(0.1*min(gs_motifs[gs_motifs>0]), 2*max(gs_motifs[gs_motifs>0]))
  )-> p_motifs
p_motifs

#degree dist
get_idegs <- function(nw) {
  d <- (intergraph::asIgraph(nw %>% as.network) %>% degree_distribution(mode="in"))*network.size(g)
  if (length(d)>100) {
    d = d[1:100]
  }
  return(c(d, rep(0, 100-length(d))))
}
dists <- data.frame(sapply(samples, get_idegs)) %>%  mutate(id=1:100)
er_dists <- data.frame(sapply(er_samples, get_idegs)) %>%  mutate(id=1:100)

degdist_df <- melt(dists, id.vars = "id")
er_degdist_df <- melt(er_dists, id.vars = "id")

degdist_df <- degdist_df %>% group_by(id) %>% 
  summarise(mean = mean(value),
            q5 = quantile(value, .05),
            q95 = quantile(value, .95)) %>% 
  select(mean, q5, q95, id) %>% 
  filter(id <= degree_dist_upper_bound) %>% 
  mutate(type="ERGM")

er_degdist_df <- er_degdist_df %>% group_by(id) %>%
  summarise(mean = mean(value),
            q5 = quantile(value, .05),
            q95 = quantile(value, .95)) %>%
  select(mean, q5, q95, id) %>%
  filter(id <= degree_dist_upper_bound)%>%
  mutate(type="Erdos-Reyni")

combined_degdist_df = rbind(er_degdist_df, degdist_df)

real_degree_dist <- (intergraph::asIgraph(g) %>% degree_distribution(mode="in"))*network.size(g)
if (length(real_degree_dist)>100) {
  real_degree_dist = real_degree_dist[1:100]
}
real_degree_dist <- c(real_degree_dist, rep(0, 100-length(real_degree_dist)))
names(real_degree_dist) <- NULL
ideg_real_df <- data.frame(data=real_degree_dist, id=1:100)

ggplot() +
  geom_ribbon(data = degdist_df, aes(x=id, ymin=q5, ymax=q95), fill=str_remove(cmap, "dark"), alpha=0.3) +
  # scale_fill_manual("", values=c("grey", cmap)) +
  geom_line(data = ideg_real_df %>% filter(id <= degree_dist_upper_bound), aes(x=id, y=data), size=0.4) +
  geom_line(data = degdist_df, aes(x=id, y=mean), col=cmap, size=1) +
  # geom_line(data = er_degdist_df, aes(x=id, y=mean), linetype = "dashed", col="grey", size=1) +
  
  xlab("Degree") + ylab("Number of neurons") +
  theme_classic() + theme(legend.position='none') -> p_ideg_dist
p_ideg_dist

ggplot() +
  geom_ribbon(data = combined_degdist_df, aes(x=id, ymin=q5, ymax=q95, group=type, fill=type), alpha=0.3) +
  scale_fill_manual("", values=c("grey", str_remove(cmap, "dark"))) +
  geom_line(data = degdist_df, aes(x=id, y=mean), col=cmap, size=1) +
  geom_line(data = er_degdist_df, aes(x=id, y=mean), col="grey", size=1) +
  geom_line(data = ideg_real_df %>% filter(id <= degree_dist_upper_bound), aes(x=id, y=data), size=0.4) +
  xlab("Degree") + ylab("Number of neurons") +
  theme_classic() + theme(legend.position='none') -> p_ideg_dist_with_er
p_ideg_dist_with_er

library(smoother)
ggplot() +
  # scale_fill_manual("", values=c("grey", cmap)) + 
  geom_line(data = degdist_df, aes(x=id, y=mean), col=cmap, size=1) +
  geom_line(data = er_degdist_df, aes(x=id, y=mean), col="grey", size=1) +
  geom_line(data = ideg_real_df %>% filter(id <= degree_dist_upper_bound) %>% mutate(data = smth.gaussian(data, tails = T)), aes(x=id, y=data), size=1) +
  scale_y_log10() +
  xlim(5,0.8*degree_dist_upper_bound) +
  xlab("Degree") + ylab("Number of neurons") +
  theme_classic() + theme(legend.position='none') -> p_ideg_dist_log
p_ideg_dist_log

#odegree dist
get_odegs <- function(nw) {
  d <- (intergraph::asIgraph(nw %>% as.network) %>% degree_distribution(mode="out"))*network.size(g)
  if (length(d)>100) {
    d = d[1:100]
  }
  return(c(d, rep(0, 100-length(d))))
}
dists <- data.frame(sapply(samples, get_odegs)) %>%  mutate(id=1:100)
er_dists <- data.frame(sapply(er_samples, get_odegs)) %>%  mutate(id=1:100)

degdist_df <- melt(dists, id.vars = "id")
er_degdist_df <- melt(er_dists, id.vars = "id")

degdist_df <- degdist_df %>% group_by(id) %>% 
  summarise(mean = mean(value),
            q5 = quantile(value, .05),
            q95 = quantile(value, .95)) %>% 
  select(mean, q5, q95, id) %>% 
  filter(id <= degree_dist_upper_bound) %>% 
  mutate(type="ERGM")

er_degdist_df <- er_degdist_df %>% group_by(id) %>%
summarise(mean = mean(value),
          q5 = quantile(value, .05),
          q95 = quantile(value, .95)) %>%
select(mean, q5, q95, id) %>%
filter(id <= degree_dist_upper_bound)%>%
mutate(type="Erdos-Reyni")

combined_degdist_df = rbind(er_degdist_df, degdist_df)

real_degree_dist <- (intergraph::asIgraph(g) %>% degree_distribution(mode="out"))*network.size(g)
if (length(real_degree_dist)>100) {
  real_degree_dist = real_degree_dist[1:100]
}
real_degree_dist <- c(real_degree_dist, rep(0, 100-length(real_degree_dist)))
names(real_degree_dist) <- NULL
odeg_real_df <- data.frame(data=real_degree_dist, id=1:100)

ggplot() +
  geom_ribbon(data = degdist_df, aes(x=id, ymin=q5, ymax=q95), fill=str_remove(cmap, "dark"), alpha=0.3) +
  # scale_fill_manual("", values=c("grey", cmap)) +
  geom_line(data = degdist_df, aes(x=id, y=mean), col=cmap, size=1) +
  # geom_line(data = er_degdist_df, aes(x=id, y=mean), linetype = "dashed", col="grey", size=1) +
  geom_line(data = odeg_real_df %>% filter(id <= degree_dist_upper_bound), aes(x=id, y=data), size=0.4) +
  xlab("Degree") + ylab("Number of neurons") +
  theme_classic() + theme(legend.position='none') -> p_odeg_dist
p_odeg_dist

ggplot() +
  geom_ribbon(data = combined_degdist_df, aes(x=id, ymin=q5, ymax=q95, group=type, fill=type), alpha=0.3) +
  scale_fill_manual("", values=c("grey", str_remove(cmap, "dark"))) +
  geom_line(data = degdist_df, aes(x=id, y=mean), col=cmap, size=1) +
  geom_line(data = er_degdist_df, aes(x=id, y=mean), col="grey", size=1) +
  geom_line(data = odeg_real_df %>% filter(id <= degree_dist_upper_bound), aes(x=id, y=data), size=0.4) +
  xlab("Degree") + ylab("Number of neurons") +
  theme_classic() + theme(legend.position='none') -> p_odeg_dist_with_er
p_odeg_dist_with_er

# ggsave(str_c("./individual_figures/", animal_path, "/", animal, "_best_model_odegree_dist.png"), p_odeg_dist)
# ggsave(str_c("./individual_figures/", animal_path, "/", animal, "_best_model_odegree_dist.png"), p_odeg_dist)

ggplot() +
  # scale_fill_manual("", values=c("grey", cmap)) + 
  geom_line(data = degdist_df, aes(x=id, y=mean), col=cmap, size=1) +
  geom_line(data = er_degdist_df, aes(x=id, y=mean), col="grey", size=1) +
  geom_line(data = odeg_real_df %>% filter(id <= degree_dist_upper_bound) %>% mutate(data = smth.gaussian(data, tails = T)), aes(x=id, y=data), size=1) +
  scale_y_log10() +
  xlim(5,0.8*degree_dist_upper_bound) +
  xlab("Degree") + ylab("Number of neurons") +
  theme_classic() + theme(legend.position='none') -> p_odeg_dist_log
p_odeg_dist_log


#distance - directed
get_dists <- function(nw) {
  if (class(nw)=="matrix") {
    nw = nw %>% as.network(directed=T)
  }
  real_dist <- ergm.geodistdist(nw)
  real_dist <- c(real_dist[1:3], sum(real_dist[4:(length(real_dist)-1)]), real_dist[length(real_dist)])
  return(real_dist)
}
dists <- data.frame(sapply(samples, get_dists)) 
er_dists <- data.frame(sapply(er_samples, get_dists))

dists_melt <- melt(dists %>% mutate(id = c(1:3, "4+", "Inf")), id.vars = "id") 
er_dists_melt <- melt(er_dists %>% mutate(id = c(1:3, "4+", "Inf")), id.vars = "id")

dists_melt <- dists_melt %>% group_by(id) %>% 
  summarise(mean = mean(value),
            q5 = quantile(value, .05),
            q95 = quantile(value, .95),
            type = "ERGM") %>% 
  select(mean, q5, q95, type, id)

er_dists_melt <- er_dists_melt %>% group_by(id) %>%
summarise(mean = mean(value),
          q5 = quantile(value, .05),
          q95 = quantile(value, .95),
          type = "Erdos-Reyni") %>%
select(mean, q5, q95, type, id)

real_dist <- get_dists(g)

Gs_df <- data.frame(mean = real_dist,
                    q5 = real_dist,
                    q95 = real_dist,
                    type = "Data", 
                    id = dists_melt$id)

library(forcats)
ggplot(rbind(dists_melt, er_dists_melt, Gs_df) %>% filter(id>=1), aes(x=id, y=mean, fill=fct_relevel(type, c("Erdos-Reyni", "ERGM", "Data")))) +
  geom_bar(stat="identity", position=position_dodge(0.8), width=0.7) +
  geom_errorbar(aes(ymin=q5, ymax=q95), width=.5,
                position=position_dodge(.8)) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  geom_hline(yintercept = 1) +
  scale_fill_manual("legend", values = c("grey", str_replace(cmap, "dark", ""), "black")) +
  theme_classic() + theme(legend.position='none') + 
# (axis.text.x = element_text(angle = 90),
#         panel.background = element_rect(fill = "transparent", colour = NA),
#         plot.background = element_rect(fill = "transparent", colour = NA)) +
  ylab("Number of pairs") + xlab("Geometric distance") -> p_dist_of_pairs_with_er
p_dist_of_pairs_with_er

ggplot(rbind(dists_melt, Gs_df) %>% filter(id>=1), aes(x=id, y=mean, fill=fct_relevel(type, c("Erdos-Reyni", "ERGM", "Data")))) +
  geom_bar(stat="identity", position=position_dodge(0.8), width=0.7) +
  geom_errorbar(aes(ymin=q5, ymax=q95), width=.5,
                position=position_dodge(.8)) +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  geom_hline(yintercept = 1) +
  scale_fill_manual("legend", values = c(str_replace(cmap, "dark", ""), "black")) +
  theme_classic(base_size = 15) + theme(legend.position='none') + 
  # (axis.text.x = element_text(angle = 90),
  #         panel.background = element_rect(fill = "transparent", colour = NA),
  #         plot.background = element_rect(fill = "transparent", colour = NA)) +
  ylab("Number of pairs") + xlab("Geometric distance") -> p_dist_of_pairs
p_dist_of_pairs