library(forcats)
source("../utils/utils.R")

load("./ob_figures.RData",  temp_env <- new.env())
ob_env <- as.list(temp_env)

new_df <- ob_env$df %>% arrange(job_id)
single_feature_model_job_ids <- new_df %>% filter(n_features==2) %>% select(job_id) %>% unlist
single_feature_model_job_lls <- c()
single_feature_model_job_aucs <- c()
all_minus_feature_model_job_ids <- c()
all_minus_feature_model_job_lls <- c()
all_minus_feature_model_job_aucs <- c()
single_feature_names <- c()
for (s in single_feature_model_job_ids) {
  comp_features <- 1-ob_env$terms_df[s,]
  single_feature_model_job_lls <- c(single_feature_model_job_lls,
                                    new_df %>% filter(job_id==s) %>% select(lls) %>% as.double)
  single_feature_model_job_aucs <- c(single_feature_model_job_aucs,
                                    new_df %>% filter(job_id==s) %>% select(auc_from_er) %>% as.double)
  single_feature_names <- c(single_feature_names,
                            names(ob_env$terms_df)[which(ob_env$terms_df[s,]==1)])
  for (j in 1:nrow(ob_env$terms_df)) {
    if (all(ob_env$terms_df[j,]==comp_features)) {
      all_minus_feature_model_job_ids = c(all_minus_feature_model_job_ids, j)
      all_minus_feature_model_job_lls <- c(all_minus_feature_model_job_lls,
                                        new_df %>% filter(job_id==j) %>% select(lls) %>% as.double)
      all_minus_feature_model_job_aucs <- c(all_minus_feature_model_job_aucs,
                                         new_df %>% filter(job_id==j) %>% select(auc_from_er) %>% as.double)
      
    }
  }
}
library(RColorBrewer)

full_model_auc <- ob_env$full_model$aucs$from_er
full_model_ll <- ob_env$full_model$ll

names(ob_env$terms_df)
cur_pal <- brewer.pal(9,"Paired")[c(8,7,9,6,5,2,3,4)]
data.frame(all_minus = full_model_ll-all_minus_feature_model_job_lls, 
           single = single_feature_model_job_lls,
           names = single_feature_names) %>% 
  # arrange(match(names, single_feature_names)) %>% 
  ggplot(aes(x=single, y=all_minus, color=names)) + geom_point(size=3) + theme_classic() + 
  scale_color_manual(values = cur_pal) + xlab("Feature importance (alone)") + 
  ylab("Feature importance (combined)") -> ob_env$p_complementary_model

ob_env$p_complementary_model_no_legend <- ob_env$p_complementary_model + theme(legend.position = 'none')
ob_env$p_complementary_model_no_legend


data.frame(all_minus = full_model_auc-all_minus_feature_model_job_aucs, 
           single = single_feature_model_job_aucs,
           names = single_feature_names) %>% 
  # arrange(match(names, single_feature_names)) %>% 
  ggplot(aes(x=single, y=all_minus, color=names)) + geom_point(size=3) + theme_classic() + 
  scale_color_manual(values = cur_pal) + xlab("Feature importance (alone)") + 
  ylab("Feature importance (combined)") -> ob_env$p_complementary_model_auc

ob_env$p_complementary_model_auc_no_legend <- ob_env$p_complementary_model_auc + theme(legend.position = 'none')
ob_env$p_complementary_model_auc_no_legend

ob_env$p_data <- create_heatmap(ob_env$best_m$g %>% as.matrix, "darkred", title="Data")
ob_env$single_feature_models
fig_counter = 1
for (p in ob_env$ps_means) {
  # ggsave(paste0("./individual_figures/ob/pdfs/fig_",fig_counter,".pdf") , plot=p + ggtitle("") + xlab("") + ylab(""))
  ggsave(paste0("./individual_figures/ob/pngs/fig_",fig_counter,".png") , plot=p + ggtitle("") + xlab("") + ylab(""), height = 7, width = 7, units = "cm")
  # ggsave(paste0("./individual_figures/ob/pngs/fig_",fig_counter,"_with_title.png") , plot=p + xlab("") + ylab(""))
  fig_counter = fig_counter + 1
}

fig_counter = 1
for (p in ob_env$ps_means_colored) {
  # ggsave(paste0("./individual_figures/ob/pdfs/colored_fig_",fig_counter,".pdf") , plot=p + ggtitle("") + xlab("") + ylab(""))
  ggsave(paste0("./individual_figures/ob/pngs/colored_fig_",fig_counter,".png") , plot=p + ggtitle("") + xlab("") + ylab(""), height = 7, width = 7, units = "cm")
  # ggsave(paste0("./individual_figures/ob/pngs/colored_fig_",fig_counter,"_with_title.png") , plot=p + xlab("") + ylab(""))
  fig_counter = fig_counter + 1
}

fig_counter = 1
for (p in ob_env$ps_single_sample_colored) {
  # ggsave(paste0("./individual_figures/ob/pdfs/sample_colored_fig_",fig_counter,".pdf") , plot=p + ggtitle("") + xlab("") + ylab(""))
  ggsave(paste0("./individual_figures/ob/pngs/sample_colored_fig_",fig_counter,".png") , plot=p + ggtitle("") + xlab("") + ylab(""), height = 7, width = 7, units = "cm")
  # ggsave(paste0("./individual_figures/ob/pngs/sample_colored_fig_",fig_counter,"_with_title.png") , plot=p + xlab("") + ylab(""))
  fig_counter = fig_counter + 1
}


s_corrs = 7
s_mses = 9
cmap = "darkred"
#ob figures
for (p in names(ob_env)) {
  if (startsWith(p, "p_")) {
    ggsave(paste0("./individual_figures/ob/pdfs/", p, ".pdf"), plot=ob_env[[p]], height = 9, width = 9, units = "cm")
    ggsave(paste0("./individual_figures/ob/pngs/", p, ".png"), plot=ob_env[[p]], height = 9, width = 9, units = "cm")
  }
}

ob_full_model_preds <- create_heatmap(ob_env$full_model$prediction_matrices$from_er, cmap, title="Synaptic probabilities", with_guide = "colourbar", lims=heatmap_lims) 
ggsave(paste0("./individual_figures/ob/pdfs/", "ob_full_model_preds", ".pdf"), plot=ob_full_model_preds, height = 9, width = 9, units = "cm")
ggsave(paste0("./individual_figures/ob/pngs/", "ob_full_model_preds", ".png"), plot=ob_full_model_preds, height = 9, width = 9, units = "cm")

ggsave(paste0("./individual_figures/ob/pdfs/", "p_complementary_model_no_legend", ".pdf"), plot=ob_env$p_complementary_model_no_legend, height = 9, width = 9, units = "cm")
ggsave(paste0("./individual_figures/ob/pngs/", "p_complementary_model_no_legend", ".png"), plot=ob_env$p_complementary_model_no_legend, height = 9, width = 9, units = "cm")

