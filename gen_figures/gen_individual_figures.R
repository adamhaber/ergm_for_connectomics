library(forcats)
source("../utils/utils.R")

load("./ob_figures.RData",  temp_env <- new.env())
ob_env <- as.list(temp_env)

load("./allen_figures.RData",  temp_env <- new.env())
allen_env <- as.list(temp_env)

load("./worm_figures.RData",  temp_env <- new.env())
worm_env <- as.list(temp_env)


ob_env$p_data <- create_heatmap(ob_env$best_m$g %>% as.matrix, "darkred", title="Data")

fig_counter = 1
for (p in ob_env$ps_means) {
  # ggsave(paste0("./individual_figures/ob/pdfs/fig_",fig_counter,".pdf") , plot=p + ggtitle("") + xlab("") + ylab(""))
  ggsave(paste0("./individual_figures/ob/pngs/fig_",fig_counter,".png") , plot=p + ggtitle("") + xlab("") + ylab(""), height = 2, width = 2, units = "cm", bg = "transparent")
  # ggsave(paste0("./individual_figures/ob/pngs/fig_",fig_counter,"_with_title.png") , plot=p + xlab("") + ylab(""))
  fig_counter = fig_counter + 1
}

fig_counter = 1
for (p in ob_env$ps_means_colored) {
  # ggsave(paste0("./individual_figures/ob/pdfs/colored_fig_",fig_counter,".pdf") , plot=p + ggtitle("") + xlab("") + ylab(""))
  ggsave(paste0("./individual_figures/ob/pngs/colored_fig_",fig_counter,".png") , plot=p + ggtitle("") + xlab("") + ylab(""), height = 2, width = 2, units = "cm", bg = "transparent")
  # ggsave(paste0("./individual_figures/ob/pngs/colored_fig_",fig_counter,"_with_title.png") , plot=p + xlab("") + ylab(""))
  fig_counter = fig_counter + 1
}

fig_counter = 1
for (p in ob_env$ps_single_sample_colored) {
  # ggsave(paste0("./individual_figures/ob/pdfs/sample_colored_fig_",fig_counter,".pdf") , plot=p + ggtitle("") + xlab("") + ylab(""))
  ggsave(paste0("./individual_figures/ob/pngs/sample_colored_fig_",fig_counter,".png") , plot=p + ggtitle("") + xlab("") + ylab(""), height = 2, width = 2, units = "cm", bg = "transparent")
  # ggsave(paste0("./individual_figures/ob/pngs/sample_colored_fig_",fig_counter,"_with_title.png") , plot=p + xlab("") + ylab(""))
  fig_counter = fig_counter + 1
}


s_corrs = 7
s_mses = 9
cmap = "darkred"
#ob figures
for (p in names(ob_env)) {
  if (startsWith(p, "p_")) {
    ggsave(paste0("./individual_figures/ob/pdfs/", p, ".pdf"), plot=ob_env[[p]], height = 5, width = 5, units = "cm")
    ggsave(paste0("./individual_figures/ob/pngs/", p, ".png"), plot=ob_env[[p]], height = 5, width = 5, units = "cm", bg = "transparent")
  }
}

ob_full_model_preds <- create_heatmap(ob_env$full_model$prediction_matrices$from_er, cmap, title="Synaptic probabilities")
ggsave(paste0("./individual_figures/ob/pdfs/", "ob_full_model_preds", ".pdf"), plot=ob_full_model_preds, height = 5, width = 5, units = "cm")
ggsave(paste0("./individual_figures/ob/pngs/", "ob_full_model_preds", ".png"), plot=ob_full_model_preds, height = 5, width = 5, units = "cm", bg = "transparent")


#allen figures
cmap = "darkgreen"
for (p in names(allen_env)) {
  if (startsWith(p, "p_")) {
    print(p)
    tryCatch(ggsave(paste0("./individual_figures/allen/pdfs/", p, ".pdf"), plot=allen_env[[p]], height = 5, width = 5, units = "cm"),
             error = function(e) print("fail"))
    tryCatch(ggsave(paste0("./individual_figures/allen/pngs/", p, ".png"), plot=allen_env[[p]], height = 5, width = 5, units = "cm", bg = "transparent"),
             error = function(e) print("fail"))
  }
}

non_sender <- sapply(worm_env$filenames, function(x) str_detect(x, "sender"))
cv_auc_err = worm_env$cv_auc_err[!non_sender]
#worm figures
cmap = "darkblue"
for (p in names(worm_env)) {
  if (startsWith(p, "p_")) {
    print(p)
    tryCatch(ggsave(paste0("./individual_figures/worm/pdfs/", p, ".pdf"), plot=worm_env[[p]], height = 5, width = 5, units = "cm"),
             error = function(e) print(e))
    tryCatch(ggsave(paste0("./individual_figures/worm/pngs/", p, ".png"), plot=worm_env[[p]], height = 5, width = 5, units = "cm", bg = "transparent"),
             error = function(e) print(e))
  }
}

