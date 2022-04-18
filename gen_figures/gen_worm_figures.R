library(forcats)
source("../utils/utils.R")

load("./worm_figures.RData",  temp_env <- new.env())
worm_env <- as.list(temp_env)

non_sender <- sapply(worm_env$filenames, function(x) str_detect(x, "sender"))
cv_auc_err = worm_env$cv_auc_err[!non_sender]
#worm figures
cmap = "darkblue"
for (p in names(worm_env)) {
  if (startsWith(p, "p_")) {
    print(p)
    tryCatch(ggsave(paste0("./individual_figures/worm/pdfs/", p, ".pdf"), plot=worm_env[[p]], height = 9, width = 9, units = "cm"),
             error = function(e) print(e))
    tryCatch(ggsave(paste0("./individual_figures/worm/pngs/", p, ".png"), plot=worm_env[[p]], height = 9, width = 9, units = "cm", bg = "transparent"),
             error = function(e) print(e))
  }
}

