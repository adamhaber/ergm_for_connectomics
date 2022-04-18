library(forcats)
source("../utils/utils.R")

load("./allen_figures.RData",  temp_env <- new.env())
allen_env <- as.list(temp_env)

#allen figures
cmap = "darkgreen"
for (p in names(allen_env)) {
  if (startsWith(p, "p_")) {
    print(p)
    tryCatch(ggsave(paste0("./individual_figures/allen/pdfs/", p, ".pdf"), plot=allen_env[[p]], height = 9, width = 9, units = "cm"),
             error = function(e) print("fail"))
    tryCatch(ggsave(paste0("./individual_figures/allen/pngs/", p, ".png"), plot=allen_env[[p]], height = 9, width = 9, units = "cm", bg = "transparent"),
             error = function(e) print("fail"))
  }
}