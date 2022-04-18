heatmap_lims = c(0,0.5)
source("./ob/ob_analyze_indep_results.R")
source("./gen_ob_figures.R")


rm(list=ls())
heatmap_lims = NULL
source("../utils/utils.R")
source("./allen/allen_analyze_indep_results.R")
# clear plots
source("./gen_allen_figures.R")

rm(list=ls())
heatmap_lims = NULL
source("../utils/utils.R")
source("./celegans/worm_analyze_indep_results.R")
source("./gen_worm_figures.R")

#gen motor results
source("./gen_motor_figures.R")
