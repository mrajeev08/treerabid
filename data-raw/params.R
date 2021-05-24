# Output parameters into list format ----

# packages
library(dplyr)
library(readr)
library(tidyr)

# data
dk_params <- read_csv("data-raw/DK_params.csv")
si_params <- read_csv("data-raw/SI_params.csv")
names(si_params)[names(si_params) %in% c("SI_ml", "SI2_ml")] <- c("SI_meanlog", "SI2_meanlog")

# output
dk_params %>%
  mutate(param_name = paste0(kernel_type, "_", par1name,
                             ifelse(dist == "weibull", "_weibull", "")),
         param_name2 = paste0(kernel_type, "_", par2name,
                              ifelse(dist == "weibull", "_weibull", ""))) %>%
 select(param_name, param_name2, param_value = par1est, par2est) %>%
 bind_rows(select(., param_name = param_name2, param_value = par2est)) %>%
 select(param_name, param_value) %>%
 pivot_wider(names_from = param_name, values_from = param_value) -> dk_pars

params_treerabid <- c(as.list(si_params), as.list(dk_pars))
usethis::use_data(params_treerabid, overwrite = TRUE)
