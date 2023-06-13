library(tidyverse)
md_components <- scan("../data/processed/meddiet_components.txt", what=character())
a <- read_tsv("../data/processed/gwis/hscrp_log_MDcomponents_gwis_chr16")
b <- filter(a, POS > (30125426 - 1e5), POS < (30125426 + 1e5))
for (mc in md_components) {
  mc_beta <- paste0("Beta_G-", mc)
  mc_var <- paste0("Var_Beta_G-", mc)
  b[[paste0(mc, "_z")]] <- b[[mc_beta]] / sqrt(b[[mc_var]])
  print(mc)
  print(summary(abs(b[[paste0(mc, "_z")]])))
}