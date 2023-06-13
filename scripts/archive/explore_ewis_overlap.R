library(tidyverse)


variant_df <- read_tsv("../data/processed/target_gene_variant_list.tsv") %>%
  select(rsID=`Variant name`, gene=`Gene name`)
ewis_df <- read_csv("../data/raw/GxE_node19_filtered.csv",
                    col_types=cols(Locus="c"))

merged_df <- inner_join(variant_df, ewis_df, by=c("rsID"="SNP"))
