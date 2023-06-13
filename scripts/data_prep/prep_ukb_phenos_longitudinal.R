library(tidyverse)
library(data.table)


### Longitudinal data retrieval ------------------------------------------------

wrangle_long_data <- function(
    df, 
    field_ids, 
    merge_array_func = function(x) mean(x, na.rm=TRUE)
) {
  df %>%
    select(id = f.eid, matches(paste0("f\\.", field_ids, "\\.", collapse = "|"))) %>%
    filter(if_any(-id, ~ !is.na(.))) %>%
    pivot_longer(-id, names_to = "field_colname", values_to = "value") %>%
    mutate(field_colname = gsub("^f\\.", "", field_colname)) %>%
    mutate(split_cols = str_split(field_colname, "\\.")) %>%
    mutate(field = sapply(split_cols, "[[", 1),
           instance = sapply(split_cols, "[[", 2),
           array_idx = sapply(split_cols, "[[", 3)) %>%
    select(-field_colname, -split_cols) %>%
    group_by(id, field, instance) %>%
    summarise(value = merge_array_func(value), .groups = "drop") %>%
    mutate(field = names(field_ids)[match(field, field_ids)]) %>%
    pivot_wider(names_from = "field", values_from = "value")
}

fetch_long_data_old <- function(df, field_ids, merge_array_func=mean) {
  df %>%
    select(id = f.eid,
           matches(paste0("f.", field_ids, collapse="|"))) %>%
    pivot_longer(-id, names_to="field_colname", values_to="value") %>%
    mutate(field_colname = gsub("^f\\.", "", field_colname)) %>%
    separate(field_colname, into=c("field", "instance", "array_idx"), sep=".") %>%
    group_by(id, field, instance) %>%
    summarise(value = mean(value, na.rm=TRUE), .groups="drop") %>%
    mutate(field = names(field_ids)[match(field, field_ids)]) %>%
    pivot_wider(names_from="field", values_from="value")
}

### Basic variables ------------------------------------------------------------
  
base_pheno_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb10528.tab.gz",
                     data.table=FALSE, stringsAsFactors=FALSE)

ac_fields <- c(ac = 54, ac_date = 53)
ac_long_df <- base_pheno_df %>%
  select(f.eid, contains("f.54."), contains("f.53.")) %>%
  mutate(across(everything(), as.character)) %>%
  wrangle_long_data(ac_fields, merge_array_func = function(x) x[1]) %>%
  mutate(id = as.integer(id))

basic_fields <- c(sex = 31, age = 21022,  bmi = 21001, 
                  sbp = 4080, dbp = 4079, fasting_hrs = 74)
basic_phenos_long_df <- base_pheno_df %>%
  wrangle_long_data(basic_fields) %>%
  mutate(age_squared = age^2,
         ageBySex = age * sex)

withdrawn_consent <- scan("/humgen/florezlab/UKBB_app27892/withdraw27892_220_29_September_2022.txt", what=character())

saveRDS(ac_long_df, "ac_long_df.rds")

### Biomarkers -----------------------------------------------------------------

bm_fields <- c(
  alt = 30620, alb = 30600, apoB = 30640, hscrp = 30710, chol = 30690, glu = 30740, 
  hba1c = 30750, hdl = 30760, ldl = 30780, shbg = 30830, tg = 30870, vitD = 30890
)

biomarker_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb28679.tab.gz", 
                          data.table=FALSE, stringsAsFactors=FALSE)
biomarker_long_df <- biomarker_df %>%
  wrangle_long_data(bm_fields)

# Collect medication data for biomarker adjustments
drug_df <- base_pheno_df %>%
  select(id=f.eid, contains("f.20003.0"), contains("f.6177.0"))

statin_ids <- c(
  1140861958, 1140888594, 1140888648, 1141146234, 1141192410, 1140861922, 1141146138
)
drug_df$num_statins <- 0
for (statin in statin_ids) {
  drug_df$num_statins <- drug_df$num_statins + rowSums(drug_df[, grep("20003", names(drug_df), value=TRUE)] == statin, na.rm=TRUE)
}
drug_df$bp_med <- rowSums(drug_df[, (
  grep("6177", names(drug_df), value=TRUE)
)] == 2, na.rm=TRUE)
drug_df <- select(drug_df, id, num_statins, bp_med)

# Make medication-based biomarker adjustments
biomarker_long_df <- left_join(biomarker_long_df, 
                               select(drug_df, id, num_statins),
                               by="id")  # Not currently joining by instance!
statin_adj_bms <- c("chol", "ldl", "apoB")
statin_adj_factors <- c(
  chol = 0.749,
  ldl = 0.684,
  apoB = 0.719
)
for (bm in statin_adj_bms) {
  adj_factor <- ifelse(biomarker_long_df$num_statins > 0, statin_adj_factors[bm], 1)
  biomarker_long_df[[paste0(bm, "_statinadj")]] <- biomarker_long_df[[bm]] / adj_factor
}

basic_phenos_long_df <- left_join(basic_phenos_long_df, 
                                  select(drug_df, id, bp_med), 
                                  by="id")
bp_adj_factors <- c(sbp = 15, dbp = 10)
for (bm in c("sbp", "dbp")) {
  adj_factor <- ifelse(basic_phenos_long_df$bp_med, bp_adj_factors[bm], 0)
  basic_phenos_long_df[[paste0(bm, "_medsadj")]] <- basic_phenos_long_df[[bm]] + adj_factor
}

saveRDS(biomarker_long_df, "biomarker_long_df.rds")
saveRDS(basic_phenos_long_df, "basic_phenos_long_df.rds")

### Outcomes for sample exclusion ----------------------------------------------

outcomes_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb38040.tab.gz",
                  data.table=FALSE, stringsAsFactors=FALSE) %>%
  select(id=f.eid, everything())
outcomes_df$diabetes <- rowSums(outcomes_df[, grepl("f\\.2443\\.0\\.", names(outcomes_df)), drop=FALSE] == 1, na.rm=T) > 0
outcomes_df$MI <- rowSums(outcomes_df[, grepl("f\\.6150\\.0\\.", names(outcomes_df)), drop=FALSE] == 1, na.rm=T) > 0
outcomes_df$angina <- rowSums(outcomes_df[, grepl("f\\.6150\\.0\\.", names(outcomes_df)), drop=FALSE] == 2, na.rm=T) > 0
cirrhosis_codes <- c(paste0("K70", 2:4), "K717", paste0("K74", 0:6))
cirrhosis_primary_ids <- c()
for (f in grep("f\\.41202\\.0\\.", names(outcomes_df), value=T)) {
  cirrhosis_primary_ids <- c(cirrhosis_primary_ids, outcomes_df$id[outcomes_df[[f]] %in% cirrhosis_codes])
}
cirrhosis_secondary_ids <- c()
for (f in grep("f\\.41204\\.0\\.", names(outcomes_df), value=T)) {
  cirrhosis_secondary_ids <- c(cirrhosis_secondary_ids, outcomes_df$id[outcomes_df[[f]] %in% cirrhosis_codes])
}
outcomes_df$pregnant <- rowSums(outcomes_df[, grepl("f\\.3140\\.0\\.", names(outcomes_df)), drop=FALSE] == 1, na.rm=T) > 0
cancer_tmp <- select(outcomes_df, 1, contains("f.40005."))
cancer <- cancer_tmp[, 1:7] %>%
  inner_join(select(filter(ac_long_df, instance == 0), id, ac_date), by="id")  # Add assessment center dates
cancer$ac_date = as.Date(cancer$ac_date)
for (i in 2:7) {
  x <- ifelse(abs(difftime(cancer[, i, drop=TRUE], cancer$ac_date, units="days")) <= 365, TRUE, FALSE)  # TRUE if cancer diagnosis within a year of assessment center visit
  cancer <- cbind(cancer, x)
}
cancer$cancer_within_1yearac = apply(cancer[, 9:14], 1, function(x) {
  ifelse(any(x == TRUE, na.rm=TRUE), TRUE, FALSE)
})
cancer[names(cancer) == "x"] <- NULL

outcomes_df <- outcomes_df %>%
  left_join(select(cancer, id, cancer_within_1yearac), by="id") %>%
  mutate(CHD = MI | angina,
         cirrhosis = id %in% c(cirrhosis_primary_ids, cirrhosis_secondary_ids)) %>%
  select(id, diabetes, CHD, cirrhosis, pregnant, cancer_within_1yearac)

saveRDS(outcomes_df, "outcomes_df.rds")

### Diet from 24HR -------------------------------------------------------------

# calc_diet_fields <- function(fieldIDs, df, coding=FALSE) {
#   # Given a list of fields constituting a food group:
#   # - Determine the set of 24HR that are valid for that food group
#   # - Recode the relevant variables based on their codings if necessary
#   # - Sum over all fields for that food group within each instance
#   # - Take the mean food group quantity over all instances
#   diet_field_df <- lapply(0:4, function(i) {
#     tcals_field <- paste0("f.100002.", i, ".0")  # Variable name for total calories in instance "i" 
#     typical_diet_field <- paste0("f.100020.", i, ".0")  # Variable name for typical diet in instance "i" 
#     valid_24hr <- (findInterval(df[[tcals_field]] / 4.18, c(600, 4800)) == 1) &
#       df[[typical_diet_field]] == 1
#     instance_fields <- paste0("f.", fieldIDs, ".", i, ".0")  # Variable names for all fields in instance "i"
#     instance_df <- df[, instance_fields, drop=FALSE]
#     if (coding) {  # Recode the variable if necessary (for food groups)
#       instance_df <- mutate_all(instance_df, ~codings[as.character(.)])
#     }
#     ifelse(valid_24hr,  # Sum over fields if valid 24HR, else NA
#            rowSums(instance_df, na.rm=TRUE),
#            NA)
#   }) %>%
#     setNames(paste0("instance", 0:4)) %>%
#     bind_cols()
#   diet_mean <- rowMeans(diet_field_df, na.rm=TRUE)
#   ifelse(is.nan(diet_mean), NA, diet_mean)
# }

calc_diet_fields <- function(field_ids, df, coding=FALSE) {
  # Given a list of fields constituting a food group:
  # - Determine the set of 24HR that are valid for that food group
  # - Recode the relevant variables based on their codings if necessary
  # - Sum over all fields for that food group
  valid_24hr <- (findInterval(df$TCALS_qc / 4.18, c(600, 4800)) == 1) &
    df$typical_diet_qc == 1
  df <- select(df, all_of(field_ids))
  if (coding) {  # Recode the variable if necessary (for food groups)
    df <- mutate_all(df, ~codings[as.character(.)])
  }
  ifelse(valid_24hr,  # Sum over fields if valid 24HR, else NA
         rowSums(df, na.rm=TRUE),
         NA)
}

winsorize <- function(x, SDs=3) {
  lims <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x[x < lims[1]] <- lims[1]
  x[x > lims[2]] <- lims[2]
  x
}

diet_qc_fields_list <- list(typical_diet = "100020", TCALS = "100002")
diet_nutrient_fields_list <- list(
  CHO = "100005", SUGARS = "100008", 
  FAT = "100004", SFA = "100006", PUFA = "100007",
  PRO = "100003", ALC = "100022", FIBER = "100009", 
  NUT_FE = "100011", NUT_K = "100016", NUT_FOL = "100014",
  NUT_MG = "100017", NUT_VITC = "100015", NUT_VITD = "100021"
)
diet_food_group_fields_list <- list(
  VEG =  as.character(seq(104060, 104380, 10)),  # Broad beans (104110), green beans (104120), and possibly others could arguably be included here instead of veg?
  LEGUMES = as.character(c(104000, 104010)),
  FRUIT = as.character(seq(104410, 104590, 10)),
  NUTS = as.character(seq(102410, 102440, 10)),
  FISH = as.character(seq(103150, 103230, 10)),
  WHGRAIN = as.character(c(101260, 102720, 102740, 100800, 100840, 100850)),
  REDPRMEAT = as.character(c(103010, 103020, 103030, 103040, 103050, 103070, 103080))
)
diet_fields_list <- c(diet_qc_fields_list, diet_nutrient_fields_list, 
                      diet_food_group_fields_list)
all_diet_vars <- c(names(diet_fields_list), "MUFA")
codings <- c(
  "1"=1, "2"=2, "3"=3, "4"=4, "5"=5, 
  "100"=1, "200"=2, "300"=3, "400"=4, "500"=5, "600"=6,
  "444"=0.25, "555"=0.5
)
# diet_vars <- c(
#   "TCALS", "CHO", "SUGARS", "FAT", "SFA", "MUFA", "PUFA", "PRO", "ALC", "FIBER",
#   "NUT_FE", "NUT_K", "NUT_FOL", "NUT_MG", "NUT_VITC", "NUT_VITD",
#   "VEG", "LEGUMES", "FRUIT", "NUTS", "FISH", "WHGRAIN", "REDPRMEAT"
# )

diet_24hr_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb22861.tab.gz", 
              data.table=FALSE, stringsAsFactors=FALSE)

# Do the 24HR wrangling group-by-group to keep memory usage under control
diet_24hr_qc_fields <- setNames(c("100020", "100002"), 
                                c("typical_diet_qc", "TCALS_qc"))
diet_24hr_qc_df <- diet_24hr_df %>%  # Get fields for "typical diet" and total calories to be used for QC
  wrangle_long_data(diet_24hr_qc_fields)
  
diet_24hr_long_list <- lapply(names(diet_fields_list), function(grp) {
  print(grp)
  grp_fields <- diet_fields_list[[grp]]
  is_food_group <- (grp %in% names(diet_food_group_fields_list))
  diet_24hr_df %>%
    filter(if_any(contains("f.100002."), ~ !is.na(.))) %>%
    wrangle_long_data(setNames(grp_fields, grp_fields)) %>%
    inner_join(diet_24hr_qc_df, by=c("id", "instance")) %>%
    mutate(!!grp := calc_diet_fields(grp_fields, ., coding = is_food_group)) %>%
    select(id, instance, !!grp)
})
diet_24hr_long_df <- reduce(diet_24hr_long_list, function(x, y) {
  full_join(x, y, by=c("id", "instance"))
}) %>%
  # inner_join(diet_24hr_qc_df, by=c("id", "instance")) %>%
  mutate(MUFA = FAT - SFA - PUFA,
         TCALS = TCALS / 4.18) %>%  # Energy from kJ to kcals
  mutate_at(vars(CHO, PRO), ~. * 4) %>%  # Nutrients from g to kcals (other than alcohol)
  mutate_at(vars(FAT, SFA, MUFA, PUFA), ~. * 9) %>%
  mutate_at(vars(all_of(all_diet_vars)), winsorize) %>%
  select(id, instance, all_of(all_diet_vars))

saveRDS(diet_24hr_long_df, "diet_24hr_long_df.rds")

# diet_fields <- c(unlist(diet_group_field_list), "100020")
# names(diet_fields) <- diet_fields
# diet_24hr_long_df <- diet_24hr_df %>%
#   filter(if_any(contains("f.100002."), ~ !is.na(.))) %>%
#   wrangle_long_data(diet_fields) %>%
#   mutate(TCALS = calc_diet_fields(diet_group_field_list$TCALS, .),
#          CHO = calc_diet_fields(diet_group_field_list$CHO, .),
#          SUGARS = calc_diet_fields(diet_group_field_list$SUGARS, .),
#          FAT = calc_diet_fields(diet_group_field_list$FAT, .),
#          SFA = calc_diet_fields(diet_group_field_list$SFA, .),
#          PUFA = calc_diet_fields(diet_group_field_list$PUFA, .),
#          PRO = calc_diet_fields(diet_group_field_list$PRO, .),
#          ALC = calc_diet_fields(diet_group_field_list$ALC, .),
#          FIBER = calc_diet_fields(diet_group_field_list$FIBER, .),
#          NUT_FE = calc_diet_fields(diet_group_field_list$NUT_FE, .),
#          NUT_K = calc_diet_fields(diet_group_field_list$NUT_K, .),
#          NUT_FOL = calc_diet_fields(diet_group_field_list$NUT_FOL, .),
#          NUT_MG = calc_diet_fields(diet_group_field_list$NUT_MG, .),
#          NUT_VITC = calc_diet_fields(diet_group_field_list$NUT_VITC, .),
#          NUT_VITD = calc_diet_fields(diet_group_field_list$NUT_VITD, .),
#          VEG = calc_diet_fields(diet_group_field_list$VEG, ., coding=TRUE),
#          LEGUMES = calc_diet_fields(diet_group_field_list$LEGUMES, ., coding=TRUE),
#          FRUIT = calc_diet_fields(diet_group_field_list$FRUIT, ., coding=TRUE),
#          NUTS = calc_diet_fields(diet_group_field_list$NUTS, ., coding=TRUE),
#          FISH = calc_diet_fields(diet_group_field_list$FISH, ., coding=TRUE),
#          WHGRAIN = calc_diet_fields(diet_group_field_list$WHGRAIN, ., coding=TRUE),
#          REDPRMEAT = calc_diet_fields(diet_group_field_list$REDPRMEAT, ., coding=TRUE)) %>%
#   mutate(MUFA = FAT - SFA - PUFA,
#          TCALS = TCALS / 4.18) %>%  # Energy from kJ to kcals
#   mutate_at(vars(CHO, PRO), ~. * 4) %>%  # Nutrients from g to kcals (other than alcohol)
#   mutate_at(vars(FAT, SFA, MUFA, PUFA), ~. * 9) %>%
#   mutate_at(vars(all_of(diet_vars)), winsorize) %>%
#   select(id, instance, all_of(names(diet_fields_list)))

ffq_cat_to_qt <- function(x) {
  case_when(  # Data-coding 100377
    x == 5 ~ 1,  # "Once or more daily"
    x == 4 ~ 5.5 / 7,  # "5-6 times a week"
    x == 3 ~ 3 / 7,  # "2-4 times a week"
    x == 2 ~ 1 / 7,  # "Once a week"
    x == 1 ~ 0.5 / 7,  # "Less than once a week"
    x == 0 ~ 0,  # "Never"
    TRUE ~ as.numeric(NA)
  )
}

ffq_fields <- c(oily_fish = 1329, nonoily_fish = 1339,
                 prmeat = 1349, poultry = 1359,
                 beef = 1369, lamb = 1379)
ffq_long_df <- base_pheno_df %>%
  wrangle_long_data(ffq_fields) %>%
  mutate(across(-c(id, instance), ffq_cat_to_qt))

saveRDS(ffq_long_df, "ffq_long_df.rds")

# ffq1_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb22861.tab.gz", 
#                  data.table=FALSE, stringsAsFactors=FALSE)
# ffq1_fields <- c(baked_beans = 104000, pulses = 104010)
# # N3FA = 26015, N6FA = 26016)
# ffq1_long_df <- ffq1_df %>%
#   wrangle_long_data(ffq1_fields) %>%
#   # mutate(across(c(N3FA, N6FA), ~. * 9)) %>%
#   mutate(across(-id, winsorize))

# supp_extra <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb38040.tab.gz", 
#                     data.table=FALSE, stringsAsFactors=FALSE) %>%
#   select(id = f.eid,
#          contains("f.6179.0")) %>%
#   mutate(fish_oil = as.integer(rowSums(across(-1, ~. == 1)) >= 1)) %>%
#   select(id, fish_oil)

supp_extra_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb22861.tab.gz", 
                       data.table=FALSE, stringsAsFactors=FALSE) 
supp_extra_long_df <- supp_extra_df %>%
  wrangle_long_data(c(fish_oil = 20084),
                    function(x) ifelse(all(is.na(x)), NA, 
                                       as.integer(any(x == 472, na.rm=TRUE))))

saveRDS(supp_extra_long_df, "supp_extra_long_df.rds")

calc_med_score <- function(diet_df) {
  # For additional details, see: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6842574/
  med_score_df <- diet_df %>%
    mutate(veg_score = VEG > median(VEG, na.rm=TRUE),
           legume_score = LEGUMES > median(LEGUMES, na.rm=TRUE),
           fruit_score = FRUIT > median(FRUIT, na.rm=TRUE),
           nut_score = NUTS > median(NUTS, na.rm=TRUE),
           fish_score = FISH > median(FISH, na.rm=TRUE),
           whgrain_score = WHGRAIN > median(WHGRAIN, na.rm=TRUE),
           mufa2sfa_score = MUFA2SFA > median(MUFA2SFA, na.rm=TRUE),
           redprmeat_score = REDPRMEAT < median(REDPRMEAT, na.rm=TRUE),
           alc_score = (ALC > 5) & (ALC < 25)) %>%
    mutate_at(vars(veg_score, legume_score, fruit_score, nut_score, fish_score, whgrain_score, mufa2sfa_score, redprmeat_score, alc_score),
              ~ifelse(is.na(.), mean(., na.rm=TRUE), .)) %>%
    mutate(mds = veg_score + legume_score + fruit_score + nut_score + fish_score + whgrain_score + mufa2sfa_score + redprmeat_score + alc_score)
  med_score_df$mds
}

full_diet_long_df <- diet_24hr_long_df %>%
  left_join(ffq_long_df, by=c("id", "instance")) %>%
  left_join(supp_extra_long_df, by=c("id", "instance")) %>%
  mutate(MUFA2SFA = MUFA / SFA) %>%
  mutate(mds = calc_med_score(.))

### Add genetic PCs and relatedness --------------------------------------------

gPC_df <- base_pheno_df %>%
  select(id=f.eid, contains("f.22009.0"), used_in_PCA=f.22020.0.0) %>%
  rename_with(.fn=~gsub("f.22009.0.", "gPC", .)) %>%
  mutate(unrelated = (used_in_PCA == 1)) %>%  # An unrelated subset was used in the central PCA
  select(id, contains("gPC"), unrelated)

### Merge and write "raw" phenotypes -------------------------------------------

# ac_long_df <- readRDS("ac_long_df.rds")
# biomarker_long_df <- readRDS("biomarker_long_df.rds")
# basic_phenos_long_df <- readRDS("basic_phenos_long_df.rds")
# diet_24hr_long_df <- readRDS("diet_24hr_long_df.rds")
# ffq_long_df <- readRDS("ffq_long_df.rds")
# supp_extra_long_df <- readRDS("supp_extra_long_df.rds")

phenos <- ac_long_df %>%
  inner_join(basic_phenos_long_df, by=c("id", "instance")) %>%
  inner_join(biomarker_long_df, by=c("id", "instance")) %>%
  left_join(outcomes_df, by="id") %>%
  left_join(full_diet_long_df, by=c("id", "instance")) %>%
  left_join(gPC_df, by="id") %>%
  filter(!(id %in% withdrawn_consent)) %>%
  mutate(id = format(id, scientific=FALSE)) %>%
  mutate(across(contains("gPC"), ~. * mds, .names="mdsBy{.col}"))

write_csv(phenos, "../data/processed/ukb_phenos_longitudinal_raw.csv")

### Phenotype processing and exclusions ----------------------------------------

logged_risk_factors <- c("alt", "tg", "hscrp")
risk_factors <- c(
  "bmi",
  "sbp", "dbp", "sbp_medsadj", "dbp_medsadj",
  "alt_log", 
  "chol", "ldl", "hdl", "apoB",
  "tg_log",
  "hba1c", "glu",
  "hscrp_log",
  "vitD",
  "chol_statinadj", "ldl_statinadj", "apoB_statinadj"
)
bin_diet_vars <- c("FISH", "LEGUMES", "NUTS", "WHGRAIN")

processed_phenos <- phenos %>%
  filter(!diabetes & !CHD & !cirrhosis & !cancer_within_1yearac & !pregnant) %>%
  mutate(across(one_of(logged_risk_factors), list(log=log))) %>%
  mutate(across(one_of(risk_factors), 
                ~ifelse(findInterval(., mean(., na.rm=TRUE) + c(-5, 5) * sd(., na.rm=TRUE)) != 1, 
                        as.numeric(NA), .))) %>%
  mutate(across(one_of(bin_diet_vars), ~as.integer(. > 0), .names="{.col}bin"))

### Write processed phenotypes -------------------------------------------------

write_csv(processed_phenos, "../data/processed/ukb_phenos_longitudinal.csv")

processed_phenos %>%
  filter(unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_longitudinal_unrelated.csv")

### Add Pan-UKBB data to generate European subset ------------------------------

anc_rel_df <- fread("/humgen/florezlab/UKBB_app27892/ukbreturn2442/all_pops_non_eur_pruned_within_pop_pc_covs_app27892.csv",
                    data.table=FALSE, stringsAsFactors=FALSE) %>%
  mutate(f.eid = as.character(f.eid),
         unrelated = !related_return2442) %>%
  select(id=f.eid, ancestry=pop_return2442, unrelated,
         one_of(paste0("PC", 1:10, "_return2442"))) %>%
  rename_at(vars(contains("PC")), ~gsub("_return2442", "", .)) %>%
  rename_at(vars(contains("PC")), ~gsub("PC", "gPC", .))

processed_phenos_panUKBB <- processed_phenos %>%
  select(-contains("gPC"), -unrelated) %>%
  inner_join(anc_rel_df, by="id") %>%
  mutate(across(contains("gPC"), ~. * mds, .names="mdsBy{.col}"))

processed_phenos_panUKBB %>%
  filter(ancestry == "EUR") %>%
  write_csv("../data/processed/ukb_phenos_longitudinal_EUR.csv")

processed_phenos_panUKBB %>%
  filter(ancestry == "EUR", unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_longitudinal_EUR_unrelated.csv")

processed_phenos_panUKBB %>%
  filter(ancestry == "AFR", unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_longitudinal_AFR_unrelated.csv")

processed_phenos_panUKBB %>%
  filter(ancestry == "AMR", unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_longitudinal_AMR_unrelated.csv")

processed_phenos_panUKBB %>%
  filter(ancestry == "CSA", unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_longitudinal_CSA_unrelated.csv")

processed_phenos_panUKBB %>%
  filter(ancestry == "EAS", unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_longitudinal_EAS_unrelated.csv")

processed_phenos_panUKBB %>%
  filter(ancestry == "MID", unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_longitudinal_MID_unrelated.csv")
