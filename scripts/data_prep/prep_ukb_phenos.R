library(tidyverse)
library(data.table)


### Basic variables ------------------------------------------------------------

base_phenos <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb10528.tab.gz", 
                     data.table=FALSE, stringsAsFactors=FALSE) %>%
  mutate(sbp = (f.4080.0.0 + f.4080.0.1) / 2,
         dbp = (f.4079.0.0 + f.4079.0.1) / 2) %>%
  select(id = f.eid,
         sex = f.31.0.0,
         age = f.21022.0.0,
         fasting_hrs = f.74.0.0,
         bmi = f.21001.0.0,
         sbp, dbp,
         smoking = f.20116.0.0,
         ac = f.54.0.0,
         ac_date = f.53.0.0) %>%
  mutate(age_squared = age ** 2,
         ageBysex = age * sex)

withdrawn_consent <- scan("/humgen/florezlab/UKBB_app27892/withdraw27892_220_29_September_2022.txt", what=character())

### Biomarkers -----------------------------------------------------------------

ukb_biomarker_fields <- c(
  alt = 30620, alb = 30600, alp = 30610, apoA = 30630, apoB = 30640,
  ast = 30650, hscrp = 30710, Ca = 30680, chol = 30690, creatinine = 30700,
  cysC = 30720, bilirubin_dir = 30660, ggt = 30730, glu = 30740, hba1c = 30750,
  hdl = 30760, igf1 = 30770, ldl = 30780, lipA = 30790, oestradiol = 30800,
  phos = 30810, rheum_factor = 30820, shbg = 30830, tes = 30850,
  bilirubin_tot = 30840, protein_tot = 30860, tg = 30870, urate = 30880,
  urea = 30670, vitD = 30890
)
baseline_biomarker_fields <- setNames(paste0("f.", ukb_biomarker_fields, ".0.0"),
                                      names(ukb_biomarker_fields))

biomarker_phenos <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb28679.tab.gz", 
                          data.table=FALSE, stringsAsFactors=FALSE) %>%
  select(id = f.eid, 
         all_of(baseline_biomarker_fields))

ukb_cbc_fields <- c(
  baso_ct = 30160, eos_ct = 30150, lymph_ct = 30120, mono_ct = 30130, 
  neut_ct = 30140, platelet_ct = 30080, wbc_ct = 30000,
  mpv = 30100
)
baseline_cbc_fields <- setNames(paste0("f.", ukb_cbc_fields, ".0.0"),
                                      names(ukb_cbc_fields))

cbc_phenos <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb10528.tab.gz", 
                          data.table=FALSE, stringsAsFactors=FALSE) %>%
  select(id = f.eid, 
         all_of(baseline_cbc_fields)) %>%
  mutate(across(contains("_ct"), ~ifelse(. == 0, NA, .))) %>%
  mutate(
    aisi = neut_ct * platelet_ct * mono_ct / lymph_ct,
    mlr = mono_ct / lymph_ct,
    mpr = mpv / platelet_ct,
    nlr = neut_ct / lymph_ct,
    nlpr = neut_ct / (lymph_ct * platelet_ct),
    plr = platelet_ct / lymph_ct,
    sii = neut_ct * platelet_ct / lymph_ct,
    siri = neut_ct * mono_ct / lymph_ct
  ) %>%
  select(id, aisi, sii, siri)

drugs <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb10528.tab.gz",
               data.table=FALSE, stringsAsFactors=FALSE) %>%
  select(id=f.eid, contains("f.20003.0")) %>%
  right_join(select(biomarker_phenos, id), by="id")  # To align rows by ID
statin_ids <- c(
  1140861958, 1140888594, 1140888648, 1141146234, 1141192410, 1140861922, 1141146138
)
statin_adj_bms <- c("chol", "ldl", "apoB")
statin_adj_factors <- c(
  chol = 0.749,
  ldl = 0.684,
  apoB = 0.719
)
drugs$num_statins <- 0
for (statin in statin_ids) {
  drugs$num_statins <- drugs$num_statins + rowSums(drugs[, grep("20003", names(drugs), value=TRUE)] == statin, na.rm=TRUE)
}
for (bm in statin_adj_bms) {
  adj_factor <- ifelse(drugs$num_statins > 0, statin_adj_factors[bm], 1)
  biomarker_phenos[[paste0(bm, "_statinadj")]] <- biomarker_phenos[[bm]] / adj_factor
}


sr_meds <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb10528.tab.gz", 
                 data.table=FALSE, stringsAsFactors=FALSE) %>%
  select(id=f.eid, contains("f.6177.0")) %>%
  right_join(select(base_phenos, id), by="id")  # To align rows by ID
bp_adj_factors <- c(
  sbp = 15,
  dbp = 10
)
sr_meds$bp_med <- rowSums(sr_meds[, (
  grep("6177", names(sr_meds), value=TRUE)
)] == 2, na.rm=TRUE)
for (bm in c("sbp", "dbp")) {
  adj_factor <- ifelse(sr_meds$bp_med, bp_adj_factors[bm], 0)
  base_phenos[[paste0(bm, "_medsadj")]] <- base_phenos[[bm]] + adj_factor
}

### Outcomes for sample exclusion ----------------------------------------------

outcomes <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb38040.tab.gz",
                  data.table=FALSE, stringsAsFactors=FALSE) %>%
  select(id=f.eid, everything())
outcomes$diabetes <- rowSums(outcomes[, grepl("f\\.2443\\.0\\.", names(outcomes)), drop=FALSE] == 1, na.rm=T) > 0
outcomes$MI <- rowSums(outcomes[, grepl("f\\.6150\\.0\\.", names(outcomes)), drop=FALSE] == 1, na.rm=T) > 0
outcomes$angina <- rowSums(outcomes[, grepl("f\\.6150\\.0\\.", names(outcomes)), drop=FALSE] == 2, na.rm=T) > 0
cirrhosis_codes <- c(paste0("K70", 2:4), "K717", paste0("K74", 0:6))
cirrhosis_primary_ids <- c()
for (f in grep("f\\.41202\\.0\\.", names(outcomes), value=T)) {
  cirrhosis_primary_ids <- c(cirrhosis_primary_ids, outcomes$id[outcomes[[f]] %in% cirrhosis_codes])
}
cirrhosis_secondary_ids <- c()
for (f in grep("f\\.41204\\.0\\.", names(outcomes), value=T)) {
  cirrhosis_secondary_ids <- c(cirrhosis_secondary_ids, outcomes$id[outcomes[[f]] %in% cirrhosis_codes])
}
outcomes$pregnant <- rowSums(outcomes[, grepl("f\\.3140\\.0\\.", names(outcomes)), drop=FALSE] == 1, na.rm=T) > 0
cancer_tmp <- select(outcomes, 1, contains("f.40005."))
cancer <- inner_join(cancer_tmp[, 1:7], select(base_phenos, id, ac_date), by="id")  # Add assessment center dates
cancer$ac_date = as.Date(cancer$ac_date)
for (i in 2:7) {
  x <- ifelse(abs(difftime(cancer[, i, drop=TRUE], cancer$ac_date, units="days")) <= 365, TRUE, FALSE)  # TRUE if cancer diagnosis within a year of assessment center visit
  cancer <- cbind(cancer, x)
}
cancer$cancer_within_1yearac = apply(cancer[, 9:14], 1, function(x) {
  ifelse(any(x == TRUE, na.rm=TRUE), TRUE, FALSE)
})
cancer[names(cancer) == "x"] <- NULL

outcomes <- outcomes %>%
  left_join(select(cancer, id, cancer_within_1yearac), by="id") %>%
  mutate(CHD = MI | angina,
         cirrhosis = id %in% c(cirrhosis_primary_ids, cirrhosis_secondary_ids)) %>%
  select(id, diabetes, CHD, cirrhosis, pregnant, cancer_within_1yearac)

### Diet from 24HR -------------------------------------------------------------

fetch_diet_fields <- function(fieldIDs, df, coding=FALSE) {
  # Given a list of fields constituting a food group:
  # - Determine the set of 24HR that are valid for that food group
  # - Recode the relevant variables based on their codings if necessary
  # - Sum over all fields for that food group within each instance
  # - Take the mean food group quantity over all instances
  diet_field_df <- lapply(0:4, function(i) {
    tcals_field <- paste0("f.100002.", i, ".0")  # Variable name for total calories in instance "i" 
    typical_diet_field <- paste0("f.100020.", i, ".0")  # Variable name for typical diet in instance "i" 
    valid_24hr <- (findInterval(df[[tcals_field]] / 4.18, c(600, 4800)) == 1) &
      df[[typical_diet_field]] == 1
    instance_fields <- paste0("f.", fieldIDs, ".", i, ".0")  # Variable names for all fields in instance "i"
    instance_df <- df[, instance_fields, drop=FALSE]
    if (coding) {  # Recode the variable if necessary (for food groups)
      instance_df <- mutate_all(instance_df, ~codings[as.character(.)])
    }
    ifelse(valid_24hr,  # Sum over fields if valid 24HR, else NA
           rowSums(instance_df, na.rm=TRUE),
           NA)
  }) %>%
    setNames(paste0("instance", 0:4)) %>%
    bind_cols()
  diet_mean <- rowMeans(diet_field_df, na.rm=TRUE)
  ifelse(is.nan(diet_mean), NA, diet_mean)
}

check_num_valid_24hr <- function(df) {
  valid_24hr_df <- lapply(0:4, function(i) {
    tcals_field <- paste0("f.100002.", i, ".0")  # Variable name for total calories in instance "i" 
    typical_diet_field <- paste0("f.100020.", i, ".0")  # Variable name for typical diet in instance "i" 
    valid_24hr <- (findInterval(df[[tcals_field]] / 4.18, c(600, 4800)) == 1) &
      df[[typical_diet_field]] == 1
    valid_24hr
  }) %>%
    setNames(paste0("instance", 0:4)) %>%
    bind_cols()
  rowSums(valid_24hr_df, na.rm=TRUE)
}

winsorize <- function(x, SDs=3) {
  lims <- mean(x, na.rm=T) + SDs * c(-1, 1) * sd(x, na.rm=T)
  x[x < lims[1]] <- lims[1]
  x[x > lims[2]] <- lims[2]
  x
}

fields_list <- list(
  VEG =  seq(104060, 104380, 10),  # Broad beans (104110), green beans (104120), and possibly others could arguably be included here instead of veg?
  LEGUMES = c(104000, 104010),
  FRUIT = seq(104410, 104590, 10),
  NUTS = seq(102410, 102440, 10),
  FISH = seq(103150, 103230, 10),
  WHGRAIN = c(101260, 102720, 102740, 100800, 100840, 100850),
  REDPRMEAT = c(103010, 103020, 103030, 103040, 103050, 103070, 103080)
)
codings <- c(
  "1"=1, "2"=2, "3"=3, "4"=4, "5"=5, 
  "100"=1, "200"=2, "300"=3, "400"=4, "500"=5, "600"=6,
  "444"=0.25, "555"=0.5
)
diet_vars <- c(
  "TCALS", "CHO", "SUGARS", "FAT", "SFA", "MUFA", "PUFA", "PRO", "ALC", "FIBER",
  "NUT_FE", "NUT_K", "NUT_FOL", "NUT_MG", "NUT_VITC", "NUT_VITD",
  "VEG", "LEGUMES", "FRUIT", "NUTS", "FISH", "WHGRAIN", "REDPRMEAT"
)

diet <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb22861.tab.gz", 
              data.table=FALSE, stringsAsFactors=FALSE) %>%
  mutate(TCALS = fetch_diet_fields("100002", .),
         CHO = fetch_diet_fields("100005", .),
         SUGARS = fetch_diet_fields("100008", .),
         FAT = fetch_diet_fields("100004", .),
         SFA = fetch_diet_fields("100006", .),
         PUFA = fetch_diet_fields("100007", .),
         PRO = fetch_diet_fields("100003", .),
         ALC = fetch_diet_fields("100022", .),
         FIBER = fetch_diet_fields("100009", .),
         NUT_FE = fetch_diet_fields("100011", .),
         NUT_K = fetch_diet_fields("100016", .),
         NUT_FOL = fetch_diet_fields("100014", .),
         NUT_MG = fetch_diet_fields("100017", .),
         NUT_VITC = fetch_diet_fields("100015", .),
         NUT_VITD = fetch_diet_fields("100021", .),
         VEG = fetch_diet_fields(fields_list$VEG, ., coding=TRUE),
         LEGUMES = fetch_diet_fields(fields_list$LEGUMES, ., coding=TRUE),
         FRUIT = fetch_diet_fields(fields_list$FRUIT, ., coding=TRUE),
         NUTS = fetch_diet_fields(fields_list$NUTS, ., coding=TRUE),
         FISH = fetch_diet_fields(fields_list$FISH, ., coding=TRUE),
         WHGRAIN = fetch_diet_fields(fields_list$WHGRAIN, ., coding=TRUE),
         REDPRMEAT = fetch_diet_fields(fields_list$REDPRMEAT, ., coding=TRUE),
         num_recalls = check_num_valid_24hr(.)) %>%
  mutate(MUFA = FAT - SFA - PUFA,
         TCALS = TCALS / 4.18) %>%  # Energy from kJ to kcals
  mutate_at(vars(CHO, PRO), ~. * 4) %>%  # Nutrients from g to kcals (other than alcohol)
  mutate_at(vars(FAT, SFA, MUFA, PUFA), ~. * 9) %>%
  mutate_at(vars(all_of(diet_vars)), winsorize) %>%
  select(id=f.eid, all_of(diet_vars), num_recalls)

diet_extra <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb22861.tab.gz", 
                    data.table=FALSE, stringsAsFactors=FALSE) %>%
  mutate(
    # N3FA = fetch_diet_fields("26015", .),
    # N6FA = fetch_diet_fields("26016", .),
    baked_beans = fetch_diet_fields("104000", .),
    pulses = fetch_diet_fields("104010", .)
  ) %>%
  # mutate(across(c(N3FA, N6FA), ~. * 9)) %>%
  select(id=f.eid, baked_beans, pulses) %>%
  mutate(across(-id, winsorize))

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

ffq_extra <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb10528.tab.gz", 
                     data.table=FALSE, stringsAsFactors=FALSE) %>%
  select(id = f.eid,
         oily_fish = f.1329.0.0,
         nonoily_fish = f.1339.0.0,
         prmeat = f.1349.0.0,
         poultry = f.1359.0.0,
         beef = f.1369.0.0,
         lamb = f.1379.0.0) %>%
  mutate(across(-id, ffq_cat_to_qt))

# supp_extra <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb38040.tab.gz", 
#                     data.table=FALSE, stringsAsFactors=FALSE) %>%
#   select(id = f.eid,
#          contains("f.6179.0")) %>%
#   mutate(fish_oil = as.integer(rowSums(across(-1, ~. == 1)) >= 1)) %>%
#   select(id, fish_oil)

supp_extra <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb22861.tab.gz", 
                    data.table=FALSE, stringsAsFactors=FALSE) %>%
  select(id = f.eid,
         contains("f.20084.0")) %>%
  mutate(fish_oil = as.integer(rowSums(across(-1, ~. == 472), na.rm=TRUE) >= 1)) %>%
  select(id, fish_oil)

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

diet_full <- diet %>%
  left_join(diet_extra, by="id") %>%
  left_join(ffq_extra, by="id") %>%
  left_join(supp_extra, by="id") %>%
  mutate(MUFA2SFA = MUFA / SFA) %>%
  mutate(mds = ifelse(num_recalls > 0, calc_med_score(.), NA))

### Add genetic PCs and relatedness --------------------------------------------

gPC_df <- fread("/humgen/florezlab/UKBB_app27892/UKBB_app27892_download_BEFORE_aug_2022/ukb10528.tab.gz",
                data.table=FALSE, stringsAsFactors=FALSE) %>%
  select(id=f.eid, contains("f.22009.0"), used_in_PCA=f.22020.0.0) %>%
  rename_with(.fn=~gsub("f.22009.0.", "gPC", .)) %>%
  mutate(unrelated = (used_in_PCA == 1)) %>%  # An unrelated subset was used in the central PCA
  select(id, contains("gPC"), unrelated)

### Merge and write "raw" phenotypes -------------------------------------------

phenos <- base_phenos %>%
  inner_join(biomarker_phenos, by="id") %>%
  left_join(cbc_phenos, by="id") %>%
  left_join(outcomes, by="id") %>%
  left_join(diet_full, by="id") %>%
  left_join(gPC_df, by="id") %>%
  filter(!(id %in% withdrawn_consent)) %>%
  mutate(id = format(id, scientific=FALSE)) %>%
  mutate(across(contains("gPC"), ~. * mds, .names="mdsBy{.col}"))

write_csv(phenos, "../data/processed/ukb_phenos_raw.csv")

### Phenotype processing and exclusions ----------------------------------------

logged_risk_factors <- c("alt", "tg", "hscrp",
                         "aisi", "sii", "siri")
risk_factors <- c(
  "bmi",
  "sbp", "dbp", "sbp_medsadj", "dbp_medsadj",
  "alt_log", 
  "chol", "ldl", "hdl", "apoB",
  "tg_log",
  "hba1c", "glu",
  "hscrp_log",
  "vitD",
  "chol_statinadj", "ldl_statinadj", "apoB_statinadj",
  "aisi_log", "sii_log", "siri_log"
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

write_csv(processed_phenos, "../data/processed/ukb_phenos.csv")

processed_phenos %>%
  filter(unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_unrelated.csv")

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
  write_csv("../data/processed/ukb_phenos_EUR.csv")

processed_phenos_panUKBB %>%
  filter(ancestry == "EUR", unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_EUR_unrelated.csv")

processed_phenos_panUKBB %>%
  filter(ancestry == "AFR", unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_AFR_unrelated.csv")

processed_phenos_panUKBB %>%
  filter(ancestry == "AMR", unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_AMR_unrelated.csv")

processed_phenos_panUKBB %>%
  filter(ancestry == "CSA", unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_CSA_unrelated.csv")

processed_phenos_panUKBB %>%
  filter(ancestry == "EAS", unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_EAS_unrelated.csv")

processed_phenos_panUKBB %>%
  filter(ancestry == "MID", unrelated == TRUE) %>%
  write_csv("../data/processed/ukb_phenos_MID_unrelated.csv")
