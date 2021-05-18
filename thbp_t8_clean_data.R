library(readr)
library(tidyr)
library(dplyr)

# read raw data from thbp database, convert names to lower case, make some transformations
d <- read_csv("Healthy Brain Assessment Export 2May19.csv", trim_ws = TRUE, guess_max = 3000) %>%
  filter(phase %in% c(1,2,3,4,6)) %>%
  bind_rows(read_csv("Healthy Brain Assessment Export T8_9thOct2019.csv", trim_ws = TRUE, guess_max = 3000)) %>%
  mutate(phase_f = factor(phase))
names(d) <- tolower(names(d))
d <- d[, !duplicated(colnames(d))]


# catch errors in gen_apoe_carrier and gen_bdnf_carrier
carrier_apoe <- function(a, b){
  ifelse(a == "e4", TRUE, ifelse(b == "e4", TRUE, FALSE))
}

het_apoe <- function(carrier, a, b){
  ifelse(carrier, 
         ifelse(a == 'e3', TRUE,
                ifelse(b == 'e3', TRUE, FALSE)),
         FALSE)
}

gen_apoe <- dplyr::select(d, idcode, starts_with("gen_apoe")) %>%
  filter(!is.na(gen_apoe_geno1), !is.na(gen_apoe_geno2)) %>%
  mutate(carrier = carrier_apoe(gen_apoe_geno1, gen_apoe_geno2)) %>%
  mutate(gen_apoe_carrier = carrier,
         gen_apoe_het = het_apoe(gen_apoe_carrier, gen_apoe_geno1, gen_apoe_geno2)) %>%
  dplyr::select(idcode, gen_apoe_carrier, gen_apoe_het)  %>%
  filter(!duplicated(idcode))

#xtabs(~gen_apoe_carrier + carrier, gen_apoe)   # 5 carriers not classified correctly

carrier_bdnf <- function(a, b){
  ifelse(a == "met", TRUE, ifelse(b == "met", TRUE, FALSE))
}

gen_bdnf <- dplyr::select(d, idcode, starts_with("gen_bdnf")) %>%
  filter(!is.na(gen_bdnf_geno1), !is.na(gen_bdnf_geno2)) %>%
  mutate(carrier = carrier_bdnf(gen_bdnf_geno1, gen_bdnf_geno2)) %>%
  mutate(gen_bdnf_carrier = carrier) %>%
  dplyr::select(idcode, gen_bdnf_carrier)  %>%
  filter(!duplicated(idcode))

#xtabs(~gen_bdnf_carrier + carrier, gen_bdnf)   # 8 carriers not classified correctly

d <- dplyr::select(d, -c(gen_apoe_carrier, gen_bdnf_carrier)) %>%
  left_join(gen_apoe) %>%
  left_join(gen_bdnf) %>%
  arrange(idcode, phase_f)

d0 <- dplyr::filter(d,
                    phase_f %in% c(1,2,3,4,6,8),
                    age_1 >= 50,
                    !idcode %in% c('THBS0259', 'THBS0137', ' THBS0224', 'THBS0188')) %>%    # selected only phase 1,2,3,4,6
  mutate(tmt_b_time = ifelse(tmt_b_time <= 5, NA, tmt_b_time),         
         tmt_b_time = ifelse(tmt_b_time >= 300, NA, tmt_b_time)) %>%  # outliers discussed with KS
  mutate(reversed_log_tmt_b_time = -1*log10(tmt_b_time),
         reversed_log_stroop_c_time = -1*log10(stroop_c_time),
         age_z = scale(age_1)[,1],
         ravlt_t15_tot_z = scale(ravlt_t15_tot),
         ravlt_arcl_raw_z = scale(ravlt_arcl_raw),
         rcft_rcl_raw_z = scale(rcft_rcl_raw),
         lmi_unit_tot_z = scale(lmi_unit_tot),
         lmii_unit_tot_z = scale(lmii_unit_tot),
         cowat_tot_z = scale(cowat_tot),
         wais_ds_raw_z = scale(wais_ds_raw),
         reversed_log_tmt_b_time_z = scale(reversed_log_tmt_b_time),
         pal_ftm_score_z = scale(pal_ftm_score),
         reversed_log_pal_te6_score_z = scale(log(51/(pal_te6_score+1))),
         wais_lns_raw_z = scale(wais_lns_raw),
         ssp_len_z = scale(ssp_len),
         reversed_log_swm_be_z = scale(-log10(swm_be+1)),
         reversed_log_stroop_c_time_z = scale(reversed_log_stroop_c_time),
         rvp_a_z = scale(rvp_a),
         reversed_log_bnt_raw_z = -scale(log10(60-bnt_raw+1)),
         wais_com_raw_z = scale(wais_com_raw),
         reversed_log_wais_voc_raw_z = -scale(log10(66-wais_voc_raw+1))) %>%
  mutate(wais_com_raw_z = ifelse(wais_com_raw_z > 5, NA, wais_com_raw_z),
         wais_com_raw_z = ifelse(wais_com_raw_z < -6, NA, wais_com_raw_z),
         reversed_log_bnt_raw_z = ifelse(reversed_log_bnt_raw_z < -6, NA, reversed_log_bnt_raw_z)) %>%
  arrange(idcode, phase_f) %>%
  dplyr::select(idcode,
                group,
                phase_f,
                phase,
                ravlt_t15_tot_z,
                ravlt_arcl_raw_z,
                rcft_rcl_raw_z,
                lmi_unit_tot_z,
                lmii_unit_tot_z,
                cowat_tot_z,
                wais_ds_raw_z,
                reversed_log_tmt_b_time_z, 
                pal_ftm_score_z,
                reversed_log_pal_te6_score_z,
                wais_lns_raw_z,
                ssp_len_z,
                reversed_log_swm_be_z,
                reversed_log_stroop_c_time_z,
                rvp_a_z,
                reversed_log_bnt_raw_z,
                wais_com_raw_z,
                reversed_log_wais_voc_raw_z,
                age_1,
                age_z)

# academic load
load <- read_csv("thbp_load.csv", trim_ws = TRUE, guess_max = 3000) %>%
  mutate(idcode = paste0("THBS", stringr::str_pad(IDCode, 4, pad = "0"))) %>%  # resolve integers to THBS codes, check with AR and KS
  arrange(idcode) %>%
  dplyr::select(-IDCode)

# for PCA without wtar_fsiq_z
pcr0 <- dplyr::filter(d, phase_f == 1) %>%
  mutate(
    wtar_fsiq_z = c(scale(as.numeric(wtar_fsiq))),
    mhq_edschool = c(scale(as.numeric(mhq_edschool))),
    leq_ya_spec = c(scale(as.numeric(leq_ya_spec))),
    leq_ml_nonspec = c(scale(as.numeric(leq_ml_nonspec))),
    leq_ya_nonspec = c(scale(as.numeric(leq_ya_nonspec))),
    leq_ml_bonus = c(scale(as.numeric(leq_ml_bonus))),
    leq_ml_spec = c(scale(as.numeric(leq_ml_spec)))
  ) %>%
  dplyr::select(idcode, wtar_fsiq_z, mhq_edschool, leq_ya_spec, leq_ml_nonspec, leq_ya_nonspec, leq_ml_bonus, leq_ml_spec) %>%
  na.omit()

pca <- prcomp(pcr0[, 3:8])
pcr0$pca <- scale(rowSums(pca$x[, 1:2]))[,1]  # pc1 & pc2 explain ~70% of variance

pca <- prcomp(pcr0[, 2:8])
pcr0$pca_wtar <- scale(rowSums(pca$x[, 1:2]))[,1]


apoe_geno <- function(carrier, het){
  ifelse(carrier,
         ifelse(het, "e3/e4", "e4/e4"),
         "e4-")
}


# prior cognitive reserve, academic load and other baseline variables
pcr <- dplyr::filter(d, phase_f == 1) %>%
  mutate(
    wtar_fsiq_z = scale(as.numeric(wtar_fsiq)),
    leq_ya_spec = scale(as.numeric(leq_ya_spec)),
    leq_ml_nonspec = scale(as.numeric(leq_ml_nonspec)),
    leq_ya_nonspec = scale(as.numeric(leq_ya_nonspec)),
    leq_ml_bonus = scale(as.numeric(leq_ml_bonus)),
    leq_ml_spec = scale(as.numeric(leq_ml_spec))
  ) %>%
  mutate(
    pcr = .370 * wtar_fsiq_z +
      .408 * mhq_edschool +
      .567 * leq_ya_spec +
      .565 * leq_ya_nonspec +
      .630 * leq_ml_nonspec +
      .875 * leq_ml_bonus +
      1.004 * leq_ml_spec
  ) %>%
  mutate(pcr = scale(pcr)[,1]) %>%
  dplyr::select(
    idcode,
    pcr,
    hads_anx,
    hads_dep,
    wtar_fsiq,
    wtar_fsiq_z,
    starts_with("mhq_"),
    gender_1,
    age_1,
    starts_with("gen_")
  ) %>%
  mutate(mhq_weight = ifelse(mhq_weight > 150, NA, mhq_weight),
         mhq_weight = ifelse(mhq_weight < 40, NA, mhq_weight),
         mhq_height = ifelse(mhq_height >= 200, mhq_height/10, mhq_height),
         mhq_height = ifelse(mhq_height < 100, NA, mhq_height)
  ) %>%
  transmute(
    idcode = idcode,
    pcr = pcr,
    mhq_edschool = mhq_edschool,
    bmi_baseline = mhq_weight/(mhq_height/100)^2,
    mhq_schoolage = mhq_schoolage,
    mhq_highblood = mhq_highblood,
    mhq_lowblood = mhq_lowblood,
    hads_anx_at_baseline_z = scale(hads_anx),
    hads_dep_at_baseline_z = scale(hads_dep),
    age_at_baseline = age_1,
    age_at_baseline_z = scale(age_1),
    wtar_fsiq = as.numeric(wtar_fsiq),
    wtar_fsiq_at_baseline_z = wtar_fsiq_z,
    gender = factor(
      gender_1,
      levels = c("F", "M"),
      labels = c("Female", "Male")
    ),
    apoe = factor(
      gen_apoe_carrier,
      levels = c(FALSE, TRUE),
      labels = c("e4-", "e4+")
    ),
    apoe_geno = factor(apoe_geno(gen_apoe_carrier, gen_apoe_het)),
    bdnf = factor(
      gen_bdnf_carrier,
      levels = c(FALSE, TRUE),
      labels = c("met-", "met+")
    )
  ) %>%
  left_join(load, by = "idcode") %>%  # Ss who aren't in load will have NA load
  dplyr::mutate(current_course_cp = if_else(is.na(current_course_cp), 0, current_course_cp),   
                all_courses_post2010_cp = if_else(is.na(all_courses_post2010_cp), 0, all_courses_post2010_cp)) %>%
  mutate(load_z = scale(log(all_courses_post2010_cp/100+1)))  ## change NA to 0 and calculate transformed load z-scores

# join baseline measures and pcr to d0
d0 <- left_join(d0, pcr, by = "idcode") %>%
  arrange(idcode, phase_f) %>%
  left_join(dplyr::select(pcr0, idcode, pca, pca_wtar))

# linear interpolation of load
d0 <- group_by(d0, idcode) %>% 
  arrange(idcode, phase_f) %>%
  mutate(max_phase = max(phase, na.rm = TRUE)) %>%
  mutate(load_interp = (all_courses_post2010_cp/(max_phase-1)) * (phase-1),
         load_interp = ifelse(is.nan(load_interp), 0, load_interp),
         load_interp_log = log(load_interp+1)) %>%
  ungroup()


# assign human-readable labels to tables
label(d0$age_1) <- "Age"
label(d0$gender) <- "Gender"
label(d0$pcr) <- "Prior Cognitive Reserve (z)"
label(d0$apoe) <- "APOE e4 carrier"
label(d0$phase_f) <- "Phase"
label(d0$bdnf) <- "BDNF Val66Met"
label(d0$current_course_cp) <- "THBP academic course load"
label(d0$all_courses_post2010_cp) <- "All academic course load since 2010" 
label(d0$load_z) <- "log(All academic course load since 2010) z-score"
label(d0$load_interp) <- "Linear interpolation of all academic course load since 2010"
label(d0$load_interp_log) <- "log(load_interp + 1)"
label(d0$mhq_edschool) <- "Education (years)"
label(d0$wtar_fsiq) <- "WTAR FSIQ"
label(d0$group) <- "Educational Intervention Group"
# make a baseline dataset with preserved labels (dplyr loses the labels)
d0_baseline <- dplyr::filter(d0, phase_f == 1) %>% sjlabelled::copy_labels(d0)

saveRDS(d0, "thbp_t8_clean_20200602.rds")

# all variables to be preserved during melt to long/case form
idvars = c("idcode", "phase_f", "phase", "group", "age_z", "age_1", "pcr", 
           "pca", "pca_wtar", "wtar_fsiq",
           "hads_anx_at_baseline_z", "hads_dep_at_baseline_z",
           "age_at_baseline", "age_at_baseline_z",
           "wtar_fsiq_at_baseline_z", "gender",
           "apoe", "bdnf", "apoe_geno",
           "mhq_edschool",
           "mhq_schoolage",
           "mhq_highblood",
           "mhq_lowblood",
           "bmi_baseline",
           "course_code",
           "gpa",
           "comment",
           "current_course_cp",
           "all_courses_post2010_cp",
           "hdr",
           "load_z",
           "load_interp",
           "load_interp_log",
           "max_phase")


long <- reshape2::melt(d0, id.var = c(idvars)) %>%
  mutate(test = factor(variable)) %>%
  dplyr::select(-variable) %>%
  arrange(idcode, phase_f, test)

saveRDS(long, file = "thbp_t8_long_20200602.rds")
saveRDS(d0_baseline, file = 'thbp_t8_baseline_20200602.rds')
