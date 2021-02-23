library(dplyr)
library(lme4)
library(glmmTMB)

long0 <- readRDS('thbp_7yr_long_20200617.rds')  # cognitive test score transformations detailed below

# Table 4
m1.un <- lmer(Score ~ Time*group + (1+Time|Participant), long0)
m1.adj <- lmer(Score ~ scale(age_1) + `Prior CR`*Time + Time*group + (1+Time|Participant), long0)
labs <- c("(Intercept)",
          "Time (years)",
          "Educational intervention",
          "Time x Intervention",
          "Age (years)",
          "Prior CR (z)",
          "Prior CR x Time")
tab_model(m1.un, m1.adj,
          digits = 3,
          pred.labels = labs)

# Table S1
m2.un <- glmmTMB::glmmTMB(value ~ test*Time*group + (1+test|Participant), long0)
m2.adj <- glmmTMB::glmmTMB(value ~ `Prior CR` + scale(age_1)*test + Time*test*group + (1+test|Participant), long0)

# Table S2
library(emmeans)
emtrends(m2.adj, pairwise ~ group | test, var = "Time", show.levels = TRUE)

# Figure 1
p1 <- expand.grid(`Prior CR` = 0,
                  age_1 = 60.3,
                  group = c("Control", "Intervention"),
                  Time = c(0,1,2,3,5,7),
                  test = unique(m2.adj$frame$test),
                  Participant = NA)
p1$age_1 <- p1$age_1+p1$Time

fit <- predict(m2.adj, p1, se.fit = TRUE)


p1$fit <- fit$fit
p1$lwr <- fit$fit - 1.96*fit$se.fit
p1$upr <- fit$fit + 1.96*fit$se.fit

test_names <- c(
  `ravlt_t15_tot_z` = "RAVLT tot",
  `ravlt_arcl_raw_z` = "RAVLT rcl",
  `rcft_rcl_raw_z` = "RCFT",
  `lmi_unit_tot_z` = "LMi",
  `lmii_unit_tot_z` = "LMii",
  `reversed_log_pal_te6_score_z` = "PAL te6",
  `cowat_tot_z` = "COWAT",
  `wais_ds_raw_z` = "WAIS ds",
  `reversed_log_tmt_b_time_z` = "TMT-B",
  `pal_ftm_score_z` = "PAL ftm",
  `ssp_len_z` = "SSP length",
  `reversed_log_swm_be_z` = "SWM be",
  `reversed_log_stroop_c_time_z` = "STROOP C time",
  `rvp_a_z` = "RVP-A",
  `reversed_log_bnt_raw_z` = "BNT",
  `wais_com_raw_z` = "WAIS comp",
  `wais_lns_raw_z` = "WAIS lns",
  `reversed_log_wais_voc_raw_z` = "WAIS voc",
  `-1.1` = "55 years of age",
  `0.3` = "65 years of age",
  `1.7` = "75 years of age"
)

library(ggplot2)
library(latex2exp)

ggplot(p1, aes(x = Time, y = fit, ymin = lwr, ymax = upr, colour = group, group = group)) +
  facet_wrap(~test, labeller = as_labeller(test_names)) +
  geom_line() +
  geom_ribbon(aes(fill = group, colour = NULL), alpha = 0.17) +
  geom_point(size = 1.5) +
  theme_bw() +
  scale_color_manual(values = c("slategrey", "orange3")) +
  scale_fill_manual(values = c("slategrey", "orange3")) +
  ylab(TeX("$\\widehat{Score\\,(z)}$")) +
  xlab("Time (years since baseline)") +
  guides(colour = guide_legend(title = NULL),
         fill = guide_legend(title = NULL)) +
  theme(legend.position = c(1, 0),
        legend.justification = c(0.93, 0)) 
ggsave("fig1.pdf")


# Table S3
m3.adj <- glmmTMB::glmmTMB(Score ~ `Prior CR` + Time*scale(age_1)*group + (1+Time|Participant), long0)

# Figure 2
p1 <- expand.grid(`Prior CR` = 0,
                  age_at_baseline = c(55,60,65,70,75),
                  Time = c(0,1,2,3,4,5,6,7),
                  group = c("Control", "Intervention"),
                  Participant = NA)
p1$age_1 <- p1$age_at_baseline + p1$Time
p1$age_at_baseline <- paste0(p1$age_at_baseline, " years at baseline")

fit <- predict(m3.adj, p1, se.fit = TRUE)
p1$fit <- fit$fit
p1$lwr <- fit$fit - 1.96*fit$se.fit
p1$upr <- fit$fit + 1.96*fit$se.fit

ggplot(p1 %>% filter(!age_at_baseline %in% c('60 years at baseline', '70 years at baseline')),
               aes(x = Time, y = fit, ymin = lwr, ymax = upr,
                   colour = group, group = group)) +
  facet_wrap(~age_at_baseline, ncol = 4) +
  geom_ribbon(aes(fill = group, colour = NULL), alpha = 0.2) +
  geom_line(size = 1) +
  theme_bw() +
  ylab(TeX("$\\widehat{Score\\,(z)}$")) +
  xlab("Time (years since baseline)") +
  scale_color_manual(values = c("slategrey", "orange3"), name = "Group") +
  scale_fill_manual(values = c("slategrey", "orange3"), name = "Group") +
  theme(legend.position = "bottom")
ggsave("fig2.pdf")

# Table S4
m4.un <- glmmTMB::glmmTMB(value ~ Time*scale(`Academic Load`) + (1+Time|Participant), long0)
m4.adj <- glmmTMB::glmmTMB(Score ~ `Prior CR` + scale(age_1) + Time*test*scale(`Academic Load`) + (1|Participant), long0)

# Table S5
m5.un <- lmer(value ~ group*Time*`Prior CR` + (1|Participant), long0)
m5.adj <- lmer(value ~ age_z*Time + Time*group*`Prior CR` + (1|Participant), long0)


# Cognitive test score transformations which were applied
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
         reversed_log_bnt_raw_z = ifelse(reversed_log_bnt_raw_z < -6, NA, reversed_log_bnt_raw_z))