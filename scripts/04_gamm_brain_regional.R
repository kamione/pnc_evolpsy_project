# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)
library(glue)
library(mgcv)
library(broom.mixed)
library(ggseg)
library(patchwork)
library(janitor)
library(gtsummary)
library(officer)

source(here("src", "R", "plot_smoothbyfactor.R"))
source(here("src", "R", "plot_ggseg_brain.R"))
source(here("src", "R", "get_gamm_partialR2.R"))
source(here("src", "R", "convert_df2ft.R"))


# Data I/O ---------------------------------------------------------------------
project_data_path <- "data/raw"
extralong_path <- "ExtraLong/data/datafreeze-2021/TabulatedQC"

preprocessed_dat <- here(
    "data", "processed", "data_cohort-extralong_status-preprocessed_decs-persistent.rds"
) %>% 
    read_rds()

aseg_labels <- here(project_data_path, extralong_path, "Longitudinal", 
                     "aseg_stats_2022-02-24.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    select(matches("Thalamus|Caudate|Putamen|Pallidum|Hippocampus|Amygdala")) %>% 
    janitor::clean_names() %>%
    colnames()

aseg_labels_orginal <- here(project_data_path, extralong_path, "Longitudinal", 
     "aseg_stats_2022-02-24.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    select(matches("Thalamus|Caudate|Putamen|Pallidum|Hippocampus|Amygdala")) %>% 
    colnames()
# replace the label here for ggseg
aseg_labels_orginal[1] <- "Left-Thalamus-Proper"
aseg_labels_orginal[7] <- "Right-Thalamus-Proper"


# Regional Cortical GMV --------------------------------------------------------
dk_labels <- dk$data %>% 
    drop_na(label) %>% 
    filter(!(label %in% c("lh_corpuscallosum", "rh_corpuscallosum"))) %>% 
    pull(label) %>% 
    unique()

fitted_models <- list()
res_tables <- list()
for (label in dk_labels) {
    confounds <- "sex + race3 + timepoint + acq_type"
    # confounds <- "sex + race3 + timepoint + acq_type + euler + n_scans + trauma_T1 + envSES + etiv"
    forumla_1 <- as.formula(glue("{label}_volume ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + s(age_at_scan, by = oisPS, k = 4, fx = TRUE) + {confounds}"))
    fit <- gamm(
        formula = forumla_1,
        random = list(bblid =~ 1),
        data = preprocessed_dat,
        REML = TRUE
    )
    fitted_models[[label]] <- fit
    res_table <- fit$gam %>% 
        tidy(parametric = TRUE) %>% 
        bind_rows(tidy(fit$gam)) %>% 
        select(term, statistic, p.value) %>% 
        mutate(label = label, .before = term)
    res_tables[[label]] = res_table
}
res_tables_group <- res_tables %>% 
    bind_rows() %>% 
    filter(term == "oisPS.L") %>% 
    mutate(p.value.adj = p.adjust(p.value, method = "fdr"))
res_tables_age <- res_tables %>%
    bind_rows() %>% 
    filter(term == "s(age_at_scan)") %>% 
    mutate(p.value.adj = p.adjust(p.value, method = "fdr"))
res_tables_agebygroup <- res_tables %>% 
    bind_rows() %>% 
    filter(term == "s(age_at_scan):oisPSrecurrent PS") %>% 
    mutate(p.value.adj = p.adjust(p.value, method = "fdr"))

res_tables_group %>% 
    select(-2) %>% 
    .convert_df2ft(
        footnote = "The model accounted for sex, race, timepoint, and sequence ID as control variables.",
        filename = here("outputs", "tables", "cgmv_group_effect.docx")
    )

res_tables_age %>% 
    select(-2) %>% 
    .convert_df2ft(
        footnote = "The model accounted for sex, race, timepoint, and sequence ID as control variables.",
        filename = here("outputs", "tables", "cgmv_age_effect.docx")
    )

res_tables_agebygroup %>% 
    select(-2) %>% 
    .convert_df2ft(
        footnote = "The model accounted for sex, race, timepoint, and sequence ID as control variables.",
        filename = here("outputs", "tables", "cgmv_agebygroup_effect.docx")
    )


# plot uncorrected result
vol_a <- res_tables_group %>% 
    plot_ggseg_brain(
        fill = "statistic", 
        title = "A. Group (TD - recurrent PS)"
    )
vol_b <- res_tables_age %>% 
    plot_ggseg_brain(
        fill = "statistic",
        title = "B. Age",
        min = 0,
        max = 250,
        break_int = 50,
        highlight_region = TRUE
    )
vol_c <- res_tables_agebygroup %>% 
    plot_ggseg_brain(
        fill = "statistic", 
        title = "C. Age by Group",
        min = 0,
        max = 9,
        break_int = 2
    )
vol_fig <- vol_a / vol_b / vol_c
vol_fig

ggsave(
    plot = vol_fig, 
    filename = here("outputs", "figs", "persistent_ps_uncorrected_gmv.pdf"), 
    width = 10, 
    height = 6
)

traj_d <- plot_smoothbyfactor(
    fitted_models$lh_bankssts, 
    series = "age_at_scan", 
    by = "oisPS", 
    subjid = "bblid",
    title = "D. Left bankssts",
    title_font_size = 18
)
traj_e <- plot_smoothbyfactor(
    fitted_models$lh_lingual, 
    series = "age_at_scan", 
    by = "oisPS", 
    subjid = "bblid",
    title = "E. Left lingual",
    title_font_size = 18
)
traj_f <- plot_smoothbyfactor(
    fitted_models$rh_superiorfrontal, 
    series = "age_at_scan", 
    by = "oisPS", 
    subjid = "bblid",
    title = "F. Right superior frontal",
    title_font_size = 18
)
traj_g <- plot_smoothbyfactor(
    fitted_models$rh_rostralmiddlefronta, 
    series = "age_at_scan", 
    by = "oisPS", 
    subjid = "bblid",
    title = "G. Right rostral middle frontal",
    title_font_size = 18
)

traj_vol <- (traj_d | traj_e) / (plot_spacer() | plot_spacer()) / (traj_f | traj_g) +
    plot_layout(heights = c(2, 0.1, 2))
traj_vol

ggsave(
    plot = traj_vol, 
    filename = here("outputs", "figs", "persistent_ps_agebypsy.pdf"), 
    width = 12, 
    height = 9
)


 # Subcortical Regional Volume --------------------------------------------------
sgmv_fitted_models <- list()
sgmv_res_tables <- list()
for (label in aseg_labels) {
    confounds <- "sex + race3 + timepoint + acq_type"
    forumla_1 <- as.formula(glue("{label} ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + s(age_at_scan, by = oisPS, k = 4, fx = TRUE) + {confounds}"))
    fit_ <- gamm(
        formula = forumla_1,
        random = list(bblid =~ 1),
        data = preprocessed_dat,
        REML = TRUE
    )
    sgmv_fitted_models[[label]] <- fit_
    sgmv_res_table <- fit_$gam %>% 
        tidy(parametric = TRUE) %>% 
        bind_rows(tidy(fit_$gam)) %>% 
        select(term, statistic, p.value) %>% 
        mutate(label = label, .before = term)
    sgmv_res_tables[[label]] = sgmv_res_table
}
sgmv_res_tables_psy <- sgmv_res_tables %>% 
    bind_rows() %>% 
    filter(term == "oisPS.L") %>% 
    mutate(label = aseg_labels_orginal) %>% 
    mutate(p.value.adj = p.adjust(p.value, method = "fdr"))
sgmv_res_tables_age <- sgmv_res_tables %>%
    bind_rows() %>% 
    filter(term == "s(age_at_scan)") %>% 
    mutate(label = aseg_labels_orginal) %>%
    mutate(p.value.adj = p.adjust(p.value, method = "fdr"))
sgmv_res_tables_agebypsy <- sgmv_res_tables %>% 
    bind_rows() %>% 
    filter(term == "s(age_at_scan):oisPSrecurrent PS") %>% 
    mutate(label = aseg_labels_orginal) %>%
    mutate(p.value.adj = p.adjust(p.value, method = "fdr"))

sgmv_res_tables_psy %>% 
    select(-2) %>% 
    .convert_df2ft(
        footnote = "The model accounted for sex, race, timepoint, and sequence ID as control variables.",
        filename = here("outputs", "tables", "sgmv_group_effect.docx")
    )

sgmv_res_tables_age %>% 
    select(-2) %>% 
    .convert_df2ft(
        footnote = "The model accounted for sex, race, timepoint, and sequence ID as control variables.",
        filename = here("outputs", "tables", "sgmv_age_effect.docx")
    )

sgmv_res_tables_agebypsy %>% 
    select(-2) %>% 
    .convert_df2ft(
        footnote = "The model accounted for sex, race, timepoint, and sequence ID as control variables.",
        filename = here("outputs", "tables", "sgmv_agebygroup_effect.docx")
    )



# plot uncorrected result
sgmv_a <- sgmv_res_tables_psy %>% 
    plot_ggseg_brain(
        fill = "statistic", 
        atlas = "aseg",
        title = "A. Group (TD - recurrent PS)"
    )
sgmv_b <- sgmv_res_tables_age %>% 
    plot_ggseg_brain(
        fill = "statistic",
        title = "B. Age",
        atlas = "aseg",
        min = 0,
        max = 42,
        break_int = 10
    )
sgmv_c <- sgmv_res_tables_agebypsy %>% 
    plot_ggseg_brain(
        fill = "statistic",
        title = "C. Age by Group",
        atlas = "aseg",
        min = 0,
        max = 6,
        break_int = 2
    )

sgmv_fig <- sgmv_a / sgmv_b / sgmv_c
sgmv_fig

ggsave(
    plot = sgmv_fig, 
    filename = here("outputs", "figs", "persistent_ps_uncorrected_sgmv.pdf"), 
    width = 5, 
    height = 6
)


traj_sgmv <- plot_smoothbyfactor(
    sgmv_fitted_models$right_thalamus, 
    series = "age_at_scan", 
    by = "oisPS", 
    subjid = "bblid",
    title = "D. Right Thalamus Proper",
    title_font_size = 12
)
traj_sgmv
ggsave(
    plot = traj_sgmv, 
    filename = here("outputs", "figs", "sgmv_agebypsy.pdf"), 
    width = 6, 
    height = 5
)

res_tables_age %>% 
    select(-2) %>% 
    .convert_df2ft(
        footnote = "The model accounted for sex, race, timepoint, and sequence ID as control variables.",
        filename = here("outputs", "tables", "age_effect.docx")
    )

res_tables_agebygroup %>% 
    select(-2) %>% 
    .convert_df2ft(
        footnote = "The model accounted for sex, race, timepoint, and sequence ID as control variables.",
        filename = here("outputs", "tables", "agebygroup_effect.docx")
    )





# Sensitive Analyses -----------------------------------------------------------
fitted_sens_models <- list()
res_sens_tables <- list()
for (label in dk_labels) {
    confounds <- "sex + race3 + timepoint + acq_type + euler + n_scans + trauma_T1 + envSES + etiv"
    forumla_1 <- as.formula(glue("{label}_volume ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + s(age_at_scan, by = oisPS, k = 4, fx = TRUE) + {confounds}"))
    fit <- gamm(
        formula = forumla_1,
        random = list(bblid =~ 1),
        data = preprocessed_dat,
        REML = TRUE
    )
    fitted_sens_models[[label]] <- fit
    res_table <- fit$gam %>% 
        tidy(parametric = TRUE) %>% 
        bind_rows(tidy(fit$gam)) %>% 
        select(term, statistic, p.value) %>% 
        mutate(label = label, .before = term)
    res_sens_tables[[label]] = res_table
}
res_sens_tables_group <- res_sens_tables %>% 
    bind_rows() %>% 
    filter(term == "oisPS.L") %>% 
    mutate(p.value.adj = p.adjust(p.value, method = "fdr"))
res_sens_tables_age <- res_sens_tables %>%
    bind_rows() %>% 
    filter(term == "s(age_at_scan)") %>% 
    mutate(p.value.adj = p.adjust(p.value, method = "fdr"))
res_sens_tables_agebygroup <- res_sens_tables %>% 
    bind_rows() %>% 
    filter(term == "s(age_at_scan):oisPSrecurrent PS") %>% 
    mutate(p.value.adj = p.adjust(p.value, method = "fdr"))

res_sens_tables_group %>% 
    select(-2) %>% 
    .convert_df2ft(
        footnote = "The model accounted for sex, race, timepoint, sequence ID, euler number, numbers of scans, traumatic stress load, envSES, and eTIV as control variables.",
        filename = here("outputs", "tables", "sens_cgmv_group_effect.docx")
    )
res_sens_tables_age %>% 
    select(-2) %>% 
    .convert_df2ft(
        footnote = "The model accounted for sex, race, timepoint, sequence ID, euler number, numbers of scans, traumatic stress load, envSES, and eTIV as control variables.",
        filename = here("outputs", "tables", "sens_cgmv_age_effect.docx")
    )
res_sens_tables_agebygroup %>% 
    select(-2) %>% 
    .convert_df2ft(
        footnote = "The model accounted for sex, race, timepoint, sequence ID, euler number, numbers of scans, traumatic stress load, envSES, and eTIV as control variables.",
        filename = here("outputs", "tables", "sens_cgmv_agebygroup_effect.docx")
    )


sgmv_sens_fitted_models <- list()
sgmv_sens_res_tables <- list()
for (label in aseg_labels) {
    confounds <- "sex + race3 + timepoint + acq_type + euler + timepoint + n_scans + ers + etiv"
    forumla_1 <- as.formula(glue("{label} ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + s(age_at_scan, by = oisPS, k = 4, fx = TRUE) + {confounds}"))
    fit_ <- gamm(
        formula = forumla_1,
        random = list(bblid =~ 1),
        data = preprocessed_dat,
        REML = TRUE
    )
    sgmv_sens_fitted_models[[label]] <- fit_
    res_table <- fit_$gam %>% 
        tidy(parametric = TRUE) %>% 
        bind_rows(tidy(fit_$gam)) %>% 
        select(term, statistic, p.value) %>% 
        mutate(label = label, .before = term)
    sgmv_sens_res_tables[[label]] = res_table
}
sgmv_res_tables_psy <- sgmv_sens_res_tables %>% 
    bind_rows() %>% 
    filter(term == "oisPS.L") %>% 
    mutate(label = aseg_labels_orginal) %>% 
    mutate(p.value.adj = p.adjust(p.value, method = "fdr"))
sgmv_res_tables_age <- sgmv_sens_res_tables %>%
    bind_rows() %>% 
    filter(term == "s(age_at_scan)") %>% 
    mutate(label = aseg_labels_orginal) %>%
    mutate(p.value.adj = p.adjust(p.value, method = "fdr"))
sgmv_res_tables_agebypsy <- sgmv_sens_res_tables %>% 
    bind_rows() %>% 
    filter(term == "s(age_at_scan):oisPSrecurrent PS") %>% 
    mutate(label = aseg_labels_orginal) %>%
    mutate(p.value.adj = p.adjust(p.value, method = "fdr"))

sgmv_res_tables_psy %>% 
    select(-2) %>% 
    .convert_df2ft(
        footnote = "The model accounted for sex, race, timepoint, sequence ID, euler number, numbers of scans, traumatic stress load, envSES, and eTIV as control variables.",
        filename = here("outputs", "tables", "sens_sgmv_group_effect.docx")
    )
sgmv_res_tables_age %>% 
    select(-2) %>% 
    .convert_df2ft(
        footnote = "The model accounted for sex, race, timepoint, sequence ID, euler number, numbers of scans, traumatic stress load, envSES, and eTIV as control variables.",
        filename = here("outputs", "tables", "sens_sgmv_age_effect.docx")
    )
sgmv_res_tables_agebypsy %>% 
    select(-2) %>% 
    .convert_df2ft(
        footnote = "The model accounted for sex, race, timepoint, sequence ID, euler number, numbers of scans, traumatic stress load, envSES, and eTIV as control variables.",
        filename = here("outputs", "tables", "sens_sgmv_agebygroup_effect.docx")
    )

# Sensitive Analyses 2 ---------------------------------------------------------
fitted_sens_models <- list()
res_sens_tables <- list()
for (label in dk_labels) {
    confounds <- "sex + race3 + timepoint + acq_type"
    forumla_1 <- as.formula(glue("{label}_volume ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + s(age_at_scan, by = oisPS, k = 4, fx = TRUE) + {confounds}"))
    fit <- gamm(
        formula = forumla_1,
        random = list(bblid =~ 1),
        data = preprocessed_dat %>% filter(!age_at_scan < 12.5),
        REML = TRUE
    )
    fitted_sens_models[[label]] <- fit
    res_table <- fit$gam %>% 
        tidy(parametric = TRUE) %>% 
        bind_rows(tidy(fit$gam)) %>% 
        select(term, statistic, p.value) %>% 
        mutate(label = label, .before = term)
    res_sens_tables[[label]] = res_table
}
res_sens_tables_group <- res_sens_tables %>% 
    bind_rows() %>% 
    filter(term == "oisPS.L") %>% 
    mutate(p.value.adj = p.adjust(p.value, method = "fdr"))
res_sens_tables_age <- res_sens_tables %>%
    bind_rows() %>% 
    filter(term == "s(age_at_scan)") %>% 
    mutate(p.value.adj = p.adjust(p.value, method = "fdr"))
res_sens_tables_agebygroup <- res_sens_tables %>% 
    bind_rows() %>% 
    filter(term == "s(age_at_scan):oisPSrecurrent PS") %>% 
    mutate(p.value.adj = p.adjust(p.value, method = "fdr"))

res_sens_tables_group %>% 
    select(-2) %>% 
    .convert_df2ft(
        footnote = "The model accounted for sex, race, timepoint, sequence ID as control variables.",
        filename = here("outputs", "tables", "sens2_cgmv_group_effect.docx")
    )
res_sens_tables_age %>% 
    select(-2) %>% 
    .convert_df2ft(
        footnote = "The model accounted for sex, race, timepoint, sequence ID as control variables.",
        filename = here("outputs", "tables", "sens2_cgmv_age_effect.docx")
    )
res_sens_tables_agebygroup %>% 
    select(-2) %>% 
    .convert_df2ft(
        footnote = "The model accounted for sex, race, timepoint, sequence ID, as control variables.",
        filename = here("outputs", "tables", "sens2_cgmv_agebygroup_effect.docx")
    )


sgmv_sens_fitted_models <- list()
sgmv_sens_res_tables <- list()
for (label in aseg_labels) {
    confounds <- "sex + race3 + timepoint + acq_type"
    forumla_1 <- as.formula(glue("{label} ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + s(age_at_scan, by = oisPS, k = 4, fx = TRUE) + {confounds}"))
    fit_ <- gamm(
        formula = forumla_1,
        random = list(bblid =~ 1),
        data = preprocessed_dat %>% filter(!age_at_scan < 12.5),
        REML = TRUE
    )
    sgmv_sens_fitted_models[[label]] <- fit_
    res_table <- fit_$gam %>% 
        tidy(parametric = TRUE) %>% 
        bind_rows(tidy(fit_$gam)) %>% 
        select(term, statistic, p.value) %>% 
        mutate(label = label, .before = term)
    sgmv_sens_res_tables[[label]] = res_table
}
sgmv_res_tables_psy <- sgmv_sens_res_tables %>% 
    bind_rows() %>% 
    filter(term == "oisPS.L") %>% 
    mutate(label = aseg_labels_orginal) %>% 
    mutate(p.value.adj = p.adjust(p.value, method = "fdr"))
sgmv_res_tables_age <- sgmv_sens_res_tables %>%
    bind_rows() %>% 
    filter(term == "s(age_at_scan)") %>% 
    mutate(label = aseg_labels_orginal) %>%
    mutate(p.value.adj = p.adjust(p.value, method = "fdr"))
sgmv_res_tables_agebypsy <- sgmv_sens_res_tables %>% 
    bind_rows() %>% 
    filter(term == "s(age_at_scan):oisPSrecurrent PS") %>% 
    mutate(label = aseg_labels_orginal) %>%
    mutate(p.value.adj = p.adjust(p.value, method = "fdr"))

sgmv_res_tables_psy %>% 
    select(-2) %>% 
    .convert_df2ft(
        footnote = "The model accounted for sex, race, timepoint, sequence ID as control variables.",
        filename = here("outputs", "tables", "sens2_sgmv_group_effect.docx")
    )
sgmv_res_tables_age %>% 
    select(-2) %>% 
    .convert_df2ft(
        footnote = "The model accounted for sex, race, timepoint, sequence ID as control variables.",
        filename = here("outputs", "tables", "sens2_sgmv_age_effect.docx")
    )
sgmv_res_tables_agebypsy %>% 
    select(-2) %>% 
    .convert_df2ft(
        footnote = "The model accounted for sex, race, timepoint, sequence ID as control variables.",
        filename = here("outputs", "tables", "sens2_sgmv_agebygroup_effect.docx")
    )
