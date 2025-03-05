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
aseg_labels_orginal[1] <- "Left-Thalamus-Proper"
aseg_labels_orginal[7] <- "Right-Thalamus-Proper"


# Global GMV -------------------------------------------------------------------
## cortical
confounds <- "sex + race3 + timepoint + acq_type"

forumla_1_full <- as.formula(
    glue("gmv ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + s(age_at_scan, by = oisPS, k = 4, fx = TRUE) + {confounds}")
)
fit_1_full <- gamm(
    formula = forumla_1_full,
    random = list(bblid =~ 1),
    data = preprocessed_dat,
    REML = TRUE
)
parameters::model_parameters(fit_1_full)

forumla_1_reduced <- as.formula(
    glue("gmv ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + {confounds}")
)
fit_1_reduced <- gamm(
    formula = forumla_1_reduced,
    random = list(bblid =~ 1),
    data = preprocessed_dat,
    REML = TRUE
)

get_gamm_partialR2(fit_1_full, fit_1_reduced)
options(scipen = 999)
gmv_traj <- plot_smoothbyfactor(
    modobj = fit_1_full, 
    series = "age_at_scan", 
    by = "oisPS", 
    subjid = "bblid",
    title = "A. Cortical gray matter volume",
    title_font_size = 18
)

## sgmv
forumla_2_full <- as.formula(
    glue("sgmv ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + s(age_at_scan, k = 4, by = oisPS, fx = TRUE) + {confounds}")
)
fit_2_full <- gamm(
    formula = forumla_2_full,
    random = list(bblid =~ 1),
    data = preprocessed_dat,
    REML = TRUE
)
parameters::model_parameters(fit_2_full)

forumla_2_reduced <- as.formula(
    glue("sgmv ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + {confounds}")
)
fit_2_reduced <- gamm(
    formula = forumla_2_reduced,
    random = list(bblid =~ 1),
    data = preprocessed_dat,
    REML = TRUE
)
get_gamm_partialR2(fit_2_full, fit_2_reduced)

sgmv_traj <- plot_smoothbyfactor(
    modobj = fit_2_full, 
    series = "age_at_scan", 
    by = "oisPS", 
    subjid = "bblid",
    title = "B. Subcortical gray matter volume",
    title_font_size = 18
)

gmv_traj | sgmv_traj
ggsave(
    plot = gmv_traj | sgmv_traj, 
    filename = here("outputs", "figs", "figure01.pdf"), 
    width = 12, 
    height = 5
)


# Sensitivity Analysis ---------------------------------------------------------
confounds_sens <- "sex + race3 + timepoint + acq_type + n_scans + euler + trauma_T1 + envSES"

cgmv_full_sens_formula <- as.formula(
    glue("gmv ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + s(age_at_scan, by = oisPS, k = 4, fx = TRUE) + {confounds_sens}")
)
cgmv_full_sens_fit <- gamm(
    formula = cgmv_full_sens_formula,
    random = list(bblid =~ 1),
    data = preprocessed_dat,
    REML = TRUE
)
parameters::model_parameters(cgmv_full_sens_fit)

cgmv_reduced_sens_formula <- as.formula(
    glue("gmv ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + {confounds_sens}")
)
cgmv_reduced_sens_fit <- gamm(
    formula = cgmv_reduced_sens_formula,
    random = list(bblid =~ 1),
    data = preprocessed_dat,
    REML = TRUE
)
get_gamm_partialR2(cgmv_full_sens_fit, cgmv_reduced_sens_fit)

cgmv_sens_traj <- plot_smoothbyfactor(
    modobj = cgmv_full_sens_fit, 
    series = "age_at_scan", 
    by = "oisPS", 
    subjid = "bblid",
    title = "A. Cortical gray matter volume",
    title_font_size = 18
)

## sgmv
sgmv_full_sens_formula <- as.formula(
    glue("sgmv ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + s(age_at_scan, k = 4, by = oisPS, fx = TRUE) + {confounds_sens}")
)
sgmv_full_sens_fit <- gamm(
    formula = sgmv_full_sens_formula,
    random = list(bblid =~ 1),
    data = preprocessed_dat,
    REML = TRUE
)
parameters::model_parameters(sgmv_full_sens_fit)

sgmv_reduced_sens_formula <- as.formula(
    glue("sgmv ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + {confounds_sens}")
)
sgmv_reduced_sens_fit <- gamm(
    formula = sgmv_reduced_sens_formula,
    random = list(bblid =~ 1),
    data = preprocessed_dat,
    REML = TRUE
)
get_gamm_partialR2(sgmv_full_sens_fit, sgmv_reduced_sens_fit)

sgmv_sens_traj <- plot_smoothbyfactor(
    modobj = sgmv_full_sens_fit, 
    series = "age_at_scan", 
    by = "oisPS", 
    subjid = "bblid",
    title = "B. Subcortical gray matter volume",
    title_font_size = 18
)

parameters::model_parameters(sgmv_full_sens_fit)


# Sensitivity Analysis 2 -------------------------------------------------------
confounds_sens <- "sex + race3 + timepoint + acq_type"

cgmv_full_sens_formula <- as.formula(
    glue("gmv ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + s(age_at_scan, by = oisPS, k = 4, fx = TRUE) + {confounds_sens}")
)
cgmv_full_sens_fit <- gamm(
    formula = cgmv_full_sens_formula,
    random = list(bblid =~ 1),
    data = preprocessed_dat %>% filter(!age_at_scan < 12.5),
    REML = TRUE
)
parameters::model_parameters(cgmv_full_sens_fit)

cgmv_reduced_sens_formula <- as.formula(
    glue("gmv ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + {confounds_sens}")
)
cgmv_reduced_sens_fit <- gamm(
    formula = cgmv_reduced_sens_formula,
    random = list(bblid =~ 1),
    data = preprocessed_dat %>% filter(!age_at_scan < 12.5),
    REML = TRUE
)
get_gamm_partialR2(cgmv_full_sens_fit, cgmv_reduced_sens_fit)

plot_smoothbyfactor(
    modobj = cgmv_full_sens_fit, 
    series = "age_at_scan", 
    by = "oisPS", 
    subjid = "bblid",
    title = "B. Subcortical gray matter volume",
    title_font_size = 18
)

## sgmv
sgmv_full_sens_formula <- as.formula(
    glue("sgmv ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + s(age_at_scan, k = 4, by = oisPS, fx = TRUE) + {confounds_sens}")
)
sgmv_full_sens_fit <- gamm(
    formula = sgmv_full_sens_formula,
    random = list(bblid =~ 1),
    data = preprocessed_dat %>% filter(!age_at_scan < 12.5),
    REML = TRUE
)
parameters::model_parameters(sgmv_full_sens_fit)

sgmv_reduced_sens_formula <- as.formula(
    glue("sgmv ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + {confounds_sens}")
)
sgmv_reduced_sens_fit <- gamm(
    formula = sgmv_reduced_sens_formula,
    random = list(bblid =~ 1),
    data = preprocessed_dat %>% filter(!age_at_scan < 12.5),
    REML = TRUE
)
get_gamm_partialR2(sgmv_full_sens_fit, sgmv_reduced_sens_fit)

plot_smoothbyfactor(
    modobj = sgmv_full_sens_fit, 
    series = "age_at_scan", 
    by = "oisPS", 
    subjid = "bblid",
    title = "B. Subcortical gray matter volume",
    title_font_size = 18
)

parameters::model_parameters(sgmv_full_sens_fit)
