# Environment ------------------------------------------------------------------
library(tidyverse)
library(here)
library(glue)
library(janitor)
library(lubridate)
library(gtsummary)
library(flextable)
library(officer)

source(here("src", "R", "plot_longitudinal_distribution.R"))


# Global Variables -------------------------------------------------------------
project_data_path <- "data/raw"
extralong_path <- "ExtraLong/data/datafreeze-2021/TabulatedQC"


# Data I/O ---------------------------------------------------------------------
fs_lh_volume <- here(project_data_path, extralong_path, "longitudinal",
               "aparc_volume_lh_2022-02-24.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    select(-c("lh.aparc.volume", "BrainSegVolNotVent", "eTIV"))

fs_rh_volume <- here(project_data_path, extralong_path, "longitudinal",
                     "aparc_volume_rh_2022-02-24.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    select(-c("rh.aparc.volume", "BrainSegVolNotVent", "eTIV"))

# using ENIGMA protocol
subcortical_volume <- here(project_data_path, extralong_path, "Longitudinal", 
                           "aseg_stats_2022-02-24.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    select(
        bblid, seslabel, CortexVol, SubCortGrayVol, EstimatedTotalIntraCranialVol,
        matches("Thalamus|Caudate|Putamen|Pallidum|Hippocampus|Amygdala|Accumbens")
    ) %>% 
    rename(
        "gmv" = "CortexVol",
        "sgmv" = "SubCortGrayVol", 
        "etiv" = "EstimatedTotalIntraCranialVol"
    ) %>% 
    janitor::clean_names()

ers <- here("data", "raw", "ers_xavier_gur_2.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    mutate(ers = ERS_Bifactor) %>% 
    select(bblid, ers, gaf001, envSES, trauma_T1, Overall_Efficiency,
           F1_Social_Cognition_Efficiency, F2_Complex_Reasoning_Efficiency,
           F3_Memory_Efficiency, F4_Executive_Efficiency,
           MOOD_CorrTraits, FEAR_CorrTraits, EXTERNALIZING_CorrTraits,
           PSYCHOSIS_CorrTraits)

imginfo <- here(project_data_path, extralong_path, 
               "imglook.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    rename("seslabel" = "SCANID") %>% 
    janitor::clean_names("lower_camel") %>% 
    select(bblid, seslabel, doscan)

demog <- here(project_data_path, extralong_path, "subject.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    janitor::clean_names("lower_camel") %>% 
    select(bblid, dobirth, sex, race)

goassess <- here(project_data_path, "n1601_dataFreeze", "clinical",
                 "n1601_diagnosis_dxpmr_20170509.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    select(bblid, goassessDxpmr7)

protocol <- here(project_data_path, extralong_path, "scan_metadata.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    select(bblid, sesid, ManufacturersModelName, RepetitionTime, EchoTime, 
           InversionTime, SliceThickness, AcquisitionMatrixPE, ReconMatrixPE,
           FlipAngle, PixelBandwidth, ReconMatrixPE, 
           BaseResolution) %>% 
    mutate(
        bblid = as.numeric(bblid), 
        sesid = as.numeric(sesid),
        RepetitionTime = RepetitionTime * 1000,
        EchoTime = round(EchoTime * 1000, 1),
        InversionTime = InversionTime * 1000
    ) %>% 
    unite("acq_parameters", ManufacturersModelName:BaseResolution, remove = FALSE) %>% 
    rename(
        "seslabel" = "sesid",
        "scanner" = "ManufacturersModelName"
    ) %>% 
    select(bblid, seslabel, acq_parameters)

euler <- here(project_data_path, "ExtraLong/data/datafreeze-2021/QC",
              "exclusion_datafreeze-2021_cutoff-212.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    rename(
        "bblid" = "subid",
        "seslabel" = "sesid"
    ) %>% 
    mutate(
        bblid = as.numeric(bblid),
        seslabel = as.numeric(seslabel),
    )

cnb_factor <- here("data", "external", "longitudinal_cnb_factor_scores_extralong.csv") %>% 
    read_csv(show_col_types = FALSE)

cnb_composite <- here("data", "external", "cnb_extralong_composites.csv") %>% 
    read_csv(show_col_types = FALSE)

cnb <- cnb_factor %>% 
    left_join(cnb_composite, by = c("bblid", "timepoint"))


# Preprocessing ----------------------------------------------------------------
preprocessed_dat <- fs_lh_volume %>% 
    left_join(fs_rh_volume, by = c("bblid", "seslabel")) %>% 
    left_join(subcortical_volume, by = c("bblid", "seslabel")) %>% 
    left_join(goassess, by = c("bblid")) %>% 
    filter(goassessDxpmr7 %in% c("TD", "PS")) %>%
    left_join(demog, by = c("bblid")) %>%
    left_join(ers, by = c("bblid")) %>%
    left_join(imginfo, by = c("bblid", "seslabel")) %>%
    left_join(protocol, by = c("bblid", "seslabel")) %>% 
    left_join(euler, by = c("bblid", "seslabel")) %>% 
    drop_na(acq_parameters) %>% 
    mutate(
        dobirth = dmy(dobirth),
        doscan = dmy(doscan),
        age_at_scan = as.numeric(interval(dobirth, doscan), unit = "year"),
        sex = factor(sex, levels = c(1, 2), labels = c("Male", "Female")),
        race3 = factor(case_when(
            race == 1 ~ "White",
            race == 2 ~ "Black",
            TRUE ~ "Others"
        ), levels = c("White", "Black", "Others"))
    ) %>%
    group_by(bblid) %>% 
    mutate(
        n_scans = n(),
        timepoint = row_number()
    ) %>% 
    ungroup() %>% 
    group_by(acq_parameters) %>% 
    mutate(n_acq_parameters = n()) %>% 
    ungroup() %>% 
    select(-c(dobirth, race)) %>% 
    rowwise() %>% 
    ungroup()

included_bblid <- preprocessed_dat %>% 
    pull(bblid) %>% 
    unique()

persistent_subjs_dat <- here(project_data_path, extralong_path,
     "img_data_matched_with_freesurfer_2023_05_01.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    janitor::clean_names("lower_camel") %>% 
    select(bblid, dxPscat, doscan, acq) %>% 
    group_by(bblid) %>% 
    arrange(doscan) %>% 
    drop_na(dxPscat) %>% 
    mutate(
        n_assess = n(),
        n_td = length(which(dxPscat == "noDSMdx")),
        n_ps = length(which(dxPscat == "psy")) + length(which(dxPscat == "pro"))
    ) %>% 
    ungroup() %>% 
    left_join(goassess, by = c("bblid")) %>% 
    filter(bblid %in% included_bblid) %>% 
    filter(n_assess == n_td | n_ps > 0) %>%
    mutate(group = case_when(
        (n_assess == n_td & goassessDxpmr7 == "TD") ~ "TD",
        (n_assess == n_ps & goassessDxpmr7 == "PS") ~ "recurrent PS",
        TRUE ~ NA
    )) %>%
    group_by(bblid) %>% 
    slice(1) %>% 
    ungroup() %>% 
    select(bblid, group, dxPscat)

preprocessed_persistent_dat <- preprocessed_dat %>% 
    left_join(persistent_subjs_dat, by = "bblid") %>% 
    mutate(bblid = factor(bblid)) %>% 
    mutate(
        oisPS = ordered(group, levels = c("TD", "recurrent PS"), labels = c("TD", "recurrent PS"))
    ) %>% 
    mutate(acq_type = factor(acq_parameters, labels = 1:length(unique(.$acq_parameters)))) %>% 
    drop_na(group) %>% 
    filter(exclude == FALSE)

cnb_persistent <- here(project_data_path, extralong_path,
     "img_data_matched_with_freesurfer_2023_05_01.csv") %>% 
    read_csv(show_col_types = FALSE) %>% 
    left_join(cnb, by = c("bblid", "timepoint")) %>% 
    select(bblid, doscan, Executive_CorrTraits:overall_acc) %>% 
    mutate(bblid = factor(bblid))

preprocessed_cnb_persistent_dat <- preprocessed_persistent_dat %>% 
    left_join(cnb_persistent, by = c("bblid", "doscan"))

n_subjs <- preprocessed_cnb_persistent_dat$bblid %>% n_distinct()
n_ps <- preprocessed_cnb_persistent_dat %>% 
    filter(oisPS == "recurrent PS") %>% 
    pull(bblid) %>%
    n_distinct()
n_scans <- dim(preprocessed_cnb_persistent_dat)[1]
n_ps_scans <- preprocessed_cnb_persistent_dat %>% 
    group_by(oisPS) %>% 
    summarize(n = n()) %>% 
    filter(oisPS == "recurrent PS") %>%
    pull(n)
n_acq_types <- n_distinct(preprocessed_cnb_persistent_dat$acq_type)

## print summary of subject information
cat(glue(
    "Included subjects: {n_subjs} (recurrent PS: {n_ps})",
    "Included scans: {n_scans} (recurrent PS: {n_ps_scans})",
    "Included acquisition parameters: {n_acq_types}",
    .sep = "\n"
))

# Distribution -----------------------------------------------------------------
longitudinal_distribution <- plot_longitudinal_distribution(
    dat = preprocessed_cnb_persistent_dat, 
    x = "age_at_scan", 
    id = "bblid", 
    group = "oisPS", 
    xlab = "Age at Scan", 
    ylab = "Participant", 
    colpal = c("grey30", "tomato3")
)
longitudinal_distribution
ggsave(
    plot = longitudinal_distribution, 
    filename = here("outputs", "figs", "longitudinal_distribution.pdf"), 
    width = 4, 
    height = 6
)


# Table 1 ----------------------------------------------------------------------
sect_properties <- prop_section(
    page_size = page_size(
        orient = "landscape",
        width = 8.3, height = 11.7
    ),
    type = "continuous",
    page_margins = page_mar(gutter = 0.5)
)

table01 <- preprocessed_cnb_persistent_dat %>% 
    filter(timepoint == 1) %>% 
    select(oisPS, n_scans, age_at_scan, sex, race3, trauma_T1, envSES) %>%
    tbl_summary(
        by = oisPS,
        label = list(
            sex ~ "Sex",
            age_at_scan ~ "Age at Scan (Baseline)",
            race3 ~ "Race",
            trauma_T1 ~ "Traumatic Stress Load (Baseline)",
            envSES ~ "Environmental SES (Baseline)",
            n_scans ~ "Number of Scans"
        ),
        statistic = list(all_continuous() ~ "{mean} ({sd})"),
        type = list(
            trauma_T1 ~ "continuous",
            n_scans ~ "continuous"
        ),
        digits = list(all_continuous() ~ c(2, 2))
    ) %>% 
    add_overall() %>% 
    add_p() %>% 
    add_q() %>% 
    bold_p(q = TRUE)
table01

## save table to docx
table01 %>%
    as_flex_table() %>% 
    padding(padding.top = 0, padding.bottom = 0, part = "all") %>% 
    autofit() %>% 
    save_as_docx(
        path = here("outputs", "tables", "table01.docx"),
        pr_section = sect_properties
    )


# save acquisition sequence table
acq_parameters <- c(
    "ModelName", "RepetitionTime", "EchoTime", "InversionTime",
    "SliceThickness", "AcquisitionMatrixPE", "ReconMatrixPE", "FlipAngle", 
    "PixelBandwidth", "BaseResolution"
)
stable01 <- preprocessed_cnb_persistent_dat$acq_parameters %>% 
    unique() %>% 
    tibble(parameters = .) %>% 
    mutate(parameters = str_replace(parameters, "Prisma_fit", "Prisma fit")) %>% 
    separate(col = parameters, into = acq_parameters, sep = "_") %>% 
    arrange(ModelName) %>% 
    as_flextable() %>% 
    hline(i = 1, part = "header", border = fp_border(color = "grey30", width = 1.5)) %>% 
    delete_rows(i = 2, part = "header") %>% 
    delete_part(part = "footer") %>% 
    autofit() %>% 
    fit_to_width(10) 
 
save_as_docx(
    stable01,
    path = here("outputs", "tables", "suppltable1.docx"),
    pr_section = sect_properties
)
save_as_image(
    stable01,
    path = here("outputs", "tables", "suppltable1.png"),
    pr_section = sect_properties
)


# Save -------------------------------------------------------------------------
write_rds(
    preprocessed_cnb_persistent_dat, 
    here("data", "processed", 
         "data_cohort-extralong_status-preprocessed_decs-persistent.rds")
)
