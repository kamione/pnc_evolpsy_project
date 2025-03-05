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
library(flextable)
library(ftExtra)

source(here("src", "R", "plot_smoothbyfactor.R"))
source(here("src", "R", "plot_ggseg_brain.R"))
source(here("src", "R", "get_gamm_partialR2.R"))

sect_properties <- prop_section(
    page_size = page_size(
        orient = "landscape",
        width = 8.3, height = 11.7
    ),
    type = "continuous",
    page_margins = page_mar(gutter = 0.5)
)



# Data I/O ---------------------------------------------------------------------
project_data_path <- "data/raw"
extralong_path <- "ExtraLong/data/datafreeze-2021/TabulatedQC"

preprocessed_dat <- here(
    "data", "processed", 
    "data_cohort-extralong_status-preprocessed_decs-persistent.rds"
) %>% 
    read_rds()


# CNB --------------------------------------------------------------------------
cog_labels <- preprocessed_dat %>% 
    select(Executive_CorrTraits:General_bifactor) %>% 
    colnames()
cog_labels_titles <- c(
    "Executive Functioning", "Memory", "Complex Cognition",
    "Social Cognition", "Motor", "General Cognition"
)

cog_dev_res <- list()
cog_dev_figs <- list()

confounds <- "sex + race3 + timepoint + acq_type"

for (iter in 1:length(cog_labels)) {
    
    forumla_cnb <- as.formula(
        glue(
            "{cog_labels[iter]} ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + 
             s(age_at_scan, k = 4, by = oisPS, fx = TRUE) + {confounds}"
        )
    )
    
    mod <- gamm(
        formula = forumla_cnb,
        random = list(bblid =~ 1),
        data = preprocessed_dat,
        REML = TRUE
    )
    
    cog_dev_res[[iter]] <- mod %>% 
        parameters::parameters() %>% 
        as_tibble() %>% 
        mutate(label = cog_labels_titles[iter], .before = "Parameter")
    
    cog_dev_figs[[iter]] <- plot_smoothbyfactor(
        modobj = mod,  
        series = "age_at_scan", 
        by = "oisPS",
        subjid = "bblid",
        title = cog_labels_titles[iter],
        show_differences = FALSE,
        title_font_size = 18,
        ylabel = "Cognitive Score (Z)",
        legend.position.y = 0.3
    ) +
        theme(
            plot.title = element_text(hjust = 0.2),
            axis.text.x = element_text(),
            axis.title.x = element_text()
        )
}

cog_dev <- wrap_plots(cog_dev_figs, guides = "collect") +
    plot_layout(axis_titles = "collect")
cog_dev
ggsave(
    plot = cog_dev, 
    filename = here("outputs", "figs", glue::glue("cog_dev.pdf")), 
    width = 13, 
    height = 7
)

cog_dev_res_table <- cog_dev_res %>% 
    bind_rows() %>% 
    filter(Parameter %in% c("oisPS.L", "s(age_at_scan)", "s(age_at_scan):oisPSrecurrent PS")) %>% 
    select(label, Parameter, `t / F`, p) %>% 
    rename(statistic = `t / F`) %>% 
    pivot_wider(names_from = Parameter, values_from = statistic:p) %>% 
    mutate(
        padj_oisPS.L = p.adjust(p_oisPS.L, "fdr"),
        `padj_s(age_at_scan)` = p.adjust(`p_s(age_at_scan)`, "fdr"),
        `padj_s(age_at_scan):oisPSrecurrent PS` = p.adjust(`p_s(age_at_scan):oisPSrecurrent PS`, "fdr")
    ) %>% 
    select(1, 2, 5, 8, 3, 6, 9, 4, 7, 10) %>% # reorder
    mutate(
        `p_s(age_at_scan)` = "< 0.001",
        `padj_s(age_at_scan)` = "< 0.001"
    ) %>% 
    flextable() %>% 
    colformat_double(
        j = c(2, 5, 8),
        big.mark = ",", 
        digits = 2
    ) %>% 
    colformat_double(
        j = c(3, 4, 9, 10),
        big.mark = ",", 
        digits = 3
    ) %>% 
    add_header_row(
        values = c("", "Group (TD - recurrent PS)", "s(Age)", "s(Age) by Group"), 
        colwidths = c(1, 3, 3, 3)
    ) %>% 
    set_header_labels(
        label = "Term",
        statistic_oisPS.L = "t",
        `statistic_s(age_at_scan)` = "F",
        `statistic_s(age_at_scan):oisPSrecurrent PS` = "F"
    ) %>% 
    compose(
        i = 2, j = c(3, 6, 9), part = "header",
        value = as_paragraph(
            "p"
        )
    ) %>% 
    compose(
        i = 2, j = c(4, 7, 10), part = "header",
        value = as_paragraph(
            "q"
        )
    ) %>% 
    compose(
        i = 3, j = 3, part = "body",
        value = as_paragraph(
            "< 0.001"
        )
    ) %>% 
    align(align = "center", part = "header") %>% 
    bold(part = "header") %>% 
    italic(i = 2, italic = TRUE, part = "header") %>% 
    bold(~ padj_oisPS.L < 0.05, c(4)) %>% 
    bold(~ `padj_s(age_at_scan)` == "< 0.001", c(7)) %>% 
    add_footer_lines(
        "^ Models controlling for sex, race and time point; significant results are bold"
    ) %>% 
    autofit()

cog_dev_res_table

cog_dev_res_table %>% 
    save_as_docx(
        path = here("outputs", "tables", "table02.docx"),
        pr_section = sect_properties
    )

# Sensitive Analyses -----------------------------------------------------------
cog_dev_sens_res <- list()
cog_dev_sens_figs <- list()

confounds_sens <- "sex + race3 + timepoint + acq_type + n_scans + trauma_T1 + envSES"

for (iter in 1:6) {
    forumla_cnb <- as.formula(
        glue(
            "{cog_labels[iter]} ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + 
             s(age_at_scan, k = 4, by = oisPS, fx = TRUE) + {confounds_sens}"
        )
    )
    mod <- gamm(
        formula = forumla_cnb,
        random = list(bblid =~ 1),
        data = preprocessed_dat,
        REML = TRUE
    )
    cog_dev_sens_res[[iter]] <- mod %>% 
        parameters::parameters() %>% 
        as_tibble() %>% 
        mutate(label = cog_labels_titles[iter], .before = "Parameter")
    
}

cog_dev_sens_res_table <- cog_dev_sens_res %>% 
    bind_rows() %>% 
    filter(Parameter %in% c("oisPS.L", "s(age_at_scan)", "s(age_at_scan):oisPSrecurrent PS")) %>% 
    select(label, Parameter, `t / F`, p) %>% 
    rename(statistic = `t / F`) %>% 
    pivot_wider(names_from = Parameter, values_from = statistic:p) %>% 
    mutate(
        padj_oisPS.L = p.adjust(p_oisPS.L, "fdr"),
        `padj_s(age_at_scan)` = p.adjust(`p_s(age_at_scan)`, "fdr"),
        `padj_s(age_at_scan):oisPSrecurrent PS` = p.adjust(`p_s(age_at_scan):oisPSrecurrent PS`, "fdr")
    ) %>% 
    select(1, 2, 5, 8, 3, 6, 9, 4, 7, 10) %>% # reorder
    mutate(
        `p_s(age_at_scan)` = "< 0.001",
        `padj_s(age_at_scan)` = "< 0.001"
    ) %>% 
    flextable() %>% 
    colformat_double(
        j = c(2, 5, 8),
        big.mark = ",", 
        digits = 2
    ) %>% 
    colformat_double(
        j = c(3, 4, 9, 10),
        big.mark = ",", 
        digits = 3
    ) %>% 
    add_header_row(
        values = c("", "Group (TD - recurrent PS)", "s(Age)", "s(Age) by Group"), 
        colwidths = c(1, 3, 3, 3)
    ) %>% 
    set_header_labels(
        label = "Term",
        statistic_oisPS.L = "t",
        `statistic_s(age_at_scan)` = "F",
        `statistic_s(age_at_scan):oisPSrecurrent PS` = "F"
    ) %>% 
    compose(
        i = 2, j = c(3, 6, 9), part = "header",
        value = as_paragraph(
            "p"
        )
    ) %>% 
    compose(
        i = 2, j = c(4, 7, 10), part = "header",
        value = as_paragraph(
            "q"
        )
    ) %>% 
    compose(
        i = 3, j = 3, part = "body",
        value = as_paragraph(
            "< 0.001"
        )
    ) %>% 
    align(align = "center", part = "header") %>% 
    bold(part = "header") %>% 
    italic(i = 2, italic = TRUE, part = "header") %>% 
    bold(~ padj_oisPS.L < 0.05, c(4)) %>% 
    bold(~ `padj_s(age_at_scan)` == "< 0.001", c(7)) %>% 
    add_footer_lines(
        "^ Models controlling for sex, race, time point, acquisition types, number of scans, traumatic stress load at baseline and envSES; significant results are bold"
    ) %>% 
    autofit()

cog_dev_sens_res_table

cog_dev_sens_res_table %>% 
    save_as_docx(
        path = here("outputs", "tables", "cog_dev_sens_res_table.docx"),
        pr_section = sect_properties
    )


# Sensitivity Analyses 2 -------------------------------------------------------
cog_dev_sens_res <- list()
cog_dev_sens_figs <- list()

confounds_sens <- "sex + race3 + timepoint"

# filter out 
for (iter in 1:6) {
    forumla_cnb <- as.formula(
        glue(
            "{cog_labels[iter]} ~ oisPS + s(age_at_scan, k = 4, fx = TRUE) + 
             s(age_at_scan, k = 4, by = oisPS, fx = TRUE) + {confounds_sens}"
        )
    )
    mod <- gamm(
        formula = forumla_cnb,
        random = list(bblid =~ 1),
        data = preprocessed_dat %>% filter(!age_at_scan < 12.5),
        REML = TRUE
    )
    cog_dev_sens_res[[iter]] <- mod %>% 
        parameters::parameters() %>% 
        as_tibble() %>% 
        mutate(label = cog_labels_titles[iter], .before = "Parameter")
    
}

cog_dev_sens_res_table <- cog_dev_sens_res %>% 
    bind_rows() %>% 
    filter(Parameter %in% c("oisPS.L", "s(age_at_scan)", "s(age_at_scan):oisPSrecurrent PS")) %>% 
    select(label, Parameter, `t / F`, p) %>% 
    rename(statistic = `t / F`) %>% 
    pivot_wider(names_from = Parameter, values_from = statistic:p) %>% 
    mutate(
        padj_oisPS.L = p.adjust(p_oisPS.L, "fdr"),
        `padj_s(age_at_scan)` = p.adjust(`p_s(age_at_scan)`, "fdr"),
        `padj_s(age_at_scan):oisPSrecurrent PS` = p.adjust(`p_s(age_at_scan):oisPSrecurrent PS`, "fdr")
    ) %>% 
    select(1, 2, 5, 8, 3, 6, 9, 4, 7, 10) %>% # reorder
    mutate(
        `p_s(age_at_scan)` = "< 0.001",
        `padj_s(age_at_scan)` = "< 0.001"
    ) %>% 
    flextable() %>% 
    colformat_double(
        j = c(2, 5, 8),
        big.mark = ",", 
        digits = 2
    ) %>% 
    colformat_double(
        j = c(3, 4, 9, 10),
        big.mark = ",", 
        digits = 3
    ) %>% 
    add_header_row(
        values = c("", "Group (TD - recurrent PS)", "s(Age)", "s(Age) by Group"), 
        colwidths = c(1, 3, 3, 3)
    ) %>% 
    set_header_labels(
        label = "Term",
        statistic_oisPS.L = "t",
        `statistic_s(age_at_scan)` = "F",
        `statistic_s(age_at_scan):oisPSrecurrent PS` = "F"
    ) %>% 
    compose(
        i = 2, j = c(3, 6, 9), part = "header",
        value = as_paragraph(
            "p"
        )
    ) %>% 
    compose(
        i = 2, j = c(4, 7, 10), part = "header",
        value = as_paragraph(
            "q"
        )
    ) %>% 
    compose(
        i = 3, j = 3, part = "body",
        value = as_paragraph(
            "< 0.001"
        )
    ) %>% 
    align(align = "center", part = "header") %>% 
    bold(part = "header") %>% 
    italic(i = 2, italic = TRUE, part = "header") %>% 
    bold(~ padj_oisPS.L < 0.05, c(4)) %>% 
    bold(~ `padj_s(age_at_scan)` == "< 0.001", c(7)) %>% 
    add_footer_lines(
        "^ Models controlling for sex, race, time point; significant results are bold"
    ) %>% 
    autofit()

cog_dev_sens_res_table

cog_dev_sens_res_table %>% 
    save_as_docx(
        path = here("outputs", "tables", "cog_dev_sens2_res_table.docx"),
        pr_section = sect_properties
    )
