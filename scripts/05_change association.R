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


.generate_fig_annotation <- function(beta, pval) {
    p_value_caption <- .generate_p_caption(pval)
    result <- paste0(
        "std beta = ", round(beta, 4), ", ", p_value_caption
    )
    return(result)
}

.generate_p_caption <- function(pval) {
    if (pval < 0.001) {
        result <- "p < 0.001"
    } else {
        result <- paste("p =", round(pval, 4))
    }
    return(result)
}



# Data I/O ---------------------------------------------------------------------
project_data_path <- "data/raw"
extralong_path <- "ExtraLong/data/datafreeze-2021/TabulatedQC"

preprocessed_dat <- here(
    "data", "processed", "data_cohort-extralong_status-preprocessed_decs-persistent.rds"
) %>% 
    read_rds()


cog_labels <- preprocessed_dat %>% 
    select(Executive_CorrTraits:General_bifactor) %>% 
    colnames()
cog_labels_titles <- c(
    "Executive Functioning", "Memory", "Complex Cognition",
    "Social Cognition", "Motor", "General Cognition"
)

diff_dat <- preprocessed_dat %>% 
    group_by(bblid) %>% 
    filter(row_number() == 1 | row_number() == n()) %>% 
    mutate(
        General_bifactor_diff = General_bifactor - lag(General_bifactor),
        Executive_CorrTraits_diff = Executive_CorrTraits - lag(Executive_CorrTraits),
        Memory_CorrTraits_diff = Memory_CorrTraits - lag(Memory_CorrTraits),
        Complex_CorrTraits_diff = Complex_CorrTraits - lag(Complex_CorrTraits),
        Social_CorrTraits_diff = Social_CorrTraits - lag(Social_CorrTraits), 
        Motor_CorrTraits_diff = Motor_CorrTraits - lag(Motor_CorrTraits),
        gmv_diff = gmv - lag(gmv),
        sgmv_diff = sgmv - lag(sgmv)
    ) %>% 
    #drop_na(General_bifactor_diff:sgmv_diff) %>% 
    ungroup() %>% 
    drop_na(gmv_diff, sgmv_diff)


res_dfs <- list()
cgmv_figs <- list()

for (ith in seq_along(cog_labels)) {
    
    diff_label <- paste0(cog_labels[ith], '_diff')
    
    f1 <- glue::glue(
        "{diff_label} ~ gmv_diff * group + age_at_scan + sex + race3"
    ) %>% 
        as.formula()
    
     res <- glm(
            formula =  f1,
            data = diff_dat
        ) %>% 
        parameters::parameters(standardize = "basic") %>% 
        as_tibble() %>% 
        slice(c(2, 3, 8)) %>% 
        mutate(label = cog_labels[ith], .before = "Parameter")
    
        
    res_dfs[[ith]] <- res
    annotation <- .generate_fig_annotation(res$Std_Coefficient[1], res$p[1])
    
    cgmv_figs[[ith]] <- diff_dat %>% 
        mutate(gmv_diff = scale(gmv_diff)) %>% 
        ggplot(aes(x = gmv_diff, y = !!sym(diff_label), color = group)) +
        geom_point(size = 3, alpha = 0.7) +
        stat_smooth(method = "lm") +
        scale_color_manual(values = c("tomato3", "grey30")) +
        labs(
            x = "Cortical GMV Changes (Z)", 
            y = paste0("Differences of ", cog_labels_titles[ith], " (Z)"),
            color = "", 
            fill = ""
        ) +
        ggthemes::theme_pander() +
        theme(
            plot.margin = unit(c(5, 5, 5, 5), "mm")
        )
}

cgmv_cog_fig <- wrap_plots(
    cgmv_figs, 
    guides = "collect", 
    axis_titles = "collect", 
    nrow = 2
)
ggsave(
    plot = cgmv_cog_fig, 
    filename = here("outputs", "figs", glue::glue("cgmv_cog_diff_corr.pdf")), 
    width = 10, 
    height = 8
)


# Subcortical

subcortical_res_dfs2 <- list()
subcortial_figs <- list()

for (ith in seq_along(cog_labels)) {
    
    diff_label <- paste0(cog_labels[ith], '_diff')
    
    f1 <- glue::glue(
        "{diff_label} ~ sgmv_diff * group + age_at_scan + sex + race3"
    ) %>% 
        as.formula()
    
    res <- glm(
        formula =  f1,
        data = diff_dat
    ) %>% 
        parameters::parameters(standardize = "basic") %>% 
        as_tibble() %>% 
        slice(c(2, 3, 8)) %>% 
        mutate(label = cog_labels[ith], .before = "Parameter")
    
    
    subcortical_res_dfs2[[ith]] <- res
    annotation <- .generate_fig_annotation(res$Std_Coefficient[1], res$p[1])
    
    subcortial_figs[[ith]] <- diff_dat %>% 
        mutate(sgmv_diff = scale(sgmv_diff)) %>% 
        ggplot(aes(x = sgmv_diff, y = !!sym(diff_label), color = group)) +
        geom_point(size = 3, alpha = 0.7) +
        stat_smooth(method = "lm") +
        scale_color_manual(values = c("tomato3", "grey30")) +
        labs(
            x = "Subcortical GMV Changes (Z)", 
            y = paste0("Differences of ", cog_labels_titles[ith], " (Z)"),
            color = "", 
            fill = ""
        ) +
        ggthemes::theme_pander() +
        theme(
            plot.margin = unit(c(5, 5, 5, 5), "mm")
        )
    
}
sgmv_cog_fig <- wrap_plots(
    subcortial_figs, 
    guides = "collect", 
    axis_titles = "collect", 
    nrow = 2
)
ggsave(
    plot = sgmv_cog_fig, 
    filename = here("outputs", "figs", glue::glue("sgmv_cog_diff_corr.pdf")), 
    width = 10, 
    height = 8
)


fig <- cgmv_figs[[6]] + subcortial_figs[[6]] +
    plot_layout(
        axis_titles = "collect", 
        guides = "collect"
    ) &
    #plot_annotation(
    #    caption = "^ The models were controlled for age, sex, and race"
    #) &
    theme(
        plot.caption = element_text(size = 10),
        legend.position = "bottom"
    )
fig
ggsave(
    plot = fig, 
    filename = here("outputs", "figs", glue::glue("figure04.pdf")), 
    width = 9, 
    height = 4
)

sect_properties <- prop_section(
    page_size = page_size(
        orient = "lanscape",
        width = 8.3, height = 11.7
    ),
    type = "continuous",
    page_margins = page_mar(gutter = 0.25)
)

res_dfs %>% 
    bind_rows() %>% 
    select(-c(5, 9)) %>% 
    flextable() %>% 
    colformat_double(digits = 3) %>% 
    autofit() %>% 
    save_as_docx(
        path = here("outputs", "tables", "diff_cgmv_cognition.docx"),
        pr_section = sect_properties
    )
subcortical_res_dfs2 %>% 
    bind_rows() %>% 
    select(-c(5, 9)) %>% 
    flextable() %>% 
    colformat_double(digits = 3) %>% 
    autofit() %>% 
    save_as_docx(
        path = here("outputs", "tables", "diff_sgmv_cognition.docx"),
        pr_section = sect_properties
    )
