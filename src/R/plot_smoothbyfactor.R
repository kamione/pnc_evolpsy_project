library(tidygam)
library(ggplot2)
library(ggthemes)
library(patchwork)


plot_smoothbyfactor <- function(
    modobj, 
    series = NULL, 
    series2 = NULL,
    by, 
    subjid = NULL, 
    show_data = TRUE, 
    n_preds = 1000,
    title = NULL,
    title_font_size = 12,
    ylabel = "Volume",
    show_differences = TRUE,
    legend.position.x = 0.95,
    legend.position.y = 0.95
) {
    # check
    if (identical(class(modobj), c("gamm", "list"))) {
        model <- modobj$gam
        is_gamm <- TRUE
    } else if (class(modobj) == "gam") {
        model <- modobj
    } else {
        stop("modobj must be a gam/gamm model object!")
    }
    
    dat <- model$model
    
    dv <- model$formula[[2]]
    ivs <- attr(model$terms, "term.labels")
    ivs_class <- attr(model$terms, "dataClasses")
    formula_text <- model$formula |>
        deparse(width.cutoff = 500) |>
        paste(collapse = "") |>
        stringr::str_replace_all(" ", "")
    # check 
    if (!grepl(x = formula_text, pattern = paste0("by=", by))) {
        stop("no by in model formula!")
    }

    confound_values <- list()
    confound_ivs <- setdiff(ivs, c(series, series2, by))
    for (ith_confound in confound_ivs) {
        switch(
            ivs_class[ith_confound],
            "numeric" = {confound_values[[ith_confound]] = mean(dat[[ith_confound]])},
            "factor" = {confound_values[[ith_confound]] = levels(dat[[ith_confound]])[1]},
            "ordered" = {confound_values[[ith_confound]] = levels(dat[[ith_confound]])[1]}
        )
    }
    
    compare_list <- list()
    
    
    if (!is.null(series2)) {
        series2_quantile <- dat |> 
            filter(timepoint == 1) %>% 
            pull(series2) |>
            quantile(probs = c(0.1, 0.9)) %>% 
            as.vector()
        
        compare_list[[series2]] <- series2_quantile
        
        pred_smooth_dat <- list()
        pred_diff_dat <- list()
        
        for (ith_tile in series2_quantile) {
            confound_values[[by]] <- levels(dat$oisPS)[1]
            pred_smooth_dat[[as.character(ith_tile)]] <- predict_gam(
                model = model, 
                series = series,
                values = confound_values,
                length_out = n_preds
            )
            pred_diff_dat[[as.character(ith_tile)]] <- get_difference(
                model = model,
                series = series,
                compare = compare_list,
                values = confound_values, 
                length_out = n_preds
            )
        }
        pred_smooth_dat <- dplyr::bind_rows(pred_smooth_dat)
        pred_diff_dat <- dplyr::bind_rows( pred_diff_dat)
        
        pred_smooth_dat[[series2]] <- factor(
            pred_smooth_dat[[series2]], label = c("Low ERS", "High ERS")
        )
        pred_diff_dat[[series2]] <- factor(
            pred_diff_dat[[series2]], label = c("Low ERS", "High ERS")
        )
    } else {
        compare_list[[by]] <- levels(dat[[by]])
        pred_smooth_dat <- predict_gam(
            model = model, 
            series = series,
            values = confound_values,
            length_out = n_preds
        )
        pred_diff_dat <- get_difference(
            model = model,
            series = series,
            compare = compare_list,
            values = confound_values, 
            length_out = n_preds
        )
    }
    
    smooth_fig <- pred_smooth_dat |>
        ggplot(aes(x = .data[[series]], y = .data[[dv]], group = .data[[by]])) +
        geom_line(aes(color = .data[[by]]), linewidth = 1.5) +
        geom_ribbon(
            aes(x = .data[[series]], ymin = lower_ci, ymax = upper_ci, fill = .data[[by]]), 
            alpha = 0.2, 
            linetype = 0
        ) +
        labs(x = "Age at Scan", y = glue::glue("{ylabel}"), color = "", 
             fill = "", title = title) +
        scale_color_manual(values = c("grey30", "tomato3")) +
        scale_fill_manual(values = c("grey30", "tomato3")) +
        ggthemes::theme_pander() +
        theme(
            legend.justification = c(1, 1), 
            legend.position = c(legend.position.x, legend.position.y),
            legend.background = element_rect(color = NA, fill = NA),
            plot.title = element_text(
                size = title_font_size, 
                face = "plain", 
                family = "mono",
            ),
            plot.title.position = "plot",
            axis.text.x = element_blank(),
            axis.title.x = element_blank()
        )
    if (show_data & is.null(series2)) {
        smooth_fig <- smooth_fig +
            geom_point(
                data = dat, 
                aes(x = .data[[series]], y = .data[[dv]], group = .data[[by]], color = .data[[by]]),
                alpha = 0.3,
                size = 2.5
            )
    }
    if (!is.null(subjid) & is.null(series2)) {
        smooth_fig <- smooth_fig +
            geom_line(
                data = dat, 
                aes(x = .data[[series]], y = .data[[dv]], group = .data[[subjid]], color = .data[[by]]),
                alpha = 0.3
            )
    }
    sig_int <- attr(pred_diff_dat, "sig_int")
    diff_fig <- pred_diff_dat |>
        ggplot(aes(x = .data[[series]], y = diff)) +
        geom_line(color = "grey30", linewidth = 0.5) +
        geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), fill = "grey30", alpha = 0.25) +
        geom_hline(yintercept = 0, alpha = 0.5) +
        labs(x = "Age at Scan (Year)", y = "Difference") +
        theme_pander()
    
    if (!is.null(sig_int)) {
        diff_fig <- diff_fig +
            annotate(
                "rect",
                xmin = sig_int$start, 
                xmax = sig_int$end,
                ymin = -Inf, 
                ymax = Inf, 
                alpha = 0.25,
                fill = "goldenrod2"
            )
    }
    
    if (!is.null(series2)) {
        smooth_fig <- ggpubr::facet(smooth_fig, facet.by = series2)
        diff_fig <- ggpubr::facet(diff_fig, facet.by = series2) +
            theme(
                strip.text.x = element_blank(), 
                strip.background = element_blank()
            )
    }
    if (show_differences) {
        combined_fig <- smooth_fig / plot_spacer() / diff_fig +
            plot_layout(heights = c(4, -0.1, 1))
    } else {
        combined_fig <- smooth_fig
    }
    
    return(combined_fig)
}
