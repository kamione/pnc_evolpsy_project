.convert_df2ft <- function(df, footnote = NULL, filename = NULL) {
    ft <- df |>
        set_names(c("Desikan-Killiany (DK) Atlas Region", "Statistic", "p", "FDR-adjusted p")) |>
        flextable() |>
        colformat_double(digits = 3) |> 
        bold(part = "header") |> 
        align(align = "center", part = "header") |>
        compose(j = 4, part = "header", value = as_paragraph("p", as_sub("FDR"))) |> 
        fontsize(size = 8, part = "all") |>
        autofit()
    
    if (!is.null(footnote)) {
        ft <- ft |>
            footnote(
                i = 1, 
                j = 2:4,
                ref_symbols = "a",
                part = "header",
                value = as_paragraph(footnote)
            ) |>
            fontsize(size = 8, part = "footer")
    }
    
    if (is.null(filename)) {
        return(ft)
    }
    
    sect_properties <- prop_section(
        page_size = page_size(
            orient = "portrait",
            width = 8.3, height = 11.7
        ),
        type = "continuous",
        page_margins = page_mar(gutter = 0.25)
    )
    
    ft |> 
        save_as_docx(
            path = filename,
            pr_section = sect_properties
        )
    return(invisible(ft))
}
