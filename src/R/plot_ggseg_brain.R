#' project values to a Desikan-Killiany brain atlas
#' 
#' @param dat A data frame.
#' @param fill A string.
#' @param title A string.
#' @param min A numeric.
#' @param max A numeric.
#' @returns A ggplot figure.
#' 
#' 
#' @importFrom magrittr "%>%"

library(ggseg)

plot_ggseg_brain <- function(
    dat, atlas = "dk", fill, title = NULL, min = -4, max = 4, break_int = 2, 
    remove_x_tick = TRUE, highlight_region = TRUE) {
    # color bar set up
    if (min < 0 & max > 0) {
        colors <- c("dodgerblue4", "white", "tomato3")
    } else {
        colors <- viridis::viridis(100)
    }
    
    switch(
        atlas,
        "dk" = {brain_data = dk},
        "aseg" = {brain_data = aseg}
    )
    
    dat <- dat %>% 
        filter(p.value < 0.05)
    
    if (highlight_region) {
        dat <- dat %>% 
            mutate(
                size = case_when(
                    p.value.adj < 0.05 ~ 1,
                    TRUE ~ 0.2
                ),
                color_group = case_when(
                    p.value.adj < 0.05 ~ "red",
                    p.value < 0.05 & p.value.adj > 0.05 ~ "grey75",
                    TRUE ~ "grey70"
                )
            )
        
        color_values <- dat %>% 
            pull(color_group) %>% 
            unique() %>% 
            sort()
    } else {
        dat <- dat %>% 
            mutate(
                size = 0.2, 
                color_group = "no"
            )
        color_values = "grey75"
    }
    
    if (atlas == "aseg") {
        fig <- dat %>%
            ggseg(
                atlas = brain_data, 
                mapping = aes(
                    fill = statistic, 
                    color = color_group,
                    size = I(size)
                ),
                view = "coronal"
            )
    } else if (atlas == "dk") {
        fig <- dat %>%
            ggseg(
                atlas = brain_data, 
                mapping = aes(
                    fill = statistic, 
                    color = color_group,
                    size = I(size)
                )
            )
    }
    fig <- fig +
        scale_fill_gradientn(
            colors = colors,
            na.value = "grey80",
            limits = c(min, max),
            breaks = seq(round(min), max, break_int)
        ) +
        scale_color_manual(values = color_values, guide = "none") +
        labs(title = title, fill = "")
    
    if (remove_x_tick) {
        fig <- fig +
            theme(
                axis.text.x = element_blank(),
                axis.title.x = element_blank()
            )
    }
    
    return(fig)
}
