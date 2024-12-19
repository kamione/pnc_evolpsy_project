
#' plot longitudinal distribution of the data
#' 
#' @param dat A data frame.
#' @param x A string.
#' @param group A string.
#' @param xlab A string.
#' @param ylab A string.
#' @param colpal A vector.
#' @returns A ggplot figure.
#' 
#' 
#' @importFrom magrittr "%>%"

plot_longitudinal_distribution <- function(
    dat, x, id, group, xlab, ylab, colpal
) {
    figure <- dat %>% 
        ggplot(
            aes(x = !!sym(x), 
                y = reorder(!!sym(id), !!sym(x), FUN = min),
                color = !!sym(group))
        ) +
        geom_point(size = 2, alpha = 0.75) +
        geom_line(alpha = 0.5) +
        labs(x = xlab, y = ylab, color = "", title = "") +
        scale_color_manual(values = colpal) +
        theme_classic() +
        scale_y_discrete(expand = expansion(add = 5)) +
        theme(
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            legend.position = c(0.8, 0.1)
        )
    
    return(figure)
}
