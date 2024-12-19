# check packages
wants <- c("tidyverse", "here", "ggpubr", "ggthemes", "gtsummary", "glue",
           "lmerTest", "jtools", "haven", "skimr", "correlation", "emmeans",
           "sjPlot", "lme4", "patchwork", "tidyquant", "janitor", "officer")

has <- wants %in% rownames(installed.packages(dependencies = TRUE))

if (any(!has) == 0) {
    cat("You have all the required pacakges")
} else {
    cat("The following packages are missing:\n\n")
    cat(paste0(unlist(wants[!has]), collapse = "\n"))
    cat("\n\n")
    
    ans <- menu(c("Yes", "No"), title = "Do you want to install the requierd packages?")
    
    if (ans == 2) {
        cat("Please make sure that you have installed the required packages before running the main sripts.")
    } else {
        install.packages(wants[!has])
    }
    
}
