library(broom.mixed)

get_gamm_partialR2 <- function(full_mod, reduced_mod) {
    
    if(!identical(class(full_mod), c("gamm", "list"))) {
        stop("This function only accepts mgcv::gamm output")
    }
    if(!identical(class(reduced_mod), c("gamm", "list"))) {
        stop("This function only accepts mgcv::gamm output")
    }
    
    y <- augment(full_mod$lme)$y
    y_resid <- augment(full_mod$lme)$.resid
    y_reduced <- augment(reduced_mod$lme)$y
    y_resid_reduced <- augment(reduced_mod$lme)$.resid
    
    SSres = sum((y_resid)^2)
    SStot = sum((y - mean(y))^2)
    SSreg = SStot - SSres
    # r2 = 1 - (SSres / SStot)
    
    SSres_reduced = sum((y_resid_reduced)^2)
    SStot_reduced = sum((y_reduced - mean(y_reduced))^2)
    SSreg_reduced = SStot_reduced - SSres_reduced
    # r2_reduced = 1 - (SSres_reduced / SStot_reduced)
    
    partial_r2 <- (SSreg - SSreg_reduced) / SSres_reduced
    
    return(partial_r2)
}
