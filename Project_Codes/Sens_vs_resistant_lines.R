## Drug sensitivity calling using waterfall plots
## Method:
## 1. Sensitivity calls were made using one of IC50, ActArea or Amax
## 2. Sort log IC50s (or ActArea or Amax) of the cell lines to generate a “waterfall distribution”
## 3. Identify cutoff:
##  3.1 If the waterfall distribution is non-linear (pearson cc to the linear fit <=0.95), estimate the major inflection point of the log IC50 curve as the point on the curve with the maximal distance to a line drawn between the start and end points of the distribution.
##  3.2 If the waterfall distribution appears linear (pearson cc to the linear fit > 0.95), then use the median IC50 instead.
## 4. Cell lines within a 4-fold IC50 (or within a 1.2-fold ActArea or 20% Amax difference) difference centered around this inflection point are classified as being “intermediate”,  cell lines with lower IC50s (or ActArea/Amax values) than this range are defined as sensitive, and those with IC50s (or ActArea/Amax) higher than this range are called “insensitive”.
## 5. Require at least x sensitive and x insensitive cell lines after applying these criteria (x=5 in our case).

## Input:
##  ic50: IC50 values in micro molar (positive values)
##  actarea: Activity Area, that is area under the drug activity curve (positive values)
##  amax: Activity at max concentration (positive values)
##  intermediate.fold: vector of fold changes used to define the intermediate sensitivities for ic50, actarea and amax respectively
sensitivity.calling.waterfall <- function(x, type=c("ic50", "actarea", "amax"), intermediate.fold=c(4, 1.2, 1.2), cor.min.linear=0.95, name="Drug", plot=FALSE) {
    
    type <- match.arg(type)
    
    if (is.null(names(x))) { names(x) <- paste("X", 1:length(x), sep=".") }
    
    xx <- x[complete.cases(x)]
    switch (type,
            "ic50" = {
                xx <- -log10(xx)
                ylabel <- "-log10(IC50)"
                ## 4 fold difference around IC50 cutoff
                interfold <- log10(intermediate.fold[1])
            },
            "actarea" = {
                ylabel <- "Activity area"
                ## 1.2 fold difference around Activity Area cutoff
                interfold <- intermediate.fold[2]
            },
            "amax" = {
                ylabel <- "Amax"
                ## 1.2 fold difference around Amax
                interfold <- intermediate.fold[3]
            }
    )
    cor.min.linear=0.95
    Drug = "A00055058"
    xx <- CCLE_Drug_screening_1 %>% .[.$broad_id == Drug,] %>% pull(Sensitivity) %>% na.omit()
    oo <- order(xx, decreasing=TRUE)
    
    A00077618 <- CCLE_Drug_screening_1 %>% .[.$broad_id == Drug,] 
    ggplot(A00077618, aes(x = stripped_cell_line_name, y = Sensitivity)) + 
    geom_col()+ 
    scale_x_discrete(limits = A00077618$stripped_cell_line_name[order(A00077618$Sensitivity, decreasing = T)])
    ## test linearity with Perason correlation
    cc <- cor.test(-xx[oo], 1:length(oo), method="pearson")
    ## line between the two extreme sensitivity values
    dd <- cbind("y"=xx[oo][c(1, length(oo))], "x"=c(1, length(oo)))
    rr <- lm(y ~ x, data=data.frame(dd))
    ## compute distance from sensitivity values and the line between the two extreme sensitivity values
    ddi <- apply(cbind(1:length(oo), xx[oo]), 1, function(x, slope, intercept) {
        return(distancePointLine(x=x[1], y=x[2], slope=slope, intercept=intercept))
    }, slope=rr$coefficients[2], intercept=rr$coefficients[1])
    if(cc$estimate > cor.min.linear){
        ## approximately linear waterfall
        cutoff <- which.min(abs(xx[oo] - median(xx[oo])))
        cutoffn <- names(cutoff)[1]
    } else {
        ## non linear waterfall
        ## identify cutoff as the maximum distance
        cutoff <- which.max(abs(ddi))
        cutoffn <- names(ddi)[cutoff]
    }
    ## identify intermediate sensitivities
    switch (type,
            "ic50" = {
                rang <- c(xx[oo][cutoff] - interfold, xx[oo][cutoff] + interfold)
            },
            "actarea" = {
                rang <- c(xx[oo][cutoff] / interfold, xx[oo][cutoff] * interfold)
            },
            "amax" = {
                rang <- c(xx[oo][cutoff] / interfold, xx[oo][cutoff] * interfold)
            }
    )
    ## compute calls
    calls <- rep(NA, length(xx))
    names(calls) <- names(xx)
    calls[xx < rang[1]] <- "resistant"
    calls[xx > rang[2]] <- "sensitive"
    calls[xx >= rang[1] & xx <= rang[2]] <- "intermediate"
    calls %>% table()
    tt <- rep(NA, length(x))
    names(tt) <- names(x)
    tt[names(calls)] <- calls
     
}

