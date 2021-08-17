print.summary.recEvFit <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    
    cat("\nCoefficients:\n")
    print(x$coefficients)

    cat("\nAIC: ", x$aic, "\n")
    
    cat("\nCOM-Poisson regression model for imputing missing values\n")
    print(x$CMP)
}
