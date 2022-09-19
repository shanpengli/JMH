##' @title Print BrierJMMLSM
##' @name print.BrierJMMLSM
##' @aliases print.BrierJMMLSM
##' @param x x of class 'BrierJMMLSM'.
##' @param ... Further arguments passed to or from other methods.
##' @return a list of matrices with conditional probabilities for subjects.
##' @author Shanpeng Li \email{lishanpeng0913@ucla.edu}
##' @seealso \code{\link{JMMLSM}, \link{survfitJMMLSM}}
##' @export
##' 

print.BrierJMMLSM <- function (x, ...) {
  if (!inherits(x, "BrierJMMLSM"))
    stop("Use only with 'BrierJMMLSM' xs.\n") 
  
  if (is.null(x$Brier.cv)) {
    stop("The cross validation fails. Please try using a different seed number.")
  } else {
    
    if (length(x$Brier.cv) == x$n.cv && sum(mapply(is.null, x$Brier.cv)) == 0) {
      sum <- 0
      for (j in 1:x$n.cv) {
        sum <- sum + x$Brier.cv[[j]]
      }
      sum <- sum/x$n.cv
      mean.Brier1 <- sum[, 1]
      mean.Brier2 <- sum[, 2]
      
      ExpectedBrier <- data.frame(x$horizon.time, mean.Brier1, mean.Brier2)
      colnames(ExpectedBrier) <- c("Horizon Time", "Failure 1", "Failure 2")
      cat("\nExpected Brier Score at the landmark time of", x$landmark.time, "\nbased on", x$n.cv, "fold cross validation\n")
      print(ExpectedBrier)
      invisible(x)
    } else {
      stop("The cross validation fails. Please try using a different seed number.")
    }
    
  }

}