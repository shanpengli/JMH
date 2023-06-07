#' Simulated longitudinal data
#'
#' @description The \code{ydata} data frame has 1353 rows and 6 columns.
#'
#' @format This data frame contains the following columns:
#'
#'   \describe{
#'
#'   \item{\code{ID}}{patient identifier.}
#'
#'   \item{\code{Y}}{response variable.}
#'
#'   \item{\code{time}}{visit time.}
#'
#'   \item{\code{Z1}}{treatment indicator. \code{0} denotes the placebo group and \code{1} the treatment group.}
#'    
#'    \item{\code{Z2}}{continuous variable..}
#'    
#'    \item{\code{Z3}}{continuous variable..}
#'   }
#' @usage data(ydata)
#'
"ydata"

#' Simulated competing risks data
#'
#' @description The \code{cdata} data frame has 200 rows and 6 columns.
#'
#' @format This data frame contains the following columns:
#'
#'   \describe{
#'
#'   \item{\code{ID}}{patient identifier.}
#
#'   \item{\code{survtime}}{event time.}
#'
#'   \item{\code{cmprsk}}{event indicator. \code{0} denotes censoring, \code{1} risk 1, 
#'   and \code{2} risk 2.}
#'
#'   \item{\code{var1}}{treatment indicator. \code{0} denotes the placebo group and \code{1} the treatment group.}
#'
#'   \item{\code{var2}}{continuous variable.}
#'
#'   \item{\code{var3}}{continuous variable.}
#'   }
#' @usage data(cdata)
#'
"cdata"