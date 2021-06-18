
sortdata <- function(cdata, ydata, ID, surv.formula, long.formula) {
  
  Tdata <- dplyr::left_join(ydata, cdata, by = ID)
  long <- all.vars(long.formula)
  survival <- all.vars(surv.formula)
  surv <- survival[1]
  
  Tdata <- Tdata[order(-Tdata[, surv], Tdata[, ID]), ]
  
  cdata <- unique(Tdata[, c(ID, survival)])
  ydata <- Tdata[, c(ID, long)]
  #ydata[, ID] <- as.character(ydata[, ID])
  #cdata[, ID] <- as.character(cdata[, ID])
  mdata <- as.data.frame(table(ydata[, ID]))
  colnames(mdata)[1] <- ID
  mdata[, ID] <- as.character(mdata[, ID])
  cdata[, ID] <- as.character(cdata[, ID])
  ydata[, ID] <- as.character(ydata[, ID])
  cmdata <- dplyr::left_join(cdata, mdata, by = ID)
  mdata <- cmdata[, c(1, ncol(cmdata))]
  colnames(mdata) <- c(ID, "ni")
  
  a <- list(ydata, cdata, mdata)
  names(a) <- c("ydata", "cdata", "mdata")
  
  return(a)
  
  
}
