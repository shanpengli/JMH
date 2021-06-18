Getriskset <- function(cdata, surv.formula) {
  
  survname <- all.vars(surv.formula)
  
  subcdata <- cdata[, survname[1:2]]
  
  subcadataRisk1 <- subcdata[subcdata[, 2] == 1, ]
  subcadataRisk2 <- subcdata[subcdata[, 2] == 2, ]
  tablerisk1 <- as.data.frame(table(subcadataRisk1[, 1]))
  tablerisk2 <- as.data.frame(table(subcadataRisk2[, 1]))
  tablerisk1$Var1 <- as.character(tablerisk1$Var1)
  tablerisk2$Var1 <- as.character(tablerisk2$Var1)
  tablerisk1$Var1 <- as.numeric(tablerisk1$Var1)
  tablerisk2$Var1 <- as.numeric(tablerisk2$Var1)
  tablerisk1 <- tablerisk1[order(tablerisk1$Var1), ]
  tablerisk2 <- tablerisk2[order(tablerisk2$Var1), ]
  tablerisk1 <- data.frame(tablerisk1, 0.0001)
  tablerisk2 <- data.frame(tablerisk2, 0.0001)
  
  subcadataRisk1 <- subcadataRisk1[order(subcadataRisk1$survtime), ]
  subcadataRisk2 <- subcadataRisk2[order(subcadataRisk2$survtime), ]
  
  tablerisk1$Var1 <- subcadataRisk1$survtime
  tablerisk2$Var1 <- subcadataRisk2$survtime
  
  colnames(tablerisk1) <- c("survtime", "d", "hazard")
  colnames(tablerisk2) <- c("survtime", "d", "hazard")
  a <- list(tablerisk1, tablerisk2)
  names(a) <- c("tablerisk1", "tablerisk2")
  return(a)
}