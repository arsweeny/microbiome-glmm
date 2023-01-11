# Output clean MCMC summaries 
require(dplyr)

#no ranef 
clean.MCMC.GLM <- function(x) {
  sols <- summary(x)$solutions  ## pull out relevant info from model summary
  Rcovs <- summary(x)$Rcovariances
  fixed <- data.frame(row.names(sols), sols, row.names = NULL)  ## convert to dataframes with the row.names as the first col
  residual <- data.frame(row.names(Rcovs), Rcovs, row.names = NULL)
  names(fixed)[names(fixed) == "row.names.sols."] <- "variable"  ## change the columns names to variable, so they all match
  names(residual)[names(residual) == "row.names.Rcovs."] <- "variable"
  fixed$effect <- "fixed"  ## add ID column for type of effect (fixed, random, residual)
  residual$effect <- "residual"
  modelTerms <- as.data.frame(bind_rows(fixed, residual))  # merge it all together
}

clean.MCMC.GLMM <- function(x) {
  sols <- summary(x)$solutions  ## pull out relevant info from model summary
  Gcovs <- summary(x)$Gcovariances
  Rcovs <- summary(x)$Rcovariances
  fixed <- data.frame(row.names(sols), sols, row.names = NULL)  ## convert to dataframes with the row.names as the first col
  random <- data.frame(row.names(Gcovs), Gcovs, row.names = NULL)
  residual <- data.frame(row.names(Rcovs), Rcovs, row.names = NULL)
  names(fixed)[names(fixed) == "row.names.sols."] <- "variable"  ## change the columns names to variable, so they all match
  names(random)[names(random) == "row.names.Gcovs."] <- "variable"
  names(residual)[names(residual) == "row.names.Rcovs."] <- "variable"
  fixed$effect <- "fixed"  ## add ID column for type of effect (fixed, random, residual)
  random$effect <- "random"
  residual$effect <- "residual"
  modelTerms <- as.data.frame(bind_rows(fixed, random, residual))  # merge it all together
}

getName.MCMC <- function(x) deparse(substitute(x))  # add the model name
