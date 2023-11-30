library(dplyr)
load("data.Rdata")


library(car)
gentest <- function(Y, x, ...) {
  ### Analysis of Variance (ANOVA)
  df <- cbind(Y, x, ...)
  ncol <- ncol(Y)
  covariates <- setdiff(names(df), c(names(Y), x))
  summaries <- purrr::map(df[, 1:ncol], ~ summary(lm(as.formula(paste0(".x ~ x +", paste(covariates, collapse = "*"))), data = df), test.statistic = "F", type = 3))

  ### Correction for Multiple Testing
  # Create a list to store the p-values
  # p_values <- list(1:ncol)
  # # Extract the p-values of the F-tests from the aov summaries list and store them in p_values
  # for (i in 1:ncol) {
  #   p_values[[i]] <- summaries[[i]][[1]][["Pr(>F)"]][[1]]
  # }
  # # Coerce p_values to dataframe and transpose it
  # p_values <- t(data.frame(p_values))

  # # Set as row names the metabolites
  # p_values <- data.frame(p_values, row.names = colnames(Y))

  # # Calculate the FDR-adjusted p-values
  # p_values$p_adj <- p.adjust(p_values[, ], method = "fdr")
  # out <- list(summaries, filter(p_values, p_adj < 0.05) )
  # Filter out the non-significant (a=0.05) FDR-adjusted p-values
  return(summaries)
}

Y <- df[, 1:230]
df$Diagnosis <- as.numeric(df$Diagnosis)
x <- as.(df$APOE)-1
sex <- as.numeric(df$sex)-1
D <- as.numeric(df$Diagnosis) -1

gentest(Y, x, sex)

Y <- subset(df, Diagnosis == 1)[,1:230]
x <- as.factor(subset(df, Diagnosis == 1)$APOE)
sex <- as.factor(subset(df, Diagnosis == 1)$sex)

gentest(Y, x)
install.packages("car")
