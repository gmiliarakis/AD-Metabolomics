library(dplyr)
load("data/data.Rdata")
load("data/clinical.Rdata")
library(mice)
library(ARTool)
library(robustbase)
imp <- mice(clinical, method = "norm", seed = 123)
clinical <- complete(imp)

gentest <- function(Y, x, ...) {
  ### Analysis of Variance (ANOVA)
  df <- cbind(Y, x, ...)
  ncol <- ncol(Y)
  covariates <- names(...)
  summaries <- purrr::map(df[, 1:ncol], ~ {
    frm <- as.formula(paste0(".x ~", paste(covariates, collapse = "+")))
    rm <- lm(frm, data = df)
    ffm <- as.formula(paste0(".x ~ x +", paste(covariates, collapse = "+")))
    fm <- lm(ffm, data = df)
    anova(rm, fm)
  })
  ## Correction for Multiple Testing
  # Create a list to store the p-values
  p_values <- list(1:ncol)
  # Extract the p-values of the F-tests from the aov summaries list and store them in p_values
  for (i in 1:ncol) {
    p_values[[i]] <- summaries[[i]][["Pr(>F)"]][[2]]
  }
  # Coerce p_values to dataframe and transpose it
  p_values <- t(data.frame(p_values))

  # Set as row names the metabolites
  p_values <- data.frame(p_values, row.names = colnames(Y))

  # Calculate the FDR-adjusted p-values
  p_values$p_adj <- p.adjust(p_values[, ], method = "fdr")
  # Filter out the non-significant (a=0.05) FDR-adjusted p-values
  return(dplyr::filter(p_values, p_adj < 0.05))
}

Y <- df[, 1:230]
df$Diagnosis <- as.numeric(df$Diagnosis)
x <- as.factor(df$APOE)
e4 <- as.factor(df$E4)

gentest(Y, e4, clinical)

Y <- subset(df, Diagnosis == "Probable AD")[,1:230]
x <- as.factor(subset(df, Diagnosis == "Probable AD")$APOE)
covariates = subset(clinical, df$Diagnosis == "Probable AD")

gentest(Y, x, covariates)
clnr <- as.data.frame(purrr::map(clinical[, 1:ncol(clinical)], ~ as.numeric(.x)))
clnr <- as.matrix(clnr)
cor <- cor(clnr)
cor <- na.omit(cor)
cor[14,] <- NULL
heatmaply::heatmaply_cor(cor)
anovarob <- function(Y, x, ...) {
  ### Analysis of Variance (ANOVA)
  df <- cbind(Y, x, ...)
  ncol <- ncol(Y)
  covariates <- names(...)
  summaries <- purrr::map(df[, 1:ncol], ~ {
    frm <- as.formula(paste0(".x ~", paste(covariates, collapse = "+")))
    rm <- lmrob(frm, data = df, setting = "KS2014")
    ffm <- as.formula(paste0(".x ~ x +", paste(covariates, collapse = "+")))
    fm <- lmrob(ffm, data = df, setting = "KS2014")
    stats::anova(rm, fm, test = "Wald")
  })

  ## Correction for Multiple Testing
  # Create a list to store the p-values
  p_values <- list(1:ncol)
  # Extract the p-values of the F-tests from the aov summaries list and store them in p_values
  for (i in 1:ncol) {
    p_values[[i]] <- summaries[[i]][["Pr(>chisq)"]][[2]]
  }
  # Coerce p_values to dataframe and transpose it
  p_values <- t(data.frame(p_values))

  # Set as row names the metabolites
  p_values <- data.frame(p_values, row.names = colnames(Y))

  # Calculate the FDR-adjusted p-values
  p_values$p_adj <- p.adjust(p_values[, ], method = "fdr")
  # Filter out the non-significant (a=0.05) FDR-adjusted p-values
  return(dplyr::filter(p_values, p_adj < 0.05))
}

anovarob(Y, x, clinical)
