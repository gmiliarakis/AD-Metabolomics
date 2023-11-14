library(caret)
library(dplyr)
# Transpose, standardize and store the metabolite data in Y
Y <- scale(ADmetabolites)

# Store the ApoE genotype in x as a factor
x <- as.factor(geno$APOE)
sex <- as.factor(geno$sex)
df <- cbind.data.frame(Y,x)

gentest <- function(Y, x, sex= NULL) {
  ### Analysis of Variance (ANOVA)
  # Store the ANOVA model summaries in summaries
  if (is.null(sex)){
    df <- cbind.data.frame(Y,x)
    summaries <- purrr::map(df[,1:230], ~summary(aov(.x ~ df$x)))
  } else {
    df <- cbind.data.frame(Y,x,sex)
    summaries <- purrr::map(df[,1:230], ~summary(aov(.x ~ df$x+df$x*df$sex)))}
  ### Correction for Multiple Testing
  #Create a list to store the p-values
  p_values <- list(1:230)
  #Extract the p-values of the F-tests from the aov summaries list and store them in p_values
  for (i in 1:230){
    p_values[[i]] <- summaries[[i]][[1]][["Pr(>F)"]][[1]]
  }
  # Coerce p_values to dataframe and transpose it 
  p_values <- t(data.frame(p_values))
  
  # Set as row names the metabolites
  p_values <- data.frame(p_values, row.names = colnames(Y))
  
  # Calculate the FDR-adjusted p-values
  p_values$p_adj <- p.adjust(p_values[,], method = "hochberg",n=230)
  
  # Filter out the non-significant (a=0.05) FDR-adjusted p-values
  return(filter(p_values, p_adj < 0.05))
}
gentest(Y, x, sex=AD$I_sex)
ADmetabolites
