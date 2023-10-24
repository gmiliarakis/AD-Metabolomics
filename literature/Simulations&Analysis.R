################################################################################
################################################################################
################################################################################
##
## Name:   Simulations&Analysis.R
## Author: Carel F.W. Peeters
##         Statistics for Omics Research Unit
##         Dept. of Epidemiology & Biostatistics
##         Amsterdam Public Health research institute
##         VU University medical center
##         Amsterdam, the Netherlands
## Email:	 cf.peeters@amsterdamumc.nl
##
## Last Update:	31/10/2019
## Description:	R script for the simulations and analyses contained in the
##              manuscript "Stable prediction with radiomics data"
##              https://arxiv.org/abs/1903.11696
##
################################################################################
################################################################################
################################################################################


#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Preliminaries**
#' **------------------------------------------------------------------------**

## Set working directory
setwd("C:/Users/cf.peeters/Desktop/Radiomics_pipeline/Code")

## Load FMradio
require("FMradio")

## Load additional packages
require("ggplot2")
require("reshape")
require("MASS")
require("Biobase")
require("rags2ridges")
require("DandEFA")
require("randomForestSRC")
require("pec")
require("CoxBoost")
require("party")
require("penalized")




#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Simulations on Dimensionality Assessment**
#' **Section 2 of the Supplementary Materials**
#' **------------------------------------------------------------------------**

#'
#'**---------------------------------**\
#'**Supporting Function*\
#'**---------------------------------**\

DimSim <- function(nsim, maxdim, seed = NULL, ## General arguments
                   p, m, n, simplestructure = TRUE, balanced = TRUE, ## FAsim arguments
                   loadingfix = TRUE, loadingnegative = TRUE,
                   loadingvalue = .8, loadingvaluelow = .2, numloadings,
                   loadinglowerH = .7, loadingupperH = .9, 
                   loadinglowerL = .1, loadingupperL = .3,
                   alpha = .05, Bartlett = FALSE){
  ##############################################################################
  # Simulation support for dimensionality assessment 
  #
  # General arguments:
  # nsim    > number of simulation repetitions
  # maxdim  > maximum number of latent factors to be assessed
  # seed    > starting seed for simulations
  #
  # Arguments to FAsim:
  # - p               > feature dimension
  # - m               > dimension of latent vector
  # - n               > number of samples
  # - simplestructure > logical indicating if factor structure should be 
  #                     factorially pure
  # - balanced        > logical indicating if the 'significant' loadings
  #                     should be divided evenly over the respective factors
  # - loadingfix      > logical indicating if the loadings should have a 
  #                     fixed value
  # - loadingnegative > logical indicating if, next to positive, also negative
  #                     loadings should be present
  # - loadingvalue    > value for high loadings, used when loadingfix = TRUE
  # - loadingvaluelow > value for low loadings, used when loadingfix = TRUE &
  #                     simplestructure = FALSE
  # - numloadings     > vector with length equalling argument m, indicating the
  #                     number of 'significant' loadings per factor. Used when
  #                     balanced = FALSE
  # - loadinglowerH   > lower-bound of 'significant' (high) loading, used when 
  #                     loadingfix = FALSE
  # - loadingupperH   > upper-bound of 'significant' (high) loading, used when 
  #                     loadingfix = FALSE
  # - loadinglowerL   > lower-bound of 'non-significant' (low) loading, used 
  #                     when loadingfix = FALSE & simplestructure = FALSE
  # - loadingupperL   > upper-bound of 'non-significant' (low) loading, used 
  #                     when loadingfix = FALSE & simplestructure = FALSE
  #
  # Arguments shared between dimIC, dimLRT, and dimGB and passed internally:
  # R       > (regularized) covariance or correlation matrix
  # verbose > logical indicating if function should run silently
  #
  # Arguments shared between dimIC, and dimLRT and passed internally:
  # graph   > logical indicating if output should also be visualized
  #
  # Arguments to dimIC, passed internally:
  # Type    > character specifying the penalty type: either BIC or AIC
  # n       > sample size
  # maxdim  > maximum number of latent factors to be assessed
  #
  # Arguments to dimLRT, passed internally:
  # X       > centered and scaled (subsetted) data matrix, observations in rows
  # maxdim  > maximum number of latent factors to be assessed
  # rankDOF > logical indicating if the degrees of freedom should be based on
  #           the rank of the raw correlation matrix
  #
  # Arguments to dimLRT:
  # alpha    > numeric giving alpha level
  # Bartlett > logical indicating if Bartlett correction should be applied
  #
  # NOTES:
  # - maxdim cannot exceed the Ledermann-bound
  # - usually, one wants to set maxdim (much) lower than the Ledermann-bound
  # - function checks if the p to n ratio will be such that negative degrees
  #   of freedom ensue
  # - note that, when loadingfix = FALSE, the loadings and uniquenesses 
  #   matrices are NOT treated as fixed in each simulation run
  ##############################################################################  
  
  # Dependencies:
  # require("base")
  # require("stats")
  # require("MASS")
  
  # Preliminaries for checks
  mmax <- floor((2*p + 1 - sqrt(8*p + 1))/2)
  
  # Checks
  if (maxdim > mmax){
    stop("Input (maxdim) is too high")
  }
  
  # Making sure there is a seed
  if (is.null(seed)){
    seed <- runif(1, 1, 1e5)
  }
  
  # Further preliminaries
  AICVEC   <- numeric()
  BICVEC   <- numeric()
  LRTVEC   <- numeric()
  GBVEC    <- numeric()
  Comm     <- rep(0, p)
  KMOs     <- 0
  verbose  <- FALSE
  graph    <- FALSE
  
  # Run simulation
  set.seed(seed = seed)
  for (i in 1:nsim){
    # Generate data and raw correlation matrix
    GD <- FAsim(p = p, m = m, n = n, simplestructure = simplestructure, 
                balanced = balanced, loadingfix = loadingfix, 
                loadingnegative = loadingnegative,
                loadingvalue = loadingvalue, 
                loadingvaluelow = loadingvaluelow, 
                numloadings = numloadings,
                loadinglowerH = loadinglowerH, 
                loadingupperH = loadingupperH, 
                loadinglowerL = loadinglowerL, 
                loadingupperL = loadingupperL)
    DATAscaled <- GD$data
    Re         <- GD$cormatrix
    
    # Redundancy filtering
    # Will, in most cases, remove little
    Re         <- RF(Re, t = .95)
    DATAscaled <- subSet(DATAscaled, Re)
    
    # Obtain optimal penalty and matrix
    OPT <- regcor(DATAscaled, fold = 5, verbose = verbose)
    Re  = OPT$optCor
    
    # Determine optimal dimension according to AIC
    AIC <- dimIC(Re, n = n, maxdim = maxdim, Type = "AIC", 
                 graph = graph, verbose = verbose)
    AICVEC[i] <- AIC[which(AIC[,2] == min(AIC[,2])),1]
    
    # Determine optimal dimension according to BIC
    BIC <- dimIC(Re, n = n, maxdim = maxdim, Type = "BIC", 
                 graph = graph, verbose = verbose)
    BICVEC[i] <- BIC[which(BIC[,2] == min(BIC[,2])),1]
    
    # Determine optimal dimension according to standard LRT
    LRT <- dimLRT(Re, X = DATAscaled, maxdim = maxdim, 
                  rankDOF = FALSE, graph = graph, 
                  Bartlett = Bartlett, verbose = verbose)
    LRTVEC[i] <- LRT[which(LRT$p.value >= alpha)[1],1]
    
    # Determine 'optimal' dimension according to Guttman bound
    GB       <- dimGB(Re, graph = graph, verbose = verbose)
    GBVEC[i] <- GB[1]
    
    # Determine KMO index in sample
    KMOs <- KMOs + SA(Re)$KMO
    
    # Determine generating communalities
    Comm <- Comm + diag((GD$loadings) %*% t(GD$loadings))
    
    # Can you get some coffee?
    cat(paste("Simulation run ", i, " done\n", sep = ""))
  }
  
  # Return
  return(list(AIC = table(AICVEC), BIC = table(BICVEC), 
              LRT = table(LRTVEC), GB = table(GBVEC),
              averageKMO = KMOs/nsim, communalities = Comm/nsim))
}



#'**---------------------------------**\
#'**Simulations*\
#'**---------------------------------**\

## Set parameters
p       <- c(100,200)
n       <- c(50,75,100,150,250)
m       <- c(5,12,20)
balance <- c(TRUE, FALSE)
cmmnlty <- c(.7,.8,.9)
div5    <- c(40, 20, 15, 15, 10)
div12   <- c(20, 10, 10, 10, 10, 10, 5, 5, 5, 5, 5, 5)
div20   <- c(10, 10, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 3, 3, 3, 3, 3)
Seed    <- 12345

## Simulation loop
for (i in 1:length(p)){
  for (j in 1:length(m)){
    for (k in 1:length(n)){
      for (l in 1:length(balance)){
        for (h in 1:length(cmmnlty)){
          if (m[j] == 5 & p[i] == 100){
            division = div5
          }
          if (m[j] == 5 & p[i] == 200){
            division = 2*div5
          }
          if (m[j] == 12 & p[i] == 100){
            division = div12
          }
          if (m[j] == 12 & p[i] == 200){
            division = 2*div12
          }
          if (m[j] == 20 & p[i] == 100){
            division = div20
          }
          if (m[j] == 20 & p[i] == 200){
            division = 2*div20
          }
          VL  <- sqrt((cmmnlty[h] - .6^2)/(m[j] - 1))
          Run <- DimSim(nsim = 100, 
                        maxdim = 2*m[j],
                        seed = Seed,
                        p = p[i],
                        m = m[j],
                        n = n[k],
                        loadingfix = TRUE,
                        balanced = balance[l],
                        numloadings = division,
                        loadingnegative = TRUE,
                        simplestructure = FALSE,
                        loadingvalue = .6, 
                        loadingvaluelow = VL)
          savename <- paste("p",p[i],".m",m[j],".n",n[k],
                            ".balance",balance[l],".C",cmmnlty[h],
                            sep = "")
          save(Run, file = paste(savename, ".Rdata", sep = ""))
        }
      }
    }
  }
}




#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Functions supporting the Illustrations**
#' **Section 3 of the Main Text**
#' **------------------------------------------------------------------------**

#'
#'**---------------------------------**\
#'**Convenience functions*\
#'**---------------------------------**\

Or2 <- function(this){
  ##############################################################################
  # Convenience function for calculating (overall) explained residual variation 
  # this > crps object from pec
  #
  # Notes:
  # - Function tailored towards the application/comparison at hand.
  #   Can be further generalized.
  ##############################################################################
  
  # Dependencies:
  # require("base")
  # require("pec")
  
  Ar2 <- 1 - this[,1]/this[1,1]
  Cr2 <- 1 - this[,2]/this[1,2]
  Tab <- cbind(Ar2,Cr2)
  colnames(Tab)    <- c("AppR2", "crossvalR2")
  rownames(Tab)[2] <- "FMradio"
  return(Tab)
}


plotR2Table <- function(R2Table, type = "CV", which = NULL, 
                        Ylab = NULL, Xlab = "Time"){
  ##############################################################################
  # Convenience function for plotting explained residual variation 
  # R2Table > R2 object from pec
  # type    > character indicating if the cross-validated or apparent plot
  #           must be returned. Must be one of "CV" or "AE"
  # which   > vector indicating which models to include next to FMradio
  # Ylab    > control over the y-axis label
  #
  # Notes:
  # - Function tailored towards the application/comparison at hand.
  #   Can be further generalized.
  ##############################################################################
  
  # Dependencies:
  # require("base")
  # require("pec")
  
  # Preliminary
  addLegend <- c("Conditional survival forest",
                 "Random survival forest",
                 "Cox boosting",
                 "Cox ridge",
                 "Cox lasso")
  if (is.null(which)){which = c(1,2,3,4,5)}
  
  # Plot
  if (type == "CV"){
    if (!is.null(Ylab)){Ylab = Ylab}
    if (is.null(Ylab)){Ylab = expression(Averaged ~ cross-validated ~ R^2)}
    MIN <- min(R2Table$crossvalErr[,c(-1)])
    MAX <- max(R2Table$crossvalErr[,c(-1)])
    COL <- 2:7
    plot(R2Table$crossvalErr$time, 
         R2Table$crossvalErr$RR.MetaCox, 
         type = "l",
         axes = FALSE,
         xlab = Xlab,
         ylab = Ylab,
         col = COL[1],
         lwd = 1.5,
         ylim = c(MIN,MAX)#,
         #main = "Explained residual variation under cross-validated error")
    )
    axis(2, col = "black", lwd = 1)
    axis(1, col = "black", lwd = 1)
    abline(h = 0)
    if(any(which == 1)){
      lines(R2Table$crossvalErr$time,
            R2Table$crossvalErr$RR.CforestCox,
            type = "l",
            col = COL[2],
            lwd = 1.5)
    }
    if(any(which == 2)){
      lines(R2Table$crossvalErr$time,
            R2Table$crossvalErr$RR.RforestCox,
            type = "l",
            col = COL[3],
            lwd = 1.5)
    }
    if(any(which == 3)){
      lines(R2Table$crossvalErr$time,
            R2Table$crossvalErr$RR.boostCox,
            type = "l",
            col = COL[4],
            lwd = 1.5)
    }
    if(any(which == 4)){
      lines(R2Table$crossvalErr$time,
            R2Table$crossvalErr$RR.ridgeCox,
            type = "l",
            col = COL[5],
            lwd = 1.5)
    }
    if(any(which == 5)){
      lines(R2Table$crossvalErr$time,
            R2Table$crossvalErr$RR.lassoCox,
            type = "l",
            col = COL[6],
            lwd = 1.5)
    }
    legend("bottomright", inset = .03,
           legend = c("FMradio", addLegend[which]),
           col = c(COL[1],COL[which  + 1]), 
           lty = 1,
           cex = .9,
           box.lty = 0,
           lwd = 2,
           y.intersp = 1.7)
  }
  if (type == "AE"){
    if (!is.null(Ylab)){Ylab = Ylab}
    if (is.null(Ylab)){Ylab = expression(Apparent ~ R^2)}
    MIN <- min(R2Table$AppErr[,c(-1)])
    MAX <- max(R2Table$AppErr[,c(-1)])
    COL <- 2:7
    plot(R2Table$AppErr$time, 
         R2Table$AppErr$RR.MetaCox, 
         type = "l",
         axes = FALSE,
         xlab = Xlab,
         ylab = Ylab,
         col = COL[1],
         lwd = 1.5,
         ylim = c(MIN,MAX)#,
         #main = "Explained residual variation under apparent error")
    )
    axis(2, col = "black", lwd = 1)
    axis(1, col = "black", lwd = 1)
    abline(h = 0)
    if(any(which == 1)){
      lines(R2Table$AppErr$time,
            R2Table$AppErr$RR.CforestCox,
            type = "l",
            col = COL[2],
            lwd = 1.5)
    }
    if(any(which == 2)){
      lines(R2Table$AppErr$time,
            R2Table$AppErr$RR.RforestCox,
            type = "l",
            col = COL[3],
            lwd = 1.5)
    }
    if(any(which == 3)){
      lines(R2Table$AppErr$time,
            R2Table$AppErr$RR.boostCox,
            type = "l",
            col = COL[4],
            lwd = 1.5)
    }
    if(any(which == 4)){
      lines(R2Table$AppErr$time,
            R2Table$AppErr$RR.ridgeCox,
            type = "l",
            col = COL[5],
            lwd = 1.5)
    }
    if(any(which == 5)){
      lines(R2Table$AppErr$time,
            R2Table$AppErr$RR.lassoCox,
            type = "l",
            col = COL[6],
            lwd = 1.5)
    }
    legend("bottomright", inset = .03,
           legend = c("FMradio", addLegend[which]),
           col = c(COL[1],COL[which  + 1]), 
           lty = 1,
           cex = .9,
           box.lty = 0,
           lwd = 2,
           y.intersp = 1.7)
  }
}



#'**---------------------------------**\
#'**Extensions pec functionality*\
#'**---------------------------------**\

penL12 <- function(formula, 
                   data, 
                   penalty = "lasso", 
                   fold = 5, 
                   steps = 100, 
                   minlambda = 1, 
                   maxlambda = 10^6,
                   ...){
  ##############################################################################
  # Function that extends the functionality of pec to L1 and L2 regression:
  # Fits L1 or L2 penalized survival model using CV for determining penalty value
  # formula   > formula object representing the model to be fitted
  # data      > A data.frame containing the response and the terms to 
  #             be penalized
  # penalty   > character indicating the type of penalized model to fit.
  #             Must be either "lasso" or "ridge"
  # fold      > numeric scalar indicating the number of folds for cross-validating
  #             either the L1 or L2 penalty on the basis of likelihood 
  #             cross-validation
  # steps     > numeric scalar indicating the maximum number of steps between 
  #             minlambda and maxlambda at which the cross-validated likelihood 
  #             is to be calculated in the global profiling of the likelihood
  # minlambda > numeric scalar indicating the minimum value of the penalty-value 
  # maxlambda > numeric scalar indicating the maximum value of the penalty-value 
  # ...       > optional arguments to be passed to penalized
  # 
  # NOTES:
  # - Uses the penalized package to fit L1 and L2 regressions
  # - Hence, basically a formula interface for the optL1/optL2 and penalized
  #   functions of the penalized package
  # - Assumes that the data object contains standardized features
  # - Assumes that all features on the right-hand side of the formula argument
  #   are subject to penalization
  # - The function first profiles the cross-validated likelihood globally,
  #   after which a local search commences in case of the lasso
  # - Based on penalizedS3 and penalizedOpt (infix functions within pec)
  # - For writing extensions to pec, see:
  #   Mogensen, U.B., Ishawaran, H., & Gerds, T.A. (2012). Evaluating Random 
  #   Forests for Survival Analysis Using Prediction Error Curves. Journal
  #   of Statistical Software, 50(11): 1-23.
  ##############################################################################  
  
  # Dependencies:
  # require("base")
  # require("penalized")
  
  ## Obtain penalized terms
  ff <- as.character(formula)
  response <- formula(paste(ff[[2]],"~1",sep=""))
  terms <- strsplit(ff[[3]],"\\+|\\-")[[1]]
  terms <- gsub("\n","",terms)
  terms <- sapply(terms,function(tt){
    gsub(" ","",tt)
  })
  
  ## Data for penalized features
  penalized <- data[,colnames(data) %in% terms]
  
  ## Call S4 method
  if(penalty == "lasso"){
    seedhere <- round(10^6*runif(1))
    set.seed(seedhere)
    if(maxlambda == 10^6){
      ## Because of potential local optima in case of L1-penalty, first search globally
      p2 <- try(profL1(response = response, penalized = penalized,
                       data = data, fold = fold, steps = steps,
                       minlambda1 = minlambda,...)) 
    } else {
      p2 <- try(profL1(response = response, penalized = penalized,
                       data = data, fold = fold, steps = steps,
                       minlambda1 = minlambda, maxlambda1 = maxlambda,...))
    }
    pen1 <- which.max(p2$cvl)
    if(pen1 == 1) maxpen <- 2*p2$lambda[1] else maxpen <- p2$lambda[pen1-1]
    el <- length(p2$cvl)
    if(pen1 == el) minpen <- p2$lambda[el]/2 else minpen <- p2$lambda[pen1+1]
    
    ## Use the same seed for optL1 as for profL1
    ## And search locally around the global value
    set.seed(seedhere) 
    opt <- try(optL1(response = response, penalized = penalized, 
                     data = data, fold = p2$fold, minlambda1 = max(minlambda,minpen),
                     maxlambda1 = maxpen,...)) 
    lambdaopt <- opt$lambda     
    fitS4 <- penalized(response = response,
                       penalized = penalized,
                       data = data,
                       lambda1 = lambdaopt, 
                       ...)
  } else {
    if(maxlambda == 10^6){
      opt <- try(optL2(response = response, penalized = penalized,
                       data = data, fold = fold, minlambda2 = minlambda,...)) 
    } else {
      opt <- try(optL2(response = response, penalized = penalized,
                       data = data, fold = fold,
                       minlambda2 = minlambda, maxlambda2 = maxlambda,...)) 
    }
    lambdaopt <- opt$lambda                 
    fitS4 <- penalized(response = response,
                       penalized = penalized,
                       data = data,
                       lambda2 = lambdaopt,  
                       ...)
  }
  
  ## Convert to S3 object
  fit <- list(fit = fitS4, call = match.call())
  class(fit) <- "penL12"
  fit
}


predictSurvProb.penL12 <- function(object,
                                   newdata,
                                   times,
                                   ...){
  ##############################################################################
  # Function that extends the functionality of pec to L1 and L2 regression:
  # Extracts survival probabilities
  # object  > fitted R object (of class penL12)
  # newdata > data frame with predictor variables
  # times   > numerical vector containing the times at which the requested 
  #           survival probabilties are evaluated
  # ...     > optional (train) arguments
  #
  # NOTES:
  # - Based on predictSurvProb.penfitS3 (infix function within pec)
  # - For writing extensions to pec, see:
  #   Mogensen, U.B., Ishawaran, H., & Gerds, T.A. (2012). Evaluating Random 
  #   Forests for Survival Analysis Using Prediction Error Curves. Journal
  #   of Statistical Software, 50(11): 1-23.
  ##############################################################################  
  
  # Dependencies:
  # require("prodlim")
  # require("penalized")
  
  penfit <- object$fit
  pCovaNames <- names(penfit@penalized)
  newPen <- newdata[,pCovaNames]
  ptemp <- penalized:::predict(object = penfit,
                               penalized = newPen,
                               data = newdata)
  pos <- sindex(jump.times = ptemp@time, eval.times = times)
  p <- cbind(1,ptemp@curves)[,c(pos+1)]
  p
}




#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Analysis First Illustration: External validation**
#' **Section 3 of the Main Text**
#' **------------------------------------------------------------------------**

#'
#'**---------------------------------**\
#'**Data*\
#'**---------------------------------**\

## Load data
load("C:/Users/cf.peeters/Desktop/Radiomics_pipeline/Data/Radiomics_nonorm_OPSCC_T1_20180317.Rdata")
## OPSCC_VUMC_T1_data: Training data
## OPSCC_UMCU_T1_data: Validation data
## clinical_OPSCC.VUMC: Clinical information training data
## clinical_OPSCC.UMCU: Clinical information validation data



#'**---------------------------------**\
#'**Running Pipeline*\
#'**---------------------------------**\

#########################################
## Step 1: Regularized correlation matrix
#########################################

## Scale data and get raw correlation matrix
DATAscaled <- scale(OPSCC_VUMC_T1_data)
R          <- cor(DATAscaled)

## Redundancy visualization
R0 <- R
R0[abs(R0) < .95] <- 0
radioHeat(R0, diag = FALSE, labelsize = .01)

## Redundancy filtering
## And subsetting data
## 51 features remain
RFs         <- RF(R, t = .95)
DATAscaledS <- subSet(DATAscaled, RFs)

## Optimal penalty
set.seed(303)
OPT <- regcor(DATAscaledS, fold = 5)

## Conditioning
CNplot(cor(DATAscaledS), 
       1.e-7, 1, 
       step = 4000, 
       Iaids = TRUE, 
       target = default.target(cor(DATAscaledS), type = "DUPV"),
       vertical = TRUE, 
       value = OPT$optPen,
       type = "ArchI")

## Obtain regularized correlation matrix
OPT$optPen
Re = OPT$optCor


#########################################
## Step 2: FA Data Compression
#########################################

## Assess dimensionality factor solution
## 10 considered upper bound
## Variance explained would suggest 8 factors
dimGB(Re)
dimVAR(Re, 15, graph = TRUE)

## 9th factor seems weak
## 8th factor seems weak
## Will keep solution at 7
## ML factor analysis with Varimax rotation
fito <- mlFA(Re, 7)
print(fito$Loadings, digits = 2, cutoff = .3, sort = TRUE)

## Visualizing solution
dandpal <- rev(rainbow(100, start = 0.4, end = 0.6))
dandelion(fito$Loadings, bound = .3, mcex = c(1,1), palet = dandpal)


#########################################
## Step 3: Obtaining factor scores
#########################################

## Factor scores
Lambda <- fito$Loadings
Psi    <- fito$Uniqueness
Scores <- facScore(DATAscaledS, Lambda, Psi)

## Determinacy factor scores
## Highly determinate
DF <- facSMC(Re, Lambda); DF



#'**---------------------------------**\
#'**Apparent performance*\
#'**---------------------------------**\

#########################################
## Data
#########################################

## Combine original (scaled) data with projected meta-features
## And include the survival information
DAT <- cbind(DATAscaled, Scores)

## Time to death
Status  <- as.numeric(clinical.OPSCC.VUMC$AorD) - 1
time    <- clinical.OPSCC.VUMC$FUT_OS
DAT     <- cbind(time, Status, DAT)


#########################################
## Comparison of prediction errors
#########################################

## Formulating the model formula's
FitRSF     <- as.formula(paste("Surv(time, Status)~", 
                               paste(colnames(DAT)[3:91], collapse="+")))
FitBoostCox<- as.formula(paste("Hist(time, Status)~", 
                               paste(colnames(DAT)[3:91], collapse="+")))
FitMetaCox <- as.formula(paste("Surv(time, Status) ~", 
                               paste(colnames(DAT[,c(92:98)]), collapse="+")))

set.seed(54321765)
models <- list("MetaCox" = coxph(FitMetaCox, data = DAT, x = TRUE, y = TRUE),
               "CforestCox" = pecCforest(FitRSF, data = DAT, 
                                         controls = cforest_classical(ntree = 1000)),
               "RforestCox" = rfsrc(FitRSF, data = DAT),
               "boostCox" = coxboost(FitBoostCox, data = DAT, cv = TRUE),
               "ridgeCox" = penL12(FitRSF, data = DAT, penalty = "ridge",
                                   minlambda = 1, maxlambda = 10^6,
                                   fold = 5),
               "lassoCox" = penL12(FitRSF, data = DAT, penalty = "lasso",
                                   minlambda = 9, maxlambda = 10^6,
                                   fold = 5, steps = 100))

## Assessing apparent prediction error
## Median follow-up time = 4.37
Mtime = median(time)
PredErrorA <- pec(object = models,
                  formula = Surv(time, Status) ~ 1,
                  data = DAT,
                  exact = TRUE,
                  maxtime = median(time),
                  cens.model = "marginal",
                  verbose = TRUE)

## Summary apparent prediction error
crps(PredErrorA)


#########################################
## Likelihood profiling
#########################################

## Profiling the cross-validated likelihood as a function of the penalty parameter
## Can help understand the performance of Cox ridge and Cox lasso

## Preliminaries
formula = FitRSF
data    = DAT
fold = 5
steps = 1000
minlambda = 1
ff <- as.character(formula)
response <- formula(paste(ff[[2]],"~1",sep=""))
terms <- strsplit(ff[[3]],"\\+|\\-")[[1]]
terms <- gsub("\n","",terms)
terms <- sapply(terms,function(tt){
  gsub(" ","",tt)
})
penalized <- data[,colnames(data) %in% terms]

## Profiling likelihood for ridge
set.seed(282822)
profileL2 <- profL2(response = response, penalized = penalized,
                    data = data, fold = fold, steps = steps,
                    minlambda2 = minlambda, maxlambda2 = 600)

## Profiling likelihood for lasso
set.seed(85487593)
profileL1 <- profL1(response = response, penalized = penalized,
                    data = data, fold = fold, steps = steps,
                    minlambda1 = 9.5, maxlambda1 = 50)

## Visualize
pdf("LikelihoodProfiles.pdf", width = 12, height = 12)
par(mfrow=c(2,1))
  plot(profileL2$lambda, profileL2$cvl, type="l",
       xlab = "Value penalty parameter",
       ylab = "Cross-validated likelihood Cox ridge")
  plot(profileL1$lambda, profileL1$cvl, type="l",
       xlab = "Value penalty parameter",
       ylab = "Cross-validated likelihood Cox lasso")
dev.off()



#'**---------------------------------**\
#'**Validation performance*\
#'**---------------------------------**\

#########################################
## Data
#########################################

## Scale data and get filtered subset
DATAscaledV  <- scale(OPSCC_UMCU_T1_data)
DATAscaledVS <- subSet(DATAscaledV, RFs)

## Obtain factor scores with factor solution from the training data
ScoresV <- facScore(DATAscaledVS, Lambda, Psi)

## Combine original (scaled) data with projected meta-features
## And include the survival information
DATV <- cbind(DATAscaledV, ScoresV)

## Time to death
Status  <- as.numeric(clinical.OPSCC.UMCU$AorD) - 1
clinical.OPSCC.UMCU$FUT_OS <- as.numeric(gsub(",", ".", 
                                              gsub("\\.", "", 
                                                   clinical.OPSCC.UMCU$FUT_OS)))
time    <- clinical.OPSCC.UMCU$FUT_OS
DATV    <- cbind(time, Status, DATV)


#########################################
## Comparison of prediction errors
#########################################

## Median follow-up time = 4.03
PredErrorV <- pec(object = models,
                  formula = Surv(time, Status) ~ 1,
                  data = DATV,
                  exact = TRUE,
                  maxtime = median(time),
                  cens.model = "marginal",
                  verbose = TRUE)

## Summary validation prediction error
crps(PredErrorV)



#'**---------------------------------**\
#'**R2 Results*\
#'**---------------------------------**\

Or2(
  cbind(
    crps(PredErrorA),
    crps(PredErrorV)
  )
)



#'**---------------------------------**\
#'**Visualizing results*\
#'**---------------------------------**\

pdf("External.pdf", width = 15, height = 15)
par(mfrow=c(2,2))
  ## Visualization apparent prediction error
  plot(PredErrorA, what = "AppErr",
       xlab = "Time (years)",
       ylab = "Apparent prediction error",
       legend.cex = .9,
       legend.lty = 1,
       legend.lwd = 2,
       legend.legend = c("Reference model", 
                         "FMradio",
                         "Conditional survival forest",
                         "Random survival forest",
                         "Cox boosting",
                         "Cox ridge",
                         "Cox lasso"),
       add.refline = TRUE,
       lwd = 1.5,
       legend.y.intersp = 1.7)
  
  ## Visualization validation prediction error
  plot(PredErrorV, what = "AppErr",
       xlab = "Time (years)",
       ylab = "Validation prediction error",
       legend.cex = .9,
       legend.lty = 1,
       legend.lwd = 2,
       legend.legend = c("Reference model", 
                         "FMradio",
                         "Conditional survival forest",
                         "Random survival forest",
                         "Cox boosting",
                         "Cox ridge",
                         "Cox lasso"),
       add.refline = TRUE,
       lwd = 1.5,
       legend.y.intersp = 1.7)
  
  
  ## Visualization apparent R2
  R2TableA <- R2(PredErrorA, times = seq(0,Mtime,.01), reference = 1)
  plotR2Table(R2TableA, type = "AE", Xlab = "Time (years)")
  
  ## Visualization Validation R2
  R2TableV <- R2(PredErrorV, times = seq(0,median(time),.01), reference = 1)
  plotR2Table(R2TableV, type = "AE", Ylab = expression(Validation ~ R^2),
              Xlab = "Time (years)")
dev.off()



#'**---------------------------------**\
#'**Recalibration performance*\
#'**---------------------------------**\

#########################################
## On recalibration
#########################################

## Note that the radioFM approach is evaluated more strictly than the other models:
## Not only are the new factor scores based on the training set factor solution.
## Also, these new factor scores are combined with the Cox parameter-estimates 
##   obtained from the training data.
## For the other models the new data are simply combined with the old model-calibrations
## Hence, in a sense, the proposed approach is treated more strictly.

## In determining the new factor scores through the training factor-solution, we
##   are recalibrating the the positioning of the validation obervations w.r.t. 
##   the latent dimensions determined in the training.
## In this respect one could argue that, with shifts in latent-trait positioning,
##   one should recalibrate the Cox estimates.
## Hence, here we perform this exercise.
## Note that the factor-analytic solution is obtained without information on the
##   core outcome measure (time to event).


#########################################
## Performance
#########################################

## Recalibrating the FMradio approach Cox model
## For comparison one could also recalibrate the Cox Boosting model
## One would then fit a Cox PH with the features selected by Cox Boosting
## The performance of the recalibrated estimates is then assessed on the same validation data
## Note that this will give overoptimistic performance, as the outcome measure is
##   used in training and calibrating the Cox Boosting model
## An analoguous argument can be made for recalibrating the Cox lasso model

## Obtaining models for boosting and lasso recalibration
## Formulating the model formula's
These <- which(models$boostCox$coxboost$coefficients[13,] != 0) + 2
fitBoostCoxR <- as.formula(paste("Surv(time, Status) ~", 
                                 paste(colnames(DAT[,These]), collapse="+")))
coefficients(models$lassoCox$fit)
These <- which(colnames(DAT) %in% c("FRACTALDIMENSION_FITTED",
                                    "FRACTALABUNDANCE",
                                    "LGLRE", "SRHGLE", "LRLGLE"))
fitlassoCoxR <- as.formula(paste("Surv(time, Status) ~", 
                                 paste(colnames(DAT[,These]), collapse="+")))

## Model list
set.seed(94643832)
modelsRC <- list("MetaCox" = coxph(FitMetaCox, data = DATV, x = TRUE, y = TRUE),
                 "boostCoxR" = coxph(fitBoostCoxR, data = DATV, x = TRUE, y = TRUE),
                 "lassoCoxR" = coxph(fitlassoCoxR, data = DATV, x = TRUE, y = TRUE))

## Assessing recalibrated prediction error
## Median follow-up time = 4.03
PredErrorVR <- pec(object = modelsRC,
                   formula = Surv(time, Status) ~ 1,
                   data = DATV,
                   exact = TRUE,
                   maxtime = median(time),
                   cens.model = "marginal",
                   verbose = TRUE)

## Summary validation prediction error under recalibration
crps(PredErrorVR)
1 - crps(PredErrorVR)/crps(PredErrorVR)[1]

## Visualizing results
pdf("ExternalRecalibrate.pdf", width = 15, height = 7.5)
par(mfrow=c(1,2))
## Visualization recalibrated prediction error
plot(PredErrorVR, what = "AppErr",
     xlab = "Time (years)",
     ylab = "Recalibrated validation prediction error",
     legend.cex = .9,
     legend.lty = 1,
     legend.lwd = 2,
     legend.legend = c("Reference model", 
                       "FMradio",
                       "Cox boosting",
                       "Cox lasso"),
     add.refline = TRUE,
     lwd = 1.5,
     legend.y.intersp = 1.7,
     col = c(1,2,5,7))

## Visualization Validation R2
R2TableVR <- R2(PredErrorVR, times = seq(0,median(time),.01), reference = 1)
plotR2Table(R2TableVR, type = "AE", 
            which <- c(3,5),
            Ylab = expression(Recalibrated ~ validation ~ R^2),
            Xlab = "Time (years)")
dev.off()


#########################################
## Save results
#########################################

#load(file = "PredErrorExt.Rdata")
save(PredErrorA, PredErrorV, PredErrorVR,
     file = "PredErrorExt.Rdata")




#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Analysis Second Illustration: Internal validation**
#' **Section 3 of the Main Text**
#' **------------------------------------------------------------------------**

#'
#'**---------------------------------**\
#'**Data*\
#'**---------------------------------**\

## Load basic radiomic data (eset object)
## This object removed 4 radiomic features that were de facto constants 
## This object also removed one person with a follow-up time of 0
load("C:/Users/cf.peeters/Desktop/Radiomics_pipeline/Data/ESetRadioMfilter.Rdata")



#'**---------------------------------**\
#'**Running Pipeline*\
#'**---------------------------------**\

#########################################
## Step 1: Regularized correlation matrix
#########################################

## Scale data and get raw correlation matrix
DATAscaled <- scale(t(exprs(ESetRadioMfilter)))
R          <- cor(DATAscaled)

## Redundancy visualization
R0 <- R
R0[abs(R0) < .95] <- 0
radioHeat(R0, diag = FALSE, labelsize = .01)

## Redundancy filtering
## And subsetting data
## 124 features remain
RFs         <- RF(R, t = .95)
DATAscaledS <- subSet(DATAscaled, RFs)

## Optimal penalty
set.seed(303)
OPT <- regcor(DATAscaledS, fold = 5)

## Conditioning
CNplot(cor(DATAscaledS), 
       1.e-7, 1, 
       step = 4000, 
       Iaids = TRUE, 
       target = default.target(cor(DATAscaledS), type = "DUPV"),
       vertical = TRUE, 
       value = OPT$optPen,
       type = "ArchI")

## Obtain regularized correlation matrix
OPT$optPen
Re = OPT$optCor


#########################################
## Step 2: FA Data Compression
#########################################

## Assess dimensionality factor solution
## 13 considered upper bound
## Variance explained would suggest 8 factors
dimGB(Re)
dimVAR(Re, 15, graph = TRUE)

## Assessing solutions around 8
## 9th factor seems weak
## Will keep solution at 8
## ML factor analysis with Varimax rotation
fito <- mlFA(Re, 8)
print(fito$Loadings, digits = 2, cutoff = .3, sort = TRUE)

## Visualizing solution
dandpal <- rev(rainbow(100, start = 0.4, end = 0.6))
dandelion(fito$Loadings, bound = .3, mcex = c(1,1), palet = dandpal)


#########################################
## Step 3: Obtaining factor scores
#########################################

## Factor scores
Lambda <- fito$Loadings
Psi    <- fito$Uniqueness
Scores <- facScore(DATAscaledS, Lambda, Psi)

## Determinacy factor scores
## Highly determinate
DF <- facSMC(Re, Lambda); DF



#'**---------------------------------**\
#'**Compare Methods*\
#'**---------------------------------**\

#########################################
## Data
#########################################

## Combine original (scaled) data with projected meta-features
## And include the survival information
DAT <- cbind(DATAscaled, Scores)

## Time to death
Status  <- as.numeric(ESetRadioMfilter$Death_yesno) - 1
time    <- ESetRadioMfilter$Death_followuptime_months
DAT     <- cbind(time, Status, DAT)


#########################################
## Comparison of prediction errors
#########################################

## Formulating the model formula's
FitRSF     <- as.formula(paste("Surv(time, Status)~", 
                               paste(colnames(DAT)[3:434], collapse="+")))
FitBoostCox<- as.formula(paste("Hist(time, Status)~", 
                               paste(colnames(DAT)[3:434], collapse="+")))
FitMetaCox <- as.formula(paste("Surv(time, Status) ~", 
                               paste(colnames(DAT[,c(435:442)]), collapse="+")))

models <- list("MetaCox" = coxph(FitMetaCox, data = DAT, x = TRUE, y = TRUE),
               "CforestCox" = pecCforest(FitRSF, data = DAT, 
                                         controls = cforest_classical(ntree = 1000)),
               "RforestCox" = rfsrc(FitRSF, data = DAT),
               "boostCox" = coxboost(FitBoostCox, data = DAT, cv = TRUE),
               "ridgeCox" = penL12(FitRSF, data = DAT, penalty = "ridge",
                                   minlambda = 1, maxlambda = 10^5,
                                   fold = 5),
               "lassoCox" = penL12(FitRSF, data = DAT, penalty = "lasso",
                                   minlambda = 30, maxlambda = 10^3,
                                   fold = 5, steps = 100))

## Assessing prediction error
## Median follow-up time = 25.7
set.seed(446464)
PredError <- pec(object = models,
                 formula = Surv(time, Status) ~ 1,
                 data = DAT,
                 exact = TRUE,
                 maxtime = median(time),
                 cens.model = "marginal",
                 splitMethod = "cv5",
                 B = 500,
                 verbose = TRUE)

## Summary results
crps(PredError)
Or2(crps(PredError))


#########################################
## Visualizing results
#########################################

pdf("Internal.pdf", width = 15, height = 15)
par(mfrow=c(2,2))
  ## Visualize apparent prediction error
  plot(PredError, what = "AppErr",
       xlab = "Time (months)",
       ylab = "Apparent prediction error",
       legend.cex = .9,
       legend.lty = 1,
       legend.lwd = 2,
       legend.legend = c("Reference model", 
                         "FMradio",
                         "Conditional survival forest",
                         "Random survival forest",
                         "Cox boosting",
                         "Cox ridge",
                         "Cox lasso"),
       add.refline = TRUE,
       lwd = 1.5,
       legend.y.intersp = 1.7)

  ## Visualize cross-validated prediction error
  plot(PredError, what = "crossvalErr",
       xlab = "Time (months)",
       ylab = "Averaged cross-validated prediction error",
       legend.cex = .9,
       legend.lty = 1,
       legend.lwd = 2,
       legend.legend = c("Reference model", 
                         "FMradio",
                         "Conditional survival forest",
                         "Random survival forest",
                         "Cox boosting",
                         "Cox ridge",
                         "Cox lasso"),
       add.refline = TRUE,
       lwd = 1.5,
       legend.y.intersp = 1.7)
  
  ## Visualize apparent residual explained variation 
  R2Table <- R2(PredError, times = seq(0,median(time),.01), reference = 1)
  plotR2Table(R2Table, "AE", Xlab = "Time (months)")
  
  ## Visualize cross-validated residual explained variation 
  plotR2Table(R2Table, "CV", Xlab = "Time (months)")
dev.off()


#########################################
## Save results
#########################################

#load(file = "PredError.Rdata")
save(PredError, file = "PredError.Rdata")



#'**---------------------------------**\
#'**Compare Methods FILTER*\
#'**---------------------------------**\

#########################################
## Note
#########################################

## - One could argue that the redundancy-filtering is 
##     responsible for success FMradio
## - Hence, below we compare the various methods when 
##    all are preceded by redundancy-filtering


#########################################
## Data
#########################################

## Combine original (scaled) data with projected meta-features
## And include the survival information
DAT <- cbind(DATAscaledS, Scores)

## Time to death
Status  <- as.numeric(ESetRadioMfilter$Death_yesno) - 1
time    <- ESetRadioMfilter$Death_followuptime_months
DAT     <- cbind(time, Status, DAT)


#########################################
## Comparison of prediction errors
#########################################

## Formulating the model formula's
FitRSF     <- as.formula(paste("Surv(time, Status)~", 
                               paste(colnames(DAT)[3:126], collapse="+")))
FitBoostCox<- as.formula(paste("Hist(time, Status)~", 
                               paste(colnames(DAT)[3:126], collapse="+")))
FitMetaCox <- as.formula(paste("Surv(time, Status) ~", 
                               paste(colnames(DAT[,c(127:134)]), collapse="+")))

models <- list("MetaCox" = coxph(FitMetaCox, data = DAT, x = TRUE, y = TRUE),
               "CforestCox" = pecCforest(FitRSF, data = DAT, 
                                         controls = cforest_classical(ntree = 1000)),
               "RforestCox" = rfsrc(FitRSF, data = DAT),
               "boostCox" = coxboost(FitBoostCox, data = DAT, cv = TRUE),
               "ridgeCox" = penL12(FitRSF, data = DAT, penalty = "ridge",
                                   minlambda = 1, maxlambda = 10^5,
                                   fold = 5),
               "lassoCox" = penL12(FitRSF, data = DAT, penalty = "lasso",
                                   minlambda = 5, maxlambda = 10^3,
                                   fold = 5, steps = 100))

## Assessing prediction error
## Median follow-up time = 25.7
set.seed(44646)
PredErrorF <- pec(object = models,
                  formula = Surv(time, Status) ~ 1,
                  data = DAT,
                  exact = TRUE,
                  maxtime = median(time),
                  cens.model = "marginal",
                  splitMethod = "cv5",
                  B = 500,
                  verbose = TRUE)

## Summary results
crps(PredErrorF)
Or2(crps(PredErrorF))


#########################################
## Visualizing results
#########################################

pdf("InternalFilter.pdf", width = 15, height = 15)
par(mfrow=c(2,2))
  ## Visualize apparent prediction error
  plot(PredErrorF, what = "AppErr",
       xlab = "Time (months)",
       ylab = "Apparent prediction error",
       legend.cex = .9,
       legend.lty = 1,
       legend.lwd = 2,
       legend.legend = c("Reference model", 
                         "FMradio",
                         "Conditional survival forest",
                         "Random survival forest",
                         "Cox boosting",
                         "Cox ridge",
                         "Cox lasso"),
       add.refline = TRUE,
       lwd = 1.5,
       legend.y.intersp = 1.7)

  ## Visualize cross-validated prediction error
  plot(PredErrorF, what = "crossvalErr",
       xlab = "Time (months)",
       ylab = "Averaged cross-validated prediction error",
       legend.cex = .9,
       legend.lty = 1,
       legend.lwd = 2,
       legend.legend = c("Reference model", 
                         "FMradio",
                         "Conditional survival forest",
                         "Random survival forest",
                         "Cox boosting",
                         "Cox ridge",
                         "Cox lasso"),
       add.refline = TRUE,
       lwd = 1.5,
       legend.y.intersp = 1.7)
  
  ## Visualize apparent residual explained variation 
  R2TableF <- R2(PredErrorF, times = seq(0,median(time),.01), reference = 1)
  plotR2Table(R2TableF, "AE", Xlab = "Time (months)")
  
  ## Visualize cross-validated residual explained variation 
  plotR2Table(R2TableF, "CV", Xlab = "Time (months)")
dev.off()


#########################################
## Save results
#########################################

#load(file = "PredErrorFilter.Rdata")
save(PredErrorF, file = "PredErrorFilter.Rdata")




#'#############################################################################
#'#############################################################################
#' **------------------------------------------------------------------------**
#' **Second Analysis Second Illustration: Binary Setting**
#' **Section S.5.4. of the Supplement**
#' **------------------------------------------------------------------------**

#########################################
## Note
#########################################

## One referee suggested that a logistic setting might be an insightful
##   addition to the survival setting
## Hence, for illustrative purposes, a logistic setting is offered:
## - The T-stadium is dichotomized
## - T1, T2, T3, and T4 indications refer to size/extent of tumor,
##   with larger numbers after the T referring to larger tumors, 
##   possibly with larger local extensions
## - T-stadium was dichotomized into a score for T1/T2 and a score for T3/T4
## - Then a calibrated AUC is assessed with the logistic counterparts
##   of the models considered in the survival setting
## - These counterparts are: logistic regression with the meta-features
##   as predictors, conditional inference random forest, random forest, 
##   boosted logistic regression, logistic ridge regression, 
##   and logistic lasso regression
## - The ensemble and regularized methods are all evaluated (a) with 
##   the original radiomic features and (b) with the meta-features as
##   input
require("caret")


#########################################
## Data preparation
#########################################

## Binarize T stadium
Status2 <- as.numeric(ESetRadioMfilter$T_stadium)
Status2[Status2 == 1] <- 1
Status2[Status2 == 2] <- 1
Status2[Status2 == 3] <- 2
Status2[Status2 == 4] <- 2
Status2[Status2 == 1] <- "T.stadium12"
Status2[Status2 == 2] <- "T.stadium34"
ThisOne <- which(Status2 == 5)
Status2 <- Status2[-ThisOne]
Status2 <- as.factor(Status2)

## Exclude missing
Scores     <- Scores[-ThisOne,]
DATAscaled <- DATAscaled[-ThisOne,]


#########################################
## Resampling control
#########################################

## Control for repeated cross-validation
ctrl <- trainControl(method = "repeatedcv",
                     number = 5,
                     repeats = 100,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)


#########################################
## Model comparisons
#########################################

## Simple logistic model with meta-features as predictors
set.seed(1635664)
metaGLM <- train(y = Status2,
                 x = Scores,
                 method = "glm",
                 metric = "ROC",
                 trControl = ctrl)


## Tuning parameters for the ensemble and regularized methods:
## RF: Number of Randomly Selected Predictors (mtry, numeric)
## CIRF: Number of Randomly Selected Predictors (mtry, numeric)
## Boosted logistic regression: Number of Boosting Iterations (nIter, numeric)
## GLMnet; ridge: alpha = 0, Regularization Parameter (lambda, numeric)
## GLMnet; lasso: alpha = 1, Regularization Parameter (lambda, numeric)

## Random forest with original radiomic features
set.seed(1635664)
RFA <- train(y = Status2,
             x = DATAscaled,
             method = "rf",
             metric = "ROC",
             tuneLength = 10,
             trControl = ctrl)

## Random forest with meta-features
set.seed(1635664)
RFF <- train(y = Status2,
             x = Scores,
             method = "rf",
             metric = "ROC",
             tuneLength = 10,
             trControl = ctrl)


## Conditional inference random forest with original radiomic features
set.seed(1635664)
CRFA <- train(y = Status2,
              x = DATAscaled,
              method = "cforest",
              metric = "ROC",
              tuneLength = 10,
              trControl = ctrl)

## Conditional inference random forest with meta-features
set.seed(1635664)
CRFF <- train(y = Status2,
              x = Scores,
              method = "cforest",
              metric = "ROC",
              tuneLength = 10,
              trControl = ctrl)


## Boosted logistic regression with original radiomic features
set.seed(1635664)
BLRA <- train(y = Status2,
              x = DATAscaled,
              method = "LogitBoost",
              metric = "ROC",
              tuneLength = 10,
              trControl = ctrl)

## Boosted logistic regression with meta-features
set.seed(1635664)
BLRF <- train(y = Status2,
              x = Scores,
              method = "LogitBoost",
              metric = "ROC",
              tuneLength = 10,
              trControl = ctrl)


## Logistic ridge regression with original radiomic features
set.seed(1635664)
lambdas <- seq(1, 10000, length = 100)
LRRA <- train(y = Status2,
              x = DATAscaled,
              method = "glmnet",
              metric = "ROC",
              tuneGrid = data.frame(alpha = 0, lambda = lambdas),
              trControl = ctrl)

## Logistic ridge regression with meta-features
set.seed(1635664)
lambdas <- seq(1e-7, 2, length = 100)
LRRF <- train(y = Status2,
              x = Scores,
              method = "glmnet",
              metric = "ROC",
              tuneGrid = data.frame(alpha = 0, lambda = lambdas),
              trControl = ctrl)


## Logistic lasso regression with original radiomic features
set.seed(1635664)
lambdas <- seq(1, 10000, length = 100)
LLRA <- train(y = Status2,
              x = DATAscaled,
              method = "glmnet",
              metric = "ROC",
              tuneGrid = data.frame(alpha = 1, lambda = lambdas),
              trControl = ctrl)

## Logistic lasso regression with meta-features
set.seed(1635664)
lambdas <- seq(1e-7, 2, length = 100)
LLRF <- train(y = Status2,
              x = Scores,
              method = "glmnet",
              metric = "ROC",
              tuneGrid = data.frame(alpha = 1, lambda = lambdas),
              trControl = ctrl)


#########################################
## Save
#########################################

#load("LogisticExample.Rdata")
save(metaGLM, RFA, RFF, CRFA, CRFF, 
     BLRA, BLRF, LRRA, LRRF, LLRA, LLRF,
     file = "LogisticExample.Rdata")

