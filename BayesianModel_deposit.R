library(rjags)

filePrefix <- "/Volumes/UNIXBACKUP/modelInput/"
dataFile <- "Classified_LQTS13_Dec_6.txt"
genePriorFile <- "GenePriors.txt"
plotFile <- "BayesPredictor_LQTS_13_Dec_caseFreqSub.png"
permutDir <- "PredictionSets_LQTS_13_Dec12"
modelFile <- "Ruklisa_AdditionalData2_HierModelLQTSFull.bug"
resultsFile <- "PredRes_LQTS_13_Dec_caseFreqSub.txt"
convergFile <- "PredConv_LQTS_13_Dec_caseFreqSub.txt"
itCount <- 5000
plotHeader <- "13 LQTS genes, full model"
modelType <- 12
isBrugada <- 0
isHCM <- 0

computeCoefficients <- 0
plotPosterior <- 0
coefPlotFile <- "PosteriorDistr_LQTS"
coefFile <- "EffectSizes_LQTS.txt"
itCountCoef <- 40000

classifiedVariants <- read.table(file=paste(filePrefix, dataFile, sep=""), sep = "\t", header=T)

knotNo <- 3
pphDim <- 1

# Column names for rare variant table
pphColumn <- "pph2_prob"
siftColumn <- "SIFT.score"
grantColumn <- "Grantham"
freqCol <- "X1000Genomes.Freq"
addFreqCol <- "ESP.EurAm.freq"
scoreColumn <- "Cons.Score"

modelColumns <- c(pphColumn, siftColumn)
coreMatr <- classifiedVariants[classifiedVariants[,"reliableClasses"] != 0, modelColumns]
coreMatr[, siftColumn] <- -coreMatr[, siftColumn] + 1

# Grantham scores
coreMatr <- cbind(coreMatr, classifiedVariants[, grantColumn])
dimnames(coreMatr)[[2]][dim(coreMatr)[[2]]] <- grantColumn
coreMatr[, grantColumn] <- replace(coreMatr[, grantColumn], is.na(coreMatr[, grantColumn]) & as.numeric(classifiedVariants[classifiedVariants[,"reliableClasses"] != 0, "Variant.Effect"]) - 1 == 1, 0)
coreMatr[, grantColumn] <- coreMatr[, grantColumn] / 205

# Add conservation score
coreMatr <- cbind(coreMatr, replace(replace(replace(classifiedVariants[, scoreColumn], classifiedVariants[, scoreColumn] == 1, 0), classifiedVariants[, scoreColumn] %in% c(2, 3, 4), 1), classifiedVariants[, scoreColumn] == 5, 2) / 2)
dimnames(coreMatr)[[2]][dim(coreMatr)[[2]]] <- scoreColumn

scorePattern <- c(1:dim(classifiedVariants)[[1]])
scorePattern[] <- 0
Zconserv <- cbind(as.matrix(replace(scorePattern, classifiedVariants[, scoreColumn] == 2, 1)), replace(scorePattern, classifiedVariants[, scoreColumn] %in% c(3, 4), 1), replace(scorePattern, classifiedVariants[, scoreColumn] == 5, 1))

# Add frequency
hasFreq <- -(1 - (is.na(as.numeric(classifiedVariants[, freqCol])) & is.na(as.numeric(classifiedVariants[, addFreqCol]))))

# Separate column to indicate inframe indels
radicalVars <- as.character(classifiedVariants[,"Effect"])
posIndex <- c(1:dim(classifiedVariants)[[1]])
radicalVars2 <- as.numeric(factor(replace(radicalVars, posIndex[!(radicalVars %in% c("inframe"))], "base"))) - 1
coreMatr <- cbind(coreMatr, radicalVars2)
dimnames(coreMatr)[[2]][dim(coreMatr)[[2]]] <- "is.inframe"

# New classification of radical variants, without inframe indels
coreMatr <- cbind(coreMatr, replace(as.numeric(classifiedVariants[, "Variant.Effect"]) - 1, coreMatr[, "is.inframe"] == 1, 0))
dimnames(coreMatr)[[2]][dim(coreMatr)[[2]]] <- "is.new.radical"

# Replace missing SIFT and PPH2 values
coreMatr[, pphColumn] <- replace(coreMatr[, pphColumn], is.na(coreMatr[, pphColumn]) & (coreMatr[, "is.new.radical"] == 1 | coreMatr[, "is.inframe"] == 1), 0)
coreMatr[, siftColumn] <- replace(coreMatr[, siftColumn], is.na(coreMatr[, siftColumn]) & (coreMatr[, "is.new.radical"] == 1 | coreMatr[, "is.inframe"] == 1), 0)
coreMatr[, scoreColumn] <- replace(coreMatr[, scoreColumn], is.na(coreMatr[, scoreColumn]), 0)

Zconserv <- replace(Zconserv, is.na(Zconserv), 0)

toExcluded <- apply(coreMatr, 1, function(x)any(is.na(x)))

# gene x radical classes
geneCount <- length(levels(factor(classifiedVariants[, "Gene"])))
geneRad <- as.numeric(factor(paste(coreMatr[toExcluded == FALSE, "is.new.radical"], classifiedVariants[toExcluded == FALSE, "Gene"])))
geneRadNames <- levels(factor(paste(coreMatr[toExcluded == FALSE, "is.new.radical"], classifiedVariants[toExcluded == FALSE, "Gene"])))

# Collect names of various domains and check if a domain effect has to be taken into account
domainGroups <- paste(classifiedVariants[toExcluded == FALSE, "Gene"], classifiedVariants[toExcluded == FALSE, "Region"], classifiedVariants[toExcluded == FALSE, "Domain"])

domainPresent <- c(1:length(domainGroups))
domainPresent[] <- 1
domainPresent[c(grep("KCNJ2", domainGroups), grep("CAV3", domainGroups), grep("SCN4B", domainGroups), grep("AKAP9", domainGroups), grep("SNTA1", domainGroups), grep("KCNJ5", domainGroups), grep("CACNA2D1", domainGroups), grep("CACNB2", domainGroups), grep("GPD1L", domainGroups), grep("SCN1B", domainGroups))] <- 0
domainPresent[c(grep("MYBPC3", domainGroups), grep("TNNT2", domainGroups))] <- 0
domainPresent[ coreMatr[, "is.new.radical"] == 1 ] <- 0

domainGroupsIndex <- as.numeric(factor(domainGroups[domainPresent == 1]))
fullDomainIndex <- c(1:length(domainGroups))
fullDomainIndex[] <- 1
fullDomainIndex[ domainPresent == 1 ] <- domainGroupsIndex
domainGroupsIndex <- fullDomainIndex

# New gene priors
geneOrder <- levels(factor(paste(coreMatr[toExcluded == FALSE, "is.new.radical"], classifiedVariants[toExcluded == FALSE, "Gene"])))

priors <- read.table(file=paste(filePrefix, genePriorFile, sep=""), sep = "\t", header=T)

if (isHCM) {
  priorValues <- priors[, "Literature"] * 100
	
  # Reorder gene priors
  priorCoefGene <- c(1:length(geneOrder))
  for (i in 1:dim(priors)[[1]]) {
    priorPos <- grep(priors[i, "Gene"], geneOrder)
    priorCoefGene[ priorPos[1] ] <- priorValues[i]
    if (length(priorPos) > 1) {
      priorCoefGene[ priorPos[2] ] <- priorValues[i]
    }
  }
  geneStore <- log(priorCoefGene)
} else {
  # Take priors from prospective cohort
  caseFreq <- c(1:dim(priors)[[1]])
  for (i in 1:dim(priors)[[1]]) {
    if (priors[i, "Case.freq"] == 0) {
      caseFreq[i] <- priors[i, "Literature"] * 0.5
    } else {
      caseFreq[i] <- priors[i, "Familion.2500"]
    }
  }
  priorValues <- (caseFreq - priors[, "ESP"]) / priors[, "ESP"]
  priorValues <- replace(priorValues, priorValues < 0.2, 0.2)

  # Radical gene priors
  radFreq <- c(1:dim(priors)[[1]])
  for (i in 1:dim(priors)[[1]]) {
    if (priors[i, "Case.freq.rad"] > 0) {
      radFreq[i] <- priors[i, "Case.freq.rad"]
    } else {
      radFreq[i] <- max(priors[i, "ESP.rad"], 0.0002)
    }
  }

  radPrior <- (radFreq - replace(priors[, "ESP.rad"], priors[, "ESP.rad"] == 0, 0.0002)) / replace(priors[, "ESP.rad"], priors[, "ESP.rad"] == 0, 0.0002)
  radPrior <- replace(radPrior, radPrior < 0.2, 0.2)

  # Reorder gene priors
  priorCoefGene <- c(1:length(geneOrder))
  for (i in 1:dim(priors)[[1]]) {
    priorPos <- grep(priors[i, "Gene"], geneOrder)
    priorCoefGene[ priorPos[1] ] <- priorValues[i]
    if (length(priorPos) > 1) {
      priorCoefGene[ priorPos[2] ] <- radPrior[i]
    }
  }
  geneStore <- log(priorCoefGene)
}

# Assign pathogenicity outcome
yBin <- classifiedVariants[classifiedVariants[,"reliableClasses"] != 0, "reliableClasses"]
yBin <- (replace(yBin, yBin == -1, 0))

if (isBrugada == 1) {
  # Priors for Brugada model
  priorCoefGene[] <- 1
  priorCoefGene[ grep("SCN5A", geneOrder) ] <- 30
  geneStore <- log(priorCoefGene)
}

# Initialize model parameters
if (modelType == 12) {
  pred.inits <- function () {
    list (mu.gene=rnorm(1, 0, 10), mu.rad=rnorm(1, 0, 10), b.pph=rnorm(pphDim, 0, 10), b.sift=rnorm(pphDim, 0, 10), b.grantham=rnorm(pphDim, 0, 10), prior.scale=rnorm(1, 0, 10), prior.rad.scale=rnorm(1, 0, 10), a.domain=rnorm(length(levels(factor(domainGroups[ domainPresent == 1 ]))), 0, 10), sigma.d=rgamma(1, 0.5, 5), b.inframe=rnorm(1, 0, 10), b.hasfreq=rnorm(1, 0, 10), b.cons=rnorm(2, 0, 10), aph=rnorm(1, 0, 10), agr=rnorm(1, 0, 10))
  }
}

if (computeCoefficients == 0) {
  # Cross validation, where pathogenicity of variants from a test set is predicted.
  png(paste(filePrefix, "Results/", plotFile, sep=""), height=550, width=800)
  par(omi=c(0,0,0,0))
  par(mar=c(3.5,4,2,0))  

  plot(1 ~ 1, ylim=c(0, 1), xlim=c(0, 101), type="n", main="", xlab="", ylab="", cex.axis=1.2)
  mtext("probability of pathogenicity", side=2, cex=1.4, line=2.5)
  mtext("permutation no.", side=1, cex=1.4, line=2)
  mtext(plotHeader, cex=1.6, line=0.5)

  for (h in 1:100) {
    print(paste("Test sample", h))

    selectedVars <- read.table(file = paste(filePrefix, permutDir, "/V_", h, ".txt", sep=""), colClasses="character")
    modelSet <- apply(selectedVars, 1, function(x)as.numeric(substring(x[1], 1:nchar(x[1]), 1:nchar(x[1]))))

    yRepl <- replace(yBin, modelSet[, 1] == 0, NA)
    yTested <- c(1:length(yRepl))[modelSet[, 1] == 0]

    pred.data <- list(y=yRepl, n=dim(coreMatr)[[1]], k=geneCount, r=(max(geneRad) - geneCount), e=length(levels(factor(domainGroups[ domainPresent == 1 ]))), pph=coreMatr[, pphColumn], sift=coreMatr[, siftColumn], grantham=coreMatr[, grantColumn], geneXrad=geneRad, inframe=coreMatr[, "is.inframe"], domain=domainGroupsIndex, genePrior=geneStore, is.radical=coreMatr[, "is.new.radical"], hasfreq=hasFreq, cons.primates=(Zconserv[, 1] + Zconserv[, 2] + Zconserv[, 3]), cons.all=Zconserv[, 3], domain.present=domainPresent, testedIndex=yTested, vv=length(yTested))
  
    jagsObject <- jags.model(paste(filePrefix, modelFile, sep=""), data=pred.data, inits=pred.inits, n.chains=10, n.adapt=1500)

    jagsRes <- jags.samples(jagsObject, c("y.tilde"), n.iter=itCount)

    selectedLimits <- c((itCount / 5 + 1):itCount)

    predTest <- apply(jagsRes$y.tilde[, selectedLimits, ], 1, mean)

    # Check convergence
    R.tilde <- c(1:dim(jagsRes$y.tilde)[[1]])
    for (i in 1:dim(jagsRes$y.tilde)[[1]]) {
      wMean <- apply(jagsRes$y.tilde[i, selectedLimits, ], 2, mean)
      totMean <- mean(wMean)
      B <- 0
      W <- 0
      for (k in 1:dim(jagsRes$y.tilde)[[3]]) {
        B <- B + (wMean[k] - totMean)^2
        W <- W + sum((jagsRes$y.tilde[i, selectedLimits, k] - wMean[k])^2) / (length(selectedLimits) - 1)
      }
      B <- B * length(selectedLimits) / (dim(jagsRes$y.tilde)[[3]] - 1)
      W <- W / dim(jagsRes$y.tilde)[[3]]

      varFi <- (length(selectedLimits) - 1) * W / length(selectedLimits) + B / length(selectedLimits)
      R.tilde[i] <- sqrt(varFi / W)
    }
 
    yTrue <- classifiedVariants[modelSet[, 1] == 0, "reliableClasses"]
    yTrue <- replace(yTrue, yTrue == -1, 0)

    points(predTest[yTrue == 1] ~ rep(h, length(predTest[yTrue == 1])), pch=20, col="magenta", bg="magenta")
    points(predTest[yTrue == 0] ~ rep(h, length(predTest[yTrue == 0])), pch=20, col="green", bg="green")

    if (h == 1) {
      resultsMatrix <- as.matrix(predTest)
      convMatrix <- R.tilde
    } else {
      resultsMatrix <- cbind(resultsMatrix, predTest)
      convMatrix <- cbind(convMatrix, R.tilde)
    }
  }
  dev.off()

  write.table(resultsMatrix, file=paste(filePrefix, "Results/", resultsFile, sep=""), sep = "\t", row.names=F, col.names=F, quote=F)
  write.table(convMatrix, file=paste(filePrefix, "Results/", convergFile, sep=""), sep = "\t", row.names=F, col.names=F, quote=F)
} else {
  # Compute effect sizes for all model terms
  parList <- c("sigma.d", "prior.scale", "prior.rad.scale", "a.gene", "a.domain", "b.pph", "b.sift", "b.grantham", "mu.gene", "mu.rad", "b.inframe", "b.hasfreq", "b.cons", "aph", "agr", "y.tilde")
  obsParList <- c("sigma.d", "prior.scale", "prior.rad.scale", "a.gene", "a.domain", "b.pph", "b.sift", "b.grantham", "mu.gene", "mu.rad", "b.inframe", "b.hasfreq", "b.cons", "aph", "agr")

  intPars <- obsParList
  totalDim <- c()
  for (i in 1:length(intPars)) {
    currDim <- intPars[i]
    if (intPars[i] == "a.gene") {
      currDim <- geneRadNames
    }
    if (intPars[i] == "a.domain") {
      currDim <- levels(factor(domainGroups[domainPresent == 1]))
    }
    if (intPars[i] == "b.cons") {
      currDim <- paste(intPars[i], c(1:2))
    }
    totalDim <- c(totalDim, currDim)
  }
  	
  yRepl <- yBin
  yTested <- c(1:length(yRepl))

  pred.data <- list(y=yRepl, n=dim(coreMatr)[[1]], k=geneCount, r=(max(geneRad) - geneCount), e=length(levels(factor(domainGroups[ domainPresent == 1 ]))), pph=coreMatr[, pphColumn], sift=coreMatr[, siftColumn], grantham=coreMatr[, grantColumn], geneXrad=geneRad, inframe=coreMatr[, "is.inframe"], domain=domainGroupsIndex, genePrior=geneStore, is.radical=coreMatr[, "is.new.radical"], hasfreq=hasFreq, cons.primates=(Zconserv[, 1] + Zconserv[, 2] + Zconserv[, 3]), cons.all=Zconserv[, 3], domain.present=domainPresent, testedIndex=yTested, vv=length(yTested))

  jagsObject <- jags.model(paste(filePrefix, modelFile, sep=""), data=pred.data, inits=pred.inits, n.chains=10, n.adapt=3000)

  jagsRes <- jags.samples(jagsObject, obsParList, n.iter=itCountCoef)

  selectedLimits <- c((itCountCoef / 5 + 1):itCountCoef)

  # Extracting of coefficients and checking of parameter convergence
  intPars <- obsParList
  currMeans <- c()
  currMedians <- c()
  R.tilde <- c(1:length(totalDim))
  R.tilde[] <- 0
  k <- 1
  coefArray <- array(0, dim=c(length(totalDim), length(selectedLimits), 10))
  dimnames(coefArray)[[1]] <- totalDim
  for (i in 1:length(intPars)) {
    chainDim <- dim(jagsRes[[intPars[i]]])[[3]]

    if (length(dim((jagsRes[[intPars[i]]])[, selectedLimits, ])) == 2) {
      currBlock <- (jagsRes[[intPars[i]]])[, selectedLimits, ]
      parMeans <- mean(currBlock)
      parMedians <- median(currBlock)

      wMean <- apply(currBlock, 2, mean)
      totMean <- mean(wMean)
      B <- 0
      W <- 0
      for (j in 1:chainDim) {
        B <- B + (wMean[j] - totMean)^2
        W <- W + sum((currBlock[, j] - wMean[j])^2) / (length(selectedLimits) - 1)
      }
      B <- B * length(selectedLimits) / (chainDim - 1)
      W <- W / chainDim

      varFi <- (length(selectedLimits) - 1) * W / length(selectedLimits) + B / length(selectedLimits)
      coefArray[k, , ] <- (jagsRes[[intPars[i]]])[, selectedLimits, ]
      R.tilde[k] <- sqrt(varFi / W)
      k <- k + 1
    } else {
      currBlock <- (jagsRes[[intPars[i]]])[, selectedLimits, ]
      parMeans <- apply(currBlock, 1, mean)
      parMedians <- apply(currBlock, 1, median)
    
      for (h in 1:dim(jagsRes[[intPars[i]]])[[1]]) {
        smallBlock <- currBlock[h, , ]
        wMean <- apply(smallBlock, 2, mean)
        totMean <- mean(wMean)
        B <- 0
        W <- 0
        for (j in 1:chainDim) {
          B <- B + (wMean[j] - totMean)^2
          W <- W + sum((smallBlock[, j] - wMean[j])^2) / (length(selectedLimits) - 1)
        }
        B <- B * length(selectedLimits) / (chainDim - 1)
        W <- W / chainDim
      
        varFi <- (length(selectedLimits) - 1) * W / length(selectedLimits) + B / length(selectedLimits)

        coefArray[k, , ] <- (jagsRes[[intPars[i]]])[h, selectedLimits, ]
        R.tilde[k] <- sqrt(varFi / W)
        k <- k + 1
      }
    }
    currMeans <- c(currMeans, parMeans) 
    currMedians <- c(currMedians, parMedians)
  }
  coefArray[1, , ] <- 1 / coefArray[1, , ]
  
  currMeans[1] <- 1 / currMeans[1]
  currMedians[1] <- 1 / currMedians[1]

  allMeans <- matrix(currMeans, ncol=1, byrow=F)
  allMedians <- matrix(currMedians, ncol=1, byrow=F)
  dimnames(allMeans)[[1]] <- totalDim
  dimnames(allMedians)[[1]] <- totalDim
  
  parSummary <- cbind(allMeans, allMedians)
  dimnames(parSummary)[[2]] <- c("coef.mean", "coef.median")
  write.table(as.data.frame(parSummary), file=paste(filePrefix, "Results/", coefFile, sep=""), sep = "\t", row.names=T, col.names=T, quote=F)

  if (plotPosterior == 1) {
    # Plot posterior distributions of all model parameters
    pageSize <- 10
    widthGr <- 5
    lengthGr <- 2
    rounds <- ceiling(dim(coefArray)[[1]] / pageSize)

    for (h in 1:rounds) {
      blockStart <- pageSize * (h - 1)
      limit <- min(blockStart + pageSize, dim(coefArray)[[1]])
      firstIm <- blockStart + 1

      png(paste(filePrefix, "Results/", coefPlotFile, "_", h, ".png", sep=""), height=550, width=1000)
      par(omi=c(0.05,0.05,0.4,0))
      par(mar=c(4,4,6,1))
      par(mfrow=c(lengthGr, widthGr))

      limit <- limit - blockStart
      for (i in 1:limit) { 
        rNo <- ceiling(i / widthGr)
        cNo <- (i - 1) %% widthGr + 1
        par(mfg=c(rNo, cNo))
        j <- i + blockStart

        if (h == 1 && i <= 1) {
          hist(as.numeric(coefArray[j, , ]), col="black", main="", ylab="", xlab="", axes=F, xlim=c(0, 40), breaks=c(c(0:40), max(coefArray[1:2, , ])))
        } else {
          hist(as.numeric(coefArray[j, , ]), col="black", main="", ylab="", xlab="", axes=F, xlim=c(-20, 20))
        }

        abline(v=c(currMeans[j]), col="magenta")
        abline(v=c(currMedians[j]), col="red")

        axis(2, labels=T, tick=T, cex.axis=1.3)
        axis(1, labels=T, tick=T, cex.axis=1.3)

        options(digits=4)
        mtext(paste(format(currMedians[j])), side=1, cex=1.4, line=3, col="dark blue")

        headerText <- dimnames(coefArray)[[1]][j]
        headerSplit <- substring(headerText, 1:nchar(headerText), 1:nchar(headerText))
        if (length(headerSplit) > 20) {
          if (length(headerSplit) <= 40) {
            headerText <- paste(paste(headerSplit[1:20], collapse=""), "\n ", paste(headerSplit[21:length(headerSplit)], collapse=""), sep="")
          } else {
            headerText <- paste(paste(headerSplit[1:20], collapse=""), "\n ", paste(headerSplit[21:40], collapse=""), "\n ", paste(headerSplit[41:length(headerSplit)], collapse=""), sep="")
          }
        }
    
        mtext(headerText, cex=1.2, line=0.5)
      }
      mtext(paste(plotHeader, ", page", h), outer=TRUE, cex=1.5, line=0.5)
      dev.off()
    }
  }	
}
