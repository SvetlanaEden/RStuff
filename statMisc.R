####################### author: Svetlana K. Eden, svetlana.eden@vanderbilt.edu
################################################
###
###
### various stat gadgets 

allEffectsForCRR = function(crrModel, varNames, varRanges, listOfKnots=NULL, splineSearchString = "'"){
  ### !!! ATTENTION: effects are not implemented when interactions are present.
  ###### Example: allEffectsForCRR(crrModel = compRiskModelScrAgeTrial, varNames=c("baselineAge", "cr"), varRanges=list(baselineAge = c(50, 51), cr=matrix(c(0, 1, 1, 5), byrow=TRUE, ncol=2)), listOfKnots = list(cr=quantile(bootData$cr, probs=c(.1, .5, .9))))
  ## returns the effects in the format of Hazard Ratio (HR)
  ### varNames - a vector of variables (spline components are considered as one variable)
  ###   if a variable was model with splines the function will need knots for each variable
  ### effectRange - a vector of length 2: from - to
  ### knots - if the variable is a spline, knots must be supplied to calculate the effect
  ### 
  ### effectForCRRByGrayAndFine(crrModel, varName="cr", effectRange=c(0, 1), knots=quantile(bootData$cr, probs=c(.1, .5, .9)))
  if (!all(varNames %in% names(varRanges))) stop("Provide varRanges for all variables in varNames")
  res = data.frame(Variable = c(), From = c(), To = c(), Effect = c(), Lower = c(), Upper = c(), EffectPvalue = c(), OverallPvalue = c())
  for (v in varNames){
    #cat("+++", v, "\n")
    if (is.null(dim(varRanges[[v]]))){
      res = rbind(res, effectForCRRByGrayAndFine(crrModel, varName = v, effectRange = varRanges[[v]], varKnots=listOfKnots[[v]]), splineSearchString = splineSearchString)
    }else{
      for (i in 1:nrow(varRanges[[v]])){
        res = rbind(res, effectForCRRByGrayAndFine(crrModel, varName = v, effectRange = varRanges[[v]][i,], varKnots=listOfKnots[[v]], splineSearchString = splineSearchString))
      }
    }
  }
  res
}

effectForCRRByGrayAndFine = function(crrModel, varName, effectRange, varKnots=NULL, splineSearchString = "'"){
  ####### Example: effectForCRRByGrayAndFine(crrModel = compRiskModelScrAgeTrial, varName="cr", effectRange=c(0, 1), varKnots=quantile(bootData$cr, probs=c(.1, .5, .9)))
  ## returns the effect in the format of Hazard Ratio (HR)
  ### varName - a name of a covariate
  ###   if the covariate was model with splines the function will look for name', name'', and so on in order to build splines for effectRange
  ### effectRange - a vector of length 2: from - to
  ### varKnots- if the variable is a spline, knots must be supplied to calculate the effect
  ### 
  ### effectForCRRByGrayAndFine(crrModel, varName="cr", effectRange=c(0, 1), varKnots=quantile(bootData$cr, probs=c(.1, .5, .9)))
  waldPvalueFromGray <- function (x, kk){  ## multiparameter Wald test
    w <- x$coef[kk]
    v <- x$var[kk, kk]
    stat <- sum(w * solve(v, w)) ## solves the equation ‘v %*% x = w’ for ‘x’
    df <- length(w)
    c(stat = stat, df = df, p = 1 - pchisq(stat, df))
  }
  varNames = names(crrModel[["coef"]])
  colnames(crrModel[["var"]]) = varNames
  rownames(crrModel[["var"]]) = varNames
  splineComponents = grep(paste("^", varName, "[",splineSearchString ,"]*", sep=""), varNames, value=TRUE)
  if (length(splineComponents) == 0) stop("Variable ", varName, " was not found in the model -  cannot provide the effect.")
  if (length(splineComponents) == 1){
    rangeDiff = diff(effectRange)
  }else{
    if (is.null(varKnots)) stop("Variable ", varName, " is a spline. To calculate the effect, components must be provided.")
    splineComponents = splineComponents[order(splineComponents)]
    varMatrixRange = components(effectRange, knots = varKnots)
    rangeDiff = matrix(apply(varMatrixRange, 2, diff), nrow=1)
  }
  betas = crrModel[["coef"]][splineComponents]
  vars = crrModel[["var"]][splineComponents, splineComponents]
  totalSE = as.numeric(sqrt(rangeDiff %*% vars %*% t(rangeDiff)))
  effect = sum(rangeDiff * betas)
  CIs = effect + c(-1, +1)*1.96*totalSE
  pval = waldPvalueFromGray(crrModel, splineComponents)["p"]
  effectPval = 1 - pchisq((effect/totalSE)^2, 1)
  #effectPval = 2*pnorm(abs(effect/totalSE)*(-1))
  res = data.frame(varName, effectRange[1], effectRange[2], exp(effect), exp(CIs[1]), exp(CIs[2]), effectPval, pval)
  names(res) = c("Variable", "From", "To", "Effect", "Lower", "Upper", "EffectPvalue", "OverallPvalue")
  res
}

contrastsFromCPH = function(model, varRangesForContrast, theRestOfVarList = NULL){
  ### returns effect in the format of Hazard Ratio (HR)
  ###### example call contrastsFromCPH(cphmodel, varRangesForContrast=list(age=matrix(c(10, 0, 40, 50), byrow=TRUE, ncol=2),
  ######                               cr = c(0,1), stam = c("less or equal to zero", "greater than zero")))
  ### varRangesForContrast: matrix
  ###     for example, if we look at the contrast for CESD, and we want two HRs:
  ###          one 0 relatively to 16 and the other
  ###              16 relatively to 55
  ### we have to supply a matrix   0  16
  ###                              16 55
  ### and the function will return two contrasts.
  ### theRestOfVarList: values of the rest of the variables to use in contrast
  ###     for example if there is an interaction with race and literacy
  ###     one can put theRestOfVarList = list(race = "White", literacy = "Low")
  ###     note that each variable should have only one value unlike in varRangesForContrast
  eachVarTimes = sapply(varRangesForContrast, function(e){if(is.null(dim(e))){1}else{dim(e)[1]}})
  varNames = unlist(sapply(names(eachVarTimes), function(n)rep(n, eachVarTimes[n])))
  resDataFrame = data.frame(Variable = varNames, From=NA, To=NA, Effect = NA, Lower = NA, Upper = NA, Pvalue = NA, AnovaPvalue = NA, stringsAsFactors = FALSE)
  index = 1
  anovaPvalues = anova(model)
  for (v in names(varRangesForContrast)){
    varRange = matrix(t(varRangesForContrast[[v]]), byrow=TRUE, ncol=2)
    for (j in 1:dim(varRange)[1]){
      argSubList1 = list(); argSubList2 = list()
      resDataFrame[index, "From"] = varRange[j,1]
      resDataFrame[index, "To"] = varRange[j,2]
      argSubList1[[v]] = resDataFrame[index, "To"]
      argSubList2[[v]] = resDataFrame[index, "From"]
      if (!is.null(theRestOfVarList)){
        theRestOfVarListUpdated = theRestOfVarList[setdiff(names(theRestOfVarList), v)]
        for (n in names(theRestOfVarListUpdated)){
          argSubList1[[n]] = theRestOfVarList[[n]]
          argSubList2[[n]] = theRestOfVarList[[n]]
        }
      }
      ## the following check is put just in case unlist() doesn't return the names in the same order as in varRangesForContrast
      if (v != resDataFrame$Variable[index]) stop("The data set is not in the same order as the list of ranges")
      argList = list(model, argSubList1, argSubList2)
      contrastsum = do.call(contrast, argList)
      contrres = c(exp(c(contrastsum[["Contrast"]], contrastsum[["Lower"]], contrastsum[["Upper"]])), contrastsum[["Pvalue"]])
      names(contrres) = c("Effect", "Lower", "Upper", "Pvalue")
      resDataFrame$Effect[index] = contrres["Effect"]
      resDataFrame$Lower[index] = contrres["Lower"]
      resDataFrame$Upper[index] = contrres["Upper"]
      resDataFrame$Pvalue[index] = contrres["Pvalue"]
      pvalueIndex = grep(paste("^", resDataFrame$Variable[index], "[ ]*\\(", sep="" ), rownames(anovaPvalues))
      if (length(pvalueIndex) == 0) pvalueIndex = grep(paste("^", resDataFrame$Variable[index], "$", sep="" ), rownames(anovaPvalues))
      resDataFrame$AnovaPvalue[index] = anovaPvalues[pvalueIndex, "P"]
      index = index + 1
    }
  }
  for (n in c("Effect", "Lower", "Upper")){
    resDataFrame[[n]] = round(resDataFrame[[n]], 3)
  }
  resDataFrame$Pvalue = pvalueStr(resDataFrame$Pvalue)
  resDataFrame
}



### effective sample size for logistic and proportional odds logistic regression
### can be also used for linear regression, but for linear regression
### computing samples size is trivial
reportEffectiveSampleSizeForPOLR = function(fundata, outcomes, covariates){
  for (i in 1:length(outcomes)){
    ### exclude missing rows:
    tmp = apply(fundata[, c(outcomes[i], covariates)], 1, function(x)all(!is.na(x)))
    modeldata = fundata[tmp,]
    nRecords = nrow(modeldata)
    ## table(data$outc)nrow(data)  -  proportion of each category in the data
    numOfCategoriesInOutcome = table(modeldata[outcomes[i]])
    ## effective sample size: n - (1/n^2)*sum(n1^3 + n2^3 +...+ nk^3), where k is #of categ.
    if (length(numOfCategoriesInOutcome) == 2){
      effectiveSampleSize = min(numOfCategoriesInOutcome)
    }else{
      effectiveSampleSize = nRecords-(1/nRecords^2)*sum(numOfCategoriesInOutcome^3)
    }
    degreesOfFreedom = effectiveSampleSize/15
    cat("\nOtucome:", outcomes[i], label(modeldata[[outcomes[i]]]), "\n")
    print(table(modeldata[[outcomes[i]]]))
    cat("\tNumber of non-missing records:", nRecords, "\n")
    cat("\tEffective sample size:", round(effectiveSampleSize), "\n")
    cat("\tNumber of variables (binary or continuous linear) to include without overfitting:", round(degreesOfFreedom, 1), "\n\n")
  }
}

hosmerlem = function (logRegModel, g = 10){
### SAS reference:
### http://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_logistic_sect039.htm
  y = logRegModel$y
  yhat = predict(logRegModel)
  cutyhat <- cut(yhat, breaks = quantile(yhat, probs = seq(0, 1, 1/g)), include.lowest = T)
  obs <- xtabs(cbind(1 - y, y) ~ cutyhat)
  expect <- xtabs(cbind(1 - yhat, yhat) ~ cutyhat)
  chisq <- sum((obs - expect)^2/expect)
  P <- 1 - pchisq(chisq, g - 2)
  c("X^2" = chisq, Df = g - 2, "P(>Chi)" = P)
}


spearmanAndCIs = function(x1, x2, confLevel=0.95){
  ### calculates Spearman Correlation and its confidence intervals based on Fisher_transformation
  ### see for references:
  ### http://www.stat.psu.edu/online/courses/stat509/18_agree/18_agree_print.htm
  ### http://en.wikipedia.org/wiki/Fisher_transformation
  ### !!! Note the first website has an error in the formula for CIs:
  ### instead of t distribution with n-3 degrees of freedom, there should be normal distribution
  ### Example 1 is from the website: http://www.stat.psu.edu/online/courses/stat509/18_agree/18_agree_print.htm
  ### x1 = c(23,23,27,27,39,41,45,49,50,53,53,54,56,57,58,58,60,61)
  ### x2 = c(9.5, 27.9, 7.8, 17.8, 31.4, 25.9, 27.4, 25.2, 31.1, 34.7, 42.0, 29.1, 32.5, 30.3, 33.0, 33.8, 41.1, 34.5)
  ### Example 2 is from the website: http://www.statsdirect.com/help/nonparametric_methods/spear.htm
  ### x1 = c(4,10,3,1,9,2,6,7,8,5)
  ### x2 = c(5,8,6,2,10,3,9,4,7,1)
  if (length(x1)!=length(x2))stop("Lengths of x1 and x2 must be the same")
  rankx1 = rank(x1)
  rankx2 = rank(x2)
  n = length(x1)
  nThing = n*((n+1)/2)^2
  rho = (sum(rankx1*rankx2) - nThing)/(sqrt(sum(rankx1*rankx1) - nThing) * sqrt(sum(rankx2*rankx2) - nThing))
  Zp = .5 * log((1+rho)/(1-rho))
  standardError = 1/sqrt(n-3)
  CIsTrans = Zp + c(qnorm((1-confLevel)/2), qnorm((1+confLevel)/2))*standardError
  CIs = (exp(2*CIsTrans)-1)/(exp(2*CIsTrans)+1)
  res = c(rho, CIs)
  names(res) = c("rho", "Lower", "Upper")
  res
}
 
 
pvalueStr = function(pvalue, roundnum=3){
### --- gives a p-value in as "<.001"
  str = rep("", length(pvalue))
  for (i in 1:length(pvalue)){
    if(!is.na(pvalue[i])){
      if (pvalue[i] < 0.0001){
        str[i] = "<0.0001"
      }else{
        if (pvalue[i] < 0.001){
          str[i] = "<0.001"
        }else{
          str[i] = as.character(round(pvalue[i], roundnum))
        }
      }
    }else{
      str[i] = ""
    }
  }
  str
}


bootMedianConfInt<-function (x, conflevel = 0.95, bootnum = 1000, narm = TRUE, seed=99) {
### used by medianAndIQR
  sampleFun = function(vector, n=length(vector)){
    sample(vector, n, replace=TRUE)
  }
  if (narm){
    x <- x[!is.na(x)]
  }
  n <- length(x)
  xmedian <- median(x)
  if (n < 2){
    return (c(Median= xmedian, Lower = NA, Upper = NA))
  }
  set.seed(seed)
  bootmedians <- sapply(1:bootnum, FUN=function(iboot, vector, n) { median(sampleFun(vector, n), na.rm = narm)}, vector=x, n=length(x))
  quant <- quantile(bootmedians, c((1 - conflevel)/2, (1 + conflevel)/2))
  names(confint) <- NULL
  res <- c(Median = xmedian, Lower = quant[1], Upper = quant[2])
  res
}

medianAndIQR = function(data, varNames, forValue=NULL){
### the function returns medians, IQRs (inter-quantile range) and number of non-missing for each variable
### in 'varNames' in 'data' (data.frame)
### or if 'varNames' contains only one name and 'forValue' is not NULL (wide format)
### then it returns means and CIs and number of non-missing for variable='varNames' 
### for each value of 'forValue' in 'data' (long format)
  mainErrorMessage = "meanAndCI: Choose long or wide format: either varNames has more than one variable name and forValue is null OR varNames has only one variable name and forValue contain another variable name."
  if (length(varNames)>1 & !is.null(forValue) | length(varNames)==1 & is.null(forValue))    stop(mainErrorMessage)
  if (length(varNames)==0)    stop("meanAndCI: 'varNames' is empty !")
  if (!is.null(forValue))    if (length(forValue)!=1)      stop("meanAndCI: The length of 'forValue' is not one !")
  if (length(setdiff(c(varNames, forValue), names(data)) != 0))    stop("meanAndCI: 'varNames' or 'forValue' are not in the 'data'.")
  na.rm=TRUE
  if (length(varNames)>1){
    res = matrix(NA, nrow=length(varNames), ncol=4)
    for (i in 1:length(varNames)){
      n = varNames[i]
      obsnum = sum(!is.na(data[[n]]))
      #confLevel=0.95
      funres = bootMedianConfInt(data[[n]], narm = TRUE)
      res[i, 1] = funres[1]
      res[i, 2:3] = funres[2:3]
      res[i, 4] = obsnum
    }
    colnames(res) = c("median", "iqrlow", "iqrhigh", "nonmissingnum")
    rownames(res) = varNames
  }else{
    values = setdiff(data[[forValue]], NA)
    values = values[order(values)]
    res = matrix(NA, nrow=length(values), ncol=4)
    for (i in 1:length(values)){
      v = values[i]
      d = data[[varNames]][[data[[forValue]]==v]]
      obsnum = sum(!is.na(d))
      #res[i, 1] = median(d, na.rm=na.rm)
      funres = bootMedianConfInt(d, narm = TRUE)
      res[i, 1] = funres[1]
      res[i, 2:3] = funres[2:3]
      res[i, 4] = obsnum
    }
    colnames(res) = c("median", "iqrlow", "iqrhigh", "nonmissingnum")
    rownames(res) = values
  }
  res
}

modelOutput = function(model, roundn=3, convertToChar=TRUE){
  ## takes the essential output from the model
  ## model - cox regression without interations
  ## roundn - number of decimals to round up the output
  ## convertToChar - if TRUE, the character matrix is returned
  tmps = summary(model)
  tmpa = anova(model)
  varnames = rownames(tmps)[seq(1,nrow(tmps), 2)]
  outmatr = matrix(tmps[seq(2,nrow(tmps), 2), c(1:4,6,7)], nrow = round(nrow(tmps)/2))
  rownames(outmatr) = varnames
  shortrownames = sapply(strsplit(rownames(outmatr), " "), function(x){x[[1]]})
  outmatr = cbind(outmatr, tmpa[shortrownames, 3])
  colnames(outmatr) = c("Low", "High", "Diff.", "HR", "95%HR.Lo", "95%HR.Up", "P-value")
  res = outmatr
  if (convertToChar){
    ##-------------------convert to character
    outstr = matrix(as.character(round(outmatr, roundn)), nrow = nrow(outmatr))
    outstr[,7] = pvalueStr(outmatr[,7], roundnum = roundn)
    colnames(outstr) = colnames(outmatr)
    rownames(outstr) = rownames(outmatr)
    res = outstr
    for (i in 1:ncol(outstr)){
      #cat(decimalAlign(res[,i]), "\n")
      res[,i] = decimalAlign(res[,i])
    }
  }
  res
}

printModelOutput = function(modelOutputStr){
  ncolchar = pmax(nchar(colnames(modelOutputStr)), apply(modelOutputStr, 2, function(x)max(nchar(x))))
  colnames(modelOutputStr) = fillWith(colnames(modelOutputStr), n = ncolchar)
  rownames(modelOutputStr) = fillWith(rownames(modelOutputStr), n = max(nchar(rownames(modelOutputStr))), align="l")
  for (i in 1:ncol(modelOutputStr)){
    modelOutputStr[, i] = fillWith(modelOutputStr[, i], n=ncolchar[i])
  }
  cat(fillWith("", nchar(rownames(modelOutputStr))[1]), paste(colnames(modelOutputStr), collapse="  "), "\n")
  for (i in 1:nrow(modelOutputStr)){
    cat(rownames(modelOutputStr)[i], paste(modelOutputStr[i,], collapse="  "), "\n")
  }
}
#printModelOutput(modelOutput(f.nomiss))
#modelOutputStr = modelOutput(f.nomiss)


decimalAlign = function(number){
  ## Centers numbers around a decimal point by adding spaces to the left and to the right of the number
  ## number can be a vector
  splitnumber = sapply(number, function(s)length(strsplit(s, "\\.")[[1]]))
  if (max(splitnumber)>2) stop("check format of 'number' and try again (too many 'periods')")
  if (max(splitnumber)<2){ #whole number
    res = fillWith(number, max(nchar(number)))
  }else{
    numberStr = as.character(number)
    numberStr[is.na(number)] = ""
    spl = strsplit(numberStr, "\\.")
    wholen = sapply(spl, function(x){x[1]})
    decimn = sapply(spl, function(x){x[2]})
    wholen[is.na(wholen)] = ""
    decimn[is.na(decimn)] = ""
    maxwholen = max(nchar(wholen))
    maxdecimn = max(nchar(decimn))
    wholenf = sapply(wholen, function(x){paste(paste(rep(" ", maxwholen-nchar(x)), collapse=""), x, sep="")})
    decimnf = sapply(decimn, function(x){paste(x, paste(rep(" ", maxdecimn-nchar(x)), collapse=""), sep="")})
    res = paste(wholenf, decimnf, sep=".")
    res = gsub("\\. |\\.$", "  ", res)
    #res = gsub("\\. ", "  ", res)
    res = gsub("^[ ]+$", NA, res)
  }
  res
}

fillWith = function(str, n, fill=" ", align="r"){
## Fills a string with character 'fill' either to the right ('align'=l)
## or to the left ('align'=l) of the string.
## Fills up to 'n' characters. If 'str' is longer than n - dose nothing.
## str and n can be vectors
  if (length(n)==1) n = rep(n, length(str))
  if (align=="r"){
    sapply(1:length(str), function(i, x, n){paste(paste(rep(fill, max(n[i]-nchar(x[i]), 0, na.rm=TRUE)), collapse=""), x[i], sep="")}, x=str, n=n)
  }else{
    sapply(1:length(str), function(i, x, n){paste(x[i], paste(rep(fill, max(n[i]-nchar(x[i]), 0, na.rm=TRUE)), collapse=""), sep="")}, x=str, n=n)
  }
}

plotSomethingAndItsInterval = function(sthgX, sthgY, lowerlimit, upperlimit, delta=min(abs(upperlimit-lowerlimit))/15, ...){
  lines(sthgX, sthgY, ...)
  segments(x0=sthgX, y0=lowerlimit, x1=sthgX, y1=upperlimit, ...)
  segments(x0=sthgX-delta/2, y0=lowerlimit, x1=sthgX+delta/2, y1=lowerlimit, ...)
  segments(x0=sthgX-delta/2, y0=upperlimit, x1=sthgX+delta/2, y1=upperlimit, ...)
}

percent = function(x, roundn=2){
  #round(100*x, signif=roundn)
  signif(100*x, digits=roundn)
}

subsuperhistForCont = function(contVar, varlabel="", textround=2, textstepY=1, cex.axis=.8, showXRange = NULL, showYRange = NULL, vertical=TRUE, binSize=NULL, roundNum = 0, printCounts=TRUE, superImpose = FALSE){
  xticknum = 10
  varTickNum = 9
  var = categorizeContinuousX(contVar, binSize=binSize, roundNum = roundNum)
  counts = table(var, exclude=FALSE)
  names(counts)[names(counts)=="NA"] = "Missing"
  #names(counts)[names(counts)=="NaN"] = "NaN"
  x = as.vector(counts/sum(counts))
  mindif = 1
  y = 1:length(x)
  names(y) = names(counts)
  ydelta = mindif/2.2
  textstartX = 2*max(x)
  textstartY = max(y)
  textstepY = textstepY*diff(range(y))/4
  if (is.null(showXRange)){
    showXRange = c(0, 1.05*max(x))
  }
  axisXLabelAt = min(showXRange) + (0:(xticknum-1))*(diff(showXRange))/(xticknum-1)
  axisXLabels = round(axisXLabelAt, textround)
  if (is.null(showYRange)){
    showYRangeActual = c(min(y)-ydelta, max(y)+ydelta)
    a = 0; k = 1
    axisYLabelAt = y
    axisYLabels = as.character(round(as.numeric(names(y)), textround))
    axisYLabels[is.na(axisYLabels)] = "Missing"
  }else{
    actualVarValues = as.numeric(names(y[setdiff(names(y), c("Missing", "NaN"))]))
    minVal = min(actualVarValues)
    maxVal = max(actualVarValues)
    minY = min(y[setdiff(names(y), c("Missing", "NaN"))])
    maxY = max(y[setdiff(names(y), c("Missing", "NaN"))])
    if (minY - maxY != 0){
      k = (minVal - maxVal)/(minY - maxY)
    }else{
      if (binSize == 0){
        binSize = min(actualVarValues, na.rm=TRUE)*.1
      }
      k = binSize/(textstepY*4)
    }
    a = minVal - k*minY
    showYRangeActual = (range(showYRange) - a)/k
    varTicks = min(showYRange) + (0:(varTickNum-1))*(diff(showYRange))/(varTickNum-1)
    axisYLabelAt = (varTicks-a)/k
    axisYLabels = round(varTicks, textround)
  }
  if (!vertical){
    if (!superImpose){
      plot(x=showXRange, y=showYRangeActual, type="n", axes=FALSE, ylab="variable values", xlab="proportion")
      text(x=textstartX, y=textstartY, labels=varlabel, pos=2, cex=cex.axis)
      axis(side=1, at = axisXLabelAt, labels = axisXLabels, las=2, cex.axis=cex.axis)
      axis(side=2, at = axisYLabelAt, line=0, labels = axisYLabels, las=2, cex.axis=cex.axis)
    }else{
      rect(xleft=0, ybottom=y-ydelta, xright=x, ytop=y+ydelta, col="gray", border=NA)
    }
    if (printCounts){
      text(x=x, y=y, labels=counts, pos=4, cex=cex.axis, col="black")
    }
  }else{
    if (!superImpose){
      plot(y=showXRange, x=showYRangeActual, type="n", axes=FALSE, xlab="variable values", ylab="proportion")
      text(y=textstartX, x=textstartY, labels=varlabel, pos=2, cex=cex.axis)
      axis(side=2, at = axisXLabelAt, labels = axisXLabels, las=2, cex.axis=cex.axis)
      axis(side=1, at = axisYLabelAt, line=0, labels = axisYLabels, las=2, cex.axis=cex.axis)
    }else{
      rect(ybottom=0, xleft=y-ydelta, ytop=x, xright=y+ydelta, col="gray", border=NA)
    }
    if (printCounts){
      text(x=y, y=x, labels=counts, pos=3, cex=cex.axis, col="black")
    }
  }
}

categorizeContinuousX = function(x, binSize=NULL, roundNum = 3){
  xRange = range(x, na.rm=TRUE)
  if (is.null(binSize)){
    binSize = diff(xRange)/10
  }
  if (binSize == 0){
    binSize = min(x, na.rm=TRUE)*.1
  }
  numOfBins = ceiling(diff(xRange)/binSize + .5)
  binStartOffset = (numOfBins*binSize - diff(xRange))/2
  lowerBinEdge = min(xRange) - binStartOffset + ((1:numOfBins)-1) * binSize
  upperBinEdge = min(xRange) - binStartOffset + (1:numOfBins) * binSize
  if (FALSE){
    plot(c(x, xRange + c(- binStartOffset, binStartOffset)), type="n")
    points(x, pch=19)
    text(1:length(x), x, round(x, 4), cex=.7, pos=3)
    abline(h=lowerBinEdge, col="red")
    abline(h=upperBinEdge, col="green")
  }
  greater = sapply(x, function(x){sum(lowerBinEdge <= x)})
  less = sapply(x, function(x){length(upperBinEdge) - sum(x < upperBinEdge) + 1})
  levelNames = (lowerBinEdge + upperBinEdge)/2
  xLevels = (lowerBinEdge[greater] + upperBinEdge[less])/2
  xLevels[is.na(x)] = "NA"
  levelNames = c(levelNames, "NA")
  xLevels[is.nan(x)] = "NaN"
  if (any(xLevels=="NaN")){
    levelNames = c(levelNames, "NaN")
  }
  res = factor(xLevels, levels = levelNames)
  res
}

subsuperhist = function(var, vartype, varlabel=label(var), textround=2, textstepY=1, cex.axis=.8, showYRange = NULL){
  nachar = function(x){x=as.character(x);for (i in 1:length(x)){if (is.na(x[i])&(!is.nan(x[i]))) x[i]="<NA>"}; x }
  naType=c(NA, NaN)
  names(naType) = as.character(naType)
  xticknum = 10
  counts = table(var, exclude=NULL)
  names(counts)[is.na(names(counts))] = "Missing"
  x = as.vector(counts/sum(counts))
  if (vartype=="ordinal"){
    if ( length(setdiff(var, c(NA, NaN)))<2 ){
      mindif = 1
    }else{
      mindif = abs( min( diff( sort( unique(var) ) ), na.rm=TRUE) )
    }
    y = sort(unique(var), na.last=TRUE)
    if (length(counts) != length(y)){
      y = c(y, NA)
    }
    y[is.nan(y)] = max(y, na.rm=TRUE)+mindif
    y[is.na(y)&(!is.nan(y))] = max(y, na.rm=TRUE)+mindif*(any(is.nan(y))+1)
  }
  if (vartype=="factor"){
    mindif = 1
    y = 1:length(x)
  }
  if (any(vartype %nin% c("factor", "ordinal"))){stop("Wrong data type")}
  names(y) = names(counts)

  ydelta = mindif/2.2
  textstartX = 2*max(x)
  textstartY = max(y)
  textstepY = textstepY*diff(range(y))/4
  if (is.null(showYRange)){
    showYRange = c(min(y)-ydelta, max(y)+ydelta)
  }
  plot(x=c(0, 2*max(x)), y=showYRange, type="n", axes=FALSE, ylab="variable values", xlab="proportion")
  text(x=textstartX, y=textstartY, labels=varlabel, pos=2, cex=cex.axis)
  rect(xleft=0, ybottom=y-ydelta, xright=x, ytop=y+ydelta, col="gray", border=NA)
  axis(side=1, at = (0:(xticknum-1))*max(x)/(xticknum-1), labels = round((0:(xticknum-1))*max(x)/(xticknum-1), textround), las=2, cex.axis=cex.axis)
  axis(side=2, at = y, line=-1, labels = names(y), las=2, cex.axis=cex.axis)
  text(x=x, y=y, labels=counts, pos=4, cex=cex.axis, col="black")
  if (vartype=="ordinal"){
    text(x=textstartX, y=textstartY-textstepY, labels=paste("median (IQR) = ", round(median(var, na.rm=TRUE), textround), " (", paste(round(quantile(var, probs=c(0.25, 0.75), na.rm=TRUE), textround), collapse=", "), ") ", sep=""), pos=2, cex=cex.axis)
    text(x=textstartX, y=textstartY-2*textstepY, labels=paste("mean=", signif(mean(var, na.rm=TRUE), digits=textround), "   SD=", signif(sd(var, na.rm=TRUE), digits=textround), sep=""), pos=2, cex=cex.axis)
    abline(h=quantile(var, probs = c(.5), na.rm = TRUE), col="gray", lty=3)
  }
}

superhist = function(data, vars, continVars, binSizeForContinVars = rep(NA, length(continVars)), graphnum=length(vars), textstepY=1, cex.axis=.8, showYRange=NULL){
  #source("/home/edensk/myStuff/mySources/statMisc.R")
  ### varTypes = varType defined for each variable
  ### binSizeForCont = bin size to convert cont into categorical for the histogram.
  ### Example:
  if (FALSE){
    set.seed(100)
    n = 59
    data = data.frame(cont = rnorm(n), ordin=sample(1:6, n, replace=TRUE), alsoord = round(rnorm(n), 1), col=factor(sample(c("red", "blue"), n, replace=TRUE)))
    data$ordin[21:25] = NA
    data$cont[18:28] = NaN
    data$cont[29:35] = NA
    data$ordin[1:10] = NaN
    data$ordin[data$ordin == 5] = NA
    varTypes = c("continuous", "ordinal", "ordinal", "factor")
    names(varTypes) = c("cont", "ordin", "alsoord", "col")
    binSizeForCont = c(.2)
    names(binSizeForCont) = c("cont")
    superhist(data, vars = c("cont", "ordin", "alsoord", "col"), varTypes = varTypes, binSizeForCont=binSizeForCont)
  }
  if(length(binSizeForContinVars) != length(continVars)){
    stop("Length of argument", binSizeForContinVars, "should be the same as length of argument", continVars)
  }
  par(oma=c(0, 0, 0, 0), mar=c(3,3,0,0), mfrow=c(graphnum, 1))
  for(n in vars){
    varClass = class(data[[n]])
    if(varClass=="character"){
      stop("Convert variable", n, "into a factor")
    }
    if(n %in% continVars){
      vartype = "continuous"
    }else{
      if (varClass == "factor"){
        vartype = "factor"
      }else{
        vartype = "ordinal"
      }
    }
    if (vartype == "continuous"){
      if (is.na(binSizeForCont[n])) binSize = NULL
      data[[n]] = categorizeContinuousX(data[[n]], binSize = binSize)
      vartype = "factor"
    }
    subsuperhist(data[[n]], vartype=vartype, textstepY=textstepY, cex.axis=cex.axis, showYRange=showYRange)
  }
}

twoAndTotal <- function(var1, var2, positVal, negatVal, var1Label=NULL, var2Label=NULL){
###
### functions builds two-by-two and totals table for sencitivity and specificity
###
### var1: true test (gold standard), contains only two unique values (positive and negative)
### var2: some other test, contains only two unique values (positive and negative)
### positVal: how positive value is coded
### negatVal: how negative value is coded
### example:
### var1 = c(1, 1, 0, NA, 1, 1, 0, 0, 0, NA, NA, NA, 0)
### var2 = c(1, 0, 0, NA, NA, 0, 0, 1, NA, 1, 0, 0, 0)
### twoAndTotal(var1, var2, positVal=1, negatVal=0)
  if (length(var1)!=length(var2)){stop("Lengths of var1 and var2 have to be the same")}
  if (length(setdiff(unique(var1), NA))>2){stop("var1 must have two unique values")}
  if (length(setdiff(unique(var2), NA))>2){stop("var2 must have two unique values")}
  if (FALSE){
    if (class(positVal)=="character"){
      var1 = as.character(var1)
      var2 = as.character(var2)
    }else{
      var1 = as.numeric(var1)
      var2 = as.numeric(var2)
    }
    var1 = factor(var1, levels=c(positVal, negatVal))
    var2 = factor(var1, levels=c(positVal, negatVal))
  }
  onlyNonMissing = !is.na(var1) & !is.na(var2)
  nmvar1 = var1[onlyNonMissing]
  nmvar2 = var2[onlyNonMissing]
  res = matrix(rep(NA, 4), nrow=2)
  res[1,1] = sum(nmvar1[nmvar2==positVal]==positVal)
  res[2,1] = sum(nmvar1[nmvar2==positVal]==negatVal)
  res[1,2] = sum(nmvar1[nmvar2==negatVal]==positVal)
  res[2,2] = sum(nmvar1[nmvar2==negatVal]==negatVal)
  #res = as.matrix(table(var1[onlyNonMissing], var2[onlyNonMissing]))
  trueStr = var1Label
  testStr = var2Label
  if (is.null(trueStr)) {trueStr="true"}
  if (is.null(testStr)) {testStr="test"}
  #res = res[,2:1]
  #res = res[2:1,]
  res = cbind(res, c(NA,NA))
  res = rbind(res, c(NA,NA,NA))
  res[,3] = res[,1] + res[,2]
  res[3,] = res[1,] + res[2,]
  rownames(res) = c(positVal, negatVal, "total")
  colnames(res) = c(positVal, negatVal, "total")
  rownames(res) = paste(c(trueStr), rownames(res))
  colnames(res) = paste(c(testStr), colnames(res))
  res
}

sensspecppvnpv <- function(var1, var2, positVal, negatVal, var1Label=NULL, var2Label=NULL, wilson=TRUE){
### calculates sensitivity and specificity (and not implemented: their confidence intervals)
### requires Hmisc for wilson conf.int.
### var1: true test (gold standard), contains only two unique values (positive and negative)
### var2: some other test, contains only two unique values (positive and negative)
### positVal: how positive value is coded
### negatVal: how negative value is coded
### returns matrix of sens., spec.,( and not implemented conf. int.) with corresponding col. and row names
### if wilson == TRUE calculates wilson confidence intervals
### otherwise - Wald CIs
  onlyRelevantValues = (var1 %in% c(positVal, negatVal)) & (var2 %in% c(positVal, negatVal))
  var1 = var1[onlyRelevantValues]
  var2 = var2[onlyRelevantValues]
  sensspec = twoAndTotal(var1, var2, positVal=positVal, negatVal=negatVal, var1Label=var1Label, var2Label=var2Label)
  sens = sensspec[1,1]/sensspec[1,3]
  spec = sensspec[2,2]/sensspec[2,3]
  ppv = sensspec[1,1]/sensspec[3,1]
  npv = sensspec[2,2]/sensspec[3,2]
  alpha = .05
  zhalfalpha = 1.96
  ### Confidence Intervals:
  ciStr = if (wilson) "Wilson CI" else "Wald CI"
  p = sens
  n = sensspec[1,3]
  ##
  sensci = if (wilson) wilsonCIsForProportion(p, n, confLevel=.95) else p + c(-1,1)*zhalfalpha*sqrt(p*(1-p))/sqrt(n)
  p = spec
  n = sensspec[2,3]
  #specci = p + c(-1,1)*zhalfalpha*sqrt(p*(1-p))/sqrt(n)
  specci = if (wilson) wilsonCIsForProportion(p, n, confLevel=.95) else p + c(-1,1)*zhalfalpha*sqrt(p*(1-p))/sqrt(n)

  p = ppv
  n = sensspec[3,1]
  #ppvci = p + c(-1,1)*zhalfalpha*sqrt(p*(1-p))/sqrt(n)
  ppvci = if (wilson) wilsonCIsForProportion(p, n, confLevel=.95) else p + c(-1,1)*zhalfalpha*sqrt(p*(1-p))/sqrt(n)
  p = npv
  n = sensspec[3,2]
  #npvci = p + c(-1,1)*zhalfalpha*sqrt(p*(1-p))/sqrt(n)
  npvci = if (wilson) wilsonCIsForProportion(p, n, confLevel=.95) else p + c(-1,1)*zhalfalpha*sqrt(p*(1-p))/sqrt(n)

  res = list()
  res$table = sensspec
  res$res = matrix(c(sens, spec, ppv, npv, sensci[1], specci[1], ppvci[1], npvci[1], sensci[2], specci[2], ppvci[2], npvci[2]), nrow=4, ncol=3)
  colnames(res$res) = c("", paste(ciStr, c("lower", "upper")))
  rownames(res$res) = c("Sens", "Spec", "PPV", "NPV")
  res
}

kappa = function(matr){
  ## matr is a matrix M with matr[1,1]=A=(test1pos|test2pos), matr[2,2]=D=(test1neg|test2neg), matr[2,1]=C=(test1pos|test2neg), matr[1,2]=B=(test1neg|test2pos)
  ## matr:
  ##              test1_Pos   test1_Neg
  ##  test2_Pos     A           B
  ##  test2_Neg     C           D
  ##
  if(nrow(matr)!=ncol(matr)) stop("matrix has to be square")
  mtrace = 0
  for(i in 1:nrow(matr)){mtrace = mtrace + matr[i,i]}
  sumofall = sum(matr)
  sumofrow = apply(matr, 1, sum)
  sumofcol = apply(matr, 2, sum)
  Aobs = mtrace/sumofall
  #Aexp = ((A+B)*(A+C) + (C+D)*(B+D))/(A+C+B+D)^2
  Aexp = (sum(sumofrow*sumofcol))/(sumofall)^2
  if (sum(mtrace)==sum(matr)){
    kappa=1
  }else{
    kappa = (Aobs - Aexp)/(1-Aexp)
  }
  kappa
}

#matrix(1 - c(c(0:3), c(1,0,1,2), c(2,1,0,1), c(3,2,1,0))/(4-1), nrow=4, byrow=TRUE)
#  ### quaWeightMatr = matrix(1 - (c(0:3, c(1,0,1,2), c(2,1,0,1), c(3,2,1,0)))^2/(4-1)^2, nrow=4, byrow=TRUE)
kappaGeneral = function(matr, weightMatr=NULL){
  ## matr is a matrix M with matr[1,1]=A=(test1pos|test2pos), matr[2,2]=D=(test1neg|test2neg), matr[2,1]=C=(test1pos|test2neg), matr[1,2]=B=(test1neg|test2pos)
  ### ref: http://support.sas.com/documentation/cdl/en/procstat/63032/HTML/default/procstat_freq_a0000000647.htm
  ### ref for this kappa's assimptotic variance: See Fleiss, Cohen, and Everitt (1969) for details. 
  ### example
  ### m = matrix(c(76, 0, 7, 10, 2, 0, 3, 1, 17, 1, 15, 8, 20, 3, 5, 11), nrow=4, byrow=TRUE)
  ### # for regular kappa
  ### kappaGeneral(m)
  ### linWeightMatr = matrix(NA, nrow=nrow(matr), ncol=ncol(matr))
  ### for (i in 1:nrow(matr)) for (j in 1:ncol(matr)) linWeightMatr[i,j] = 1-(abs(i-j))/(nrow(matr)-1)
  ### kappaGeneral(m, linWeightMatr)
  ### quaWeightMatr = matrix(NA, nrow=nrow(matr), ncol=ncol(matr))
  ### for (i in 1:nrow(matr)) for (j in 1:ncol(matr)) quaWeightMatr[i,j] = 1-(abs(i-j))^2/((nrow(matr)-1)^2)
  ### kappaGeneral(m, quaWeightMatr)
  ### CIs: kappa +- z(alpha/2)*sqrt(var)
  if(nrow(matr)!=ncol(matr)) stop("matrix has to be square")
  if (is.null(weightMatr)){
    weightMatr = matrix(0, nrow=dim(matr)[1], ncol=dim(matr)[1])
    for (i in 1:dim(matr)[1]){weightMatr[i,i]=1}
  }
  pMatr = matr/sum(matr)
  Pobs = sum(weightMatr*pMatr)
  sumofrowMatr = matrix(apply(pMatr, 1, sum), nrow=nrow(pMatr), ncol=ncol(pMatr))
  sumofcolMatr = matrix(apply(pMatr, 2, sum), nrow=nrow(pMatr), ncol=ncol(pMatr), byrow=TRUE)
  Pexp = sum(weightMatr*sumofrowMatr*sumofcolMatr)
  ### kappa
  if (sum(diag(matr))==sum(matr)){
    kappa=1
  }else{
    kappa = (Pobs - Pexp)/(1-Pexp)
  }
  ### kappa variance
  weightedRow = apply(sumofcolMatr * weightMatr, 1, sum)
  weightedCol = apply(sumofrowMatr * weightMatr, 2, sum)
  var = 0
  for (i in 1:nrow(matr)){
    for(j in 1:ncol(matr)){
      var = var + pMatr[i, j]*(weightMatr[i, j] - (weightedRow[i] + weightedCol[j]) * (1 - kappa))^2
    }
  }
  var = (var -  (kappa - Pexp*(1-kappa))^2)/((1-Pexp)^2 * sum(matr))
  c(kappa, var)
}

ChronbachsAlpha = function(data, conf = .95){
  ### written for Sunil's grant about three questions in health literacy (project/PILLCVD/othergrant)
  ### different conf. levels are not implemented yet
  ### references:
  ### A. Duhackek, A.T. Coughlan, D. Iacobucci, Results on the Standard Error of the Coefficient of Alpha Index of Reliability, Marketing Science, Vol. 24, No.2 Sprint 2005, pp.294-301
  ### J.M.Van, Zyl, H. Neudecher, D.G. Nel, On the distribution of the Maximum Likelihood Estimator of Cronbach's Alpha, Psychometrika, Vol. 65, No. 3, 271-280, September 2000.
  ### L.J.Cronbach, Coefficient Alpha and the Internal Structure of Tests, Psychometrika, Bol.16, No.3, September 1951
    ## data is a data.frame or matrix with rows, that stand for subjects
    ## and columns that stand for different items (numerical) of the test.
  nonmissingRows = apply(data, 1, function(x){all(!is.na(x))})
  data = data[nonmissingRows, ]
  V = cov(data)
  n = nrow(data)
  p = length(data)
 
  j = matrix(rep(1, p), nrow=p)
  jVj = t(j) %*% V %*% j
  Vsq = V %*% V
  
  ## trace
  #trV = function(matr){s = 0; for (i in 1:min(dim(matr))){s = s + matr[i, i]}; s}
  trV = function(matr){sum(diag(matr))}
  #alpha
  alpha = (p/(p-1))*(1-trV(V)/jVj)
  #Q
  Q = ( (2*(p^2)/(p-1)^2)/(jVj^3) )   *   ( jVj*(trV(Vsq)+(trV(V))^2)     -     2*trV(V)*(t(j)%*%Vsq%*%j))
  #ci
  ci = qnorm(conf+(1-conf)/2)*c(-1,1)*sqrt(Q/n)
  res = c(p, n, alpha, alpha + ci, conf)
  names(res) = c("ques.Num", "subj.Num", "alpha", "ci.lower", "ci.upper", "conf.level")
  res
}

percentAgreement = function(matr){
  ## matr is a matrix M with matr[1,1]=A=(test1pos|test2pos), matr[2,2]=D=(test1neg|test2neg), matr[2,1]=C=(test1pos|test2neg), matr[1,2]=B=(test1neg|test2pos)
  ## returns percent agreement of tests
  if(nrow(matr)!=ncol(matr)) stop("matrix has to be square")
  mtrace = 0
  for(i in 1:nrow(matr)){mtrace = mtrace + matr[i,i]}
  100*mtrace/(sum(matr))
}

###-------------------------------------------------------------------------------------------
###-------------------------------statistics related code-------------------------------------
###-------------------------------------------------------------------------------------------

convertIntoTimeVariant = function(dataList, idName, varNames, varTimeNames, outcomeDataName, outcomeName, outcomeTimeName,
                                  lastObsCarriedForward = FALSE){
  ### converts data into format that can be used in time variant cox regression
  ### dataList  -  is a list of timevariant values and times they were measured, and also outcome and times of outcome
  ### example:
  ### dataList = list(var1 = data.frame(id=c(1,1, 2,2), time = c(.5, 3, 39, 40), v = c(111, 112, 331, 332)),
  ###                 var2 = data.frame(id=c(1,1,1,1, 2,2), time = c(.5, 5, 6, 7, 39, 40), w=c(21, 22, 23, 24, 444, 445)),
  ###                 outcome = data.frame(id=c(1, 2), time=c(5.5, 40), o=c(1, 0)))
  ### idName  -  name of ID variable, should be the same in all data frames in dataList
  ### varNames - names of timevariant variables in dataList (see the example below)
  ### varTimeNames - names of time variables for each timevariant variable (see the example below)
  ### outcomeDataName - name of the outcome data frame in dataList (see the example below)
  ### outcomeName - name of the outcome variable in the outcome data frame (see the example below)
  ### lastObsCarriedForward - please read carefully and also see the example. If a subject is missing the first measurement
  ###                         it remains missing. If the missing measurement is not the first, it can be filled with the 
  ###                         previous measurement if lastObsCarriedForward = TRUE
  ### The function gives an error if the following doesn't hold
  ### outcomeTimeName - name of the time of the outcome variable in the outcome data frame (see the example below)
  ### all data frames in dataList should have the same ID sets (the same IDs)

  ### NOTE:
  ### If a variable was measured at the same time as the outcome it is assumed that this measure happened after the outcome. The reason for this is that variable measures are recorded at a start with start, while outcomes are recorded at an end time. If a variable and an outcome are recored at the same time, start and end times are equal, which is not allowed in Surv() function.
  ### Observations are discarded if recored after the outcome.
  ### The function does not merge baseline variables, it is very easy to do this afterwards, using function merge.

  if(FALSE){
    idName = "id"
    outcomeDataName = "outcome"
    varNames = c("v", "w")
    varTimeNames = c("time", "time")
    outcomeName = "o"
    outcomeTimeName = "time"
    lastObsCarriedForward = FALSE
    dataList = list(var1 = data.frame(id=c(1,1, 2,2, 3, 4,4, 5,5,5), time = c(.5,3, 39,40, .5, .5,1, .5,10,15), v = c(111,112, 331,332, 10, NA,999, 551,552,553)),
                    var2 = data.frame(id=c(1,1,1,1, 2,2, 3,3,3, 4, 5,5,5), time = c(.5, 5, 6, 7, 39, 40, 7,8,9, 1, .5,10,15), w=c(21, 22, 23, 24, 444, 445, 10.1,NA,10.3, 990, 501,502,503)),
                    outcome = data.frame(id=c(1, 2, 3, 4, 5), time=c(5.5, 40, 10, 11, 15), o=c(1, 0, 1, 0, 1)))
    convertIntoTimeVariant(dataList, idName, varNames, varTimeNames, outcomeDataName, outcomeName, outcomeTimeName,
                           lastObsCarriedForward = TRUE)
  }
  varTimeNamesWithOutcome = c(varTimeNames, outcomeTimeName)
  variableDataNames = setdiff(names(dataList), outcomeDataName)

  ### all data frames in dataList should have the same ID sets (the same IDs)
  ### ids in all data frames should NOT be class factor
  uniqueIDs = unique(unlist(sapply(dataList, function(x){x[[idName]]})))
  for (n in names(dataList)){
    if (!setequal(dataList[[n]][[idName]], uniqueIDs)) stop(paste("Dataframe", n, "has different set of IDs. Please make sure that all dataframes in dataList have the same sets of IDs."))
    if (class(dataList[[n]][[idName]]) == "factor") stop(paste("It dataframe", n, "IDs are of class 'factor'. Please convert them either to numeric or character."))
  }
  ### the data frame with the outcome has to have one record per ID.
    if (length(unique(dataList[[outcomeDataName]][[idName]])) != length(dataList[[outcomeDataName]][[idName]])){stop("The data frame with the outcome has to have one record per ID.")}

  ### define a value that none of the variables get:
  if (!lastObsCarriedForward){
    naValue = max(unlist(sapply(1:length(varNames), function(i){dataList[[variableDataNames[i]]][[varNames[i]]]})), na.rm=TRUE) + 999
  }

  ### define the length of the resulting data frame:
  dataLengths = sapply(dataList, nrow)
  resNrow = sum(dataLengths)
  resNrow
  eachDataStart = c(1, sapply(1:length(dataList), function(i){1+sum(dataLengths[1:i])})[1:(length(dataList)-1)])
  eachDataStart
  ### create a dataframe containing ids and start times for ALL data in dataList
  timeAndID = data.frame(id=rep(NA, resNrow), startTime=rep(NA, resNrow))
  names(timeAndID) = gsub("id", idName, names(timeAndID))

  ### fill out that dataframe with ids and start times
  for (i in 1:length(dataList)){
    startIndex = eachDataStart[i]
    endIndex = eachDataStart[i] + nrow(dataList[[i]]) - 1
    timeAndID[[idName]][startIndex:endIndex] = dataList[[i]][[idName]]
    timeAndID$startTime[startIndex:endIndex] = dataList[[i]][[varTimeNamesWithOutcome[i]]]
  }
  timeAndID = timeAndID[order(timeAndID[[idName]], timeAndID$startTime),]
  timeAndID = unique(timeAndID)

  ### add end time
  timeAndID$endTime = c(timeAndID$startTime[2:nrow(timeAndID)], NA)
  #timeAndID
  currID = timeAndID[[idName]][nrow(timeAndID)]
  for (i in (nrow(timeAndID)-1):1){
    if (currID != timeAndID[[idName]][i]){
      timeAndID$endTime[i] = NA
      currID = timeAndID[[idName]][i]
    }
  }
  #timeAndID

  #### merge all timevariant variables
  resData = timeAndID
  resData$idStartTime = paste(resData[[idName]], resData$startTime)
  resData$idEndTime = paste(resData[[idName]], resData$endTime)
  for (i in 1:(length(dataList)-1)){
    mergingData = dataList[[i]][, c(idName, varNames[i], varTimeNames[i])]
    if (!lastObsCarriedForward){
      mergingData[[varNames[i]]][is.na(mergingData[[varNames[i]]])] = naValue
    }
    mergingData$idStartTime = paste(mergingData[[idName]], mergingData[[varTimeNames[i]]])
    resData = merge(resData, mergingData[, setdiff(names(mergingData), c(idName, varTimeNames[i]))], by = "idStartTime", all.x=TRUE)
  }

  #### merge outcome
  mergingData = dataList[[outcomeDataName]][, c(idName, outcomeName, outcomeTimeName)]
  mergingData$idEndTime = paste(mergingData[[idName]], mergingData[[outcomeTimeName]])
  resData = merge(resData, mergingData[, setdiff(names(mergingData), c(idName, outcomeTimeName))], by = "idEndTime", all.x=TRUE)

  ### set undefined vars with the previous value
  resData = resData[order(resData[[idName]], resData$startTime),]
  currID = resData[[idName]][1]
  currVars = resData[1, varNames]
  for (i in 2:nrow(resData)){
    if (currID != resData[[idName]][i]){
      currID = resData[[idName]][i]
      currVars = resData[i, varNames]
    }else{
      resData[i, varNames][is.na(resData[i, varNames])] = currVars[is.na(resData[i, varNames])]
      currVars = resData[i, varNames]
    }
  }
  #resData[, c(idName, "startTime", "endTime", varNames, outcomeName)]

  ### set undefined outcome with 0's
  currID = resData[[idName]][nrow(resData)]
  passedTheOutcome = FALSE
  for (i in (nrow(resData)-1):1){
    if (currID != resData[[idName]][i]){
      currID = resData[[idName]][i]
      passedTheOutcome = !is.na(resData[[outcomeName]][i])
    }else{
      if (passedTheOutcome){
        resData[[outcomeName]][i][is.na(resData[[outcomeName]][i])] = 0
      }else{
        passedTheOutcome = !is.na(resData[[outcomeName]][i])
      }
    }
  }
  #resData

  ### if lastObsCarriedForward = FALSE, return all the missing values back to missing
  if (!lastObsCarriedForward){
    for (n in varNames){
      resData[[n]][resData[[n]] == naValue] = NA
    }
  }
  resData = resData[!is.na(resData[[outcomeName]]), c(idName, "startTime", "endTime", varNames, outcomeName)]
  resData
}

# bootstrapByID = function(data, idVar){
#   ## the function bootstraps rows per each id and returns
#   ## a data.frame with the bootstrapped rows (the same # of rows per id as in the original data) and attached ids
#   lengthPerID = tapply(data[[idVar]], data[[idVar]], length)
#   #data$lengthPerID = lengthPerID[as.character(data[[idVar]])]
#   someorder = sapply(lengthPerID, function(x){1:x})
#   data$someorder = NA
#   data = data[order(data[[idVar]]),]
#   for (i in names(someorder)){
#     data[as.character(data[[idVar]])==i,][["someorder"]] = someorder[[i]]
#   }
#   #data$idAndSomeorder = paste(data[[idVar]], data$someorder)
#   rownames(data) = paste(data[[idVar]], data$someorder)
#   
#   bootstrappedData = data.frame(id=data[[idVar]])
#   names(bootstrappedData) = idVar
#   for (n in setdiff(names(data), idVar)){
#     bootstrappedData[[n]] = rep(NA, nrow(data))
#   }
#   
#   freeIndex = 1
#   for (i in unique(data[[idVar]])){
#     len = lengthPerID[i]
#     #bootIDAndSomeOrder = sample(data$idAndSomeorder[data[[idVar]]==i], len, replace = TRUE)
#     bootIDAndSomeOrder = sample(rownames(data[data[[idVar]]==i,]), len, replace = TRUE)
#     bootRows = data[c(bootIDAndSomeOrder),]
#     bootstrappedData[freeIndex:(freeIndex+len-1),] = bootRows
#     freeIndex = freeIndex + len
#   }
#   bootstrappedData$lengthPerID=NULL
#   bootstrappedData$someorder=NULL
#   #bootstrappedData$idAndSomeorder=NULL
#   rownames(bootstrappedData) = 1:nrow(bootstrappedData)
#   bootstrappedData
# }
# bootstrapByID = function(data, idVar){
#   ## the function bootstraps rows per each id and returns
#   ## a data.frame with the bootstrapped rows (the same # of rows per id as in the original data) and attached ids
#   lengthPerID = tapply(data[[idVar]], data[[idVar]], length)
#   uniqueIDs = setdiff(unique(data[[idVar]]), NA)
#   bootstrappedIDs = sample(x=uniqueIDs, size=length(uniqueIDs), replace=TRUE)
#   bootstrappedData = data[data[[idVar]]==bootstrappedIDs[1],]
#   for (i in 2: length(bootstrappedIDs)){
#     bootstrappedData = rbind(bootstrappedData, data[data[[idVar]]==bootstrappedIDs[i],])
#   }
#   bootstrappedData
# }

bootstrapSubjectsByIndex = function(data, idVar, bootTimes=200){
  ## the function bootstraps rows per each id and returns
  ## a data.frame with the bootstrapped rows (the same # of rows per id as in the original data) and attached ids
  data = data[!is.na(data[[idVar]]),]
  indPerID = tapply(1:nrow(data), data[[idVar]], function(x){x})
  lengthPerID = sapply(indPerID, length)
  for (bi in 1:bootTimes){
    bootstrappedIDIndices = sample(x=1:length(indPerID), size=length(indPerID),
                                   replace=TRUE)
    newDataLength = sum(lengthPerID[bootstrappedIDIndices], na.rm=TRUE)
    allBootstrappedIndices = rep(NA, newDataLength)
    endInd = 0
    for (i in 1:length(bootstrappedIDIndices)){
      startInd = endInd + 1
      endInd = startInd + lengthPerID[bootstrappedIDIndices[i]] - 1
      allBootstrappedIndices[startInd:endInd] = indPerID[[bootstrappedIDIndices[i]]]
    }
    res = data[allBootstrappedIndices, ]
  }
}

normalCIsForMean = function(vec, confLevel=.95){
  ### After removing NAs, returns CIs based on Normal distribution
  ### example:
  ### vec = c(2.2, 2.4, 4.9, 2.5, 3.7, 4.3)
  ### normalCIsForMean(vec) returns c(2.429953, 4.236714)
  vecNoMiss = vec[!is.na(vec)]
  n = length(vecNoMiss)
  sampleMean = mean(vecNoMiss)
  sampleVar = var(vecNoMiss)
  if (n<60) cat("The number of non-missing values is less than 60.\nConsider using gossettCIsForMean() or bootstrapCIsForMean()\n")
  res = sampleMean + c(qnorm((1-confLevel)/2), qnorm((1+confLevel)/2)) * sqrt(sampleVar/n)
  res
}

gossettCIsForMean = function(vec, confLevel=.95){
  ### After removing NAs, returns CIs based on student distribution
  ### example (see Medical Statistics by Kirkwood & Sterne):
  ### vec = c(2.2, 2.4, 4.9, 2.5, 3.7, 4.3)
  ### gossettCIsForMean(vec) returns c(2.148509 4.518158)
  vecNoMiss = vec[!is.na(vec)]
  n = length(vecNoMiss)
  sampleMean = mean(vecNoMiss)
  sampleVar = var(vecNoMiss)
  if (n<30) cat("The number of non-missing values is less than 30.\nConsider using bootstrapCIsForMean()\n")
  res = sampleMean + c(qt((1-confLevel)/2, n-1), qt((1+confLevel)/2, n-1)) * sqrt(sampleVar/n)
  res
}

bootstrapCIsForMean = function(vec, confLevel=.95, replicationsNum = 999){
  ### After removing NAs, returns CIs based on non-parametric bootstrap
  ### See "Bootsrap methods and their applications" by Davison, Hinkley, p.28
  ### example:
  ### vec = c(2.2, 2.4, 4.9, 2.5, 3.7, 4.3)
  ### set.seed(1)
  ### bootstrapCIsForMean(vec) returns c(2.516667 4.100000)
  halfAlpha = (1 - confLevel)/2
  vecNoMiss = vec[!is.na(vec)]
  sampleMean = mean(vecNoMiss)
  n = length(vecNoMiss)
  repMean = rep(NA, replicationsNum)
  for (i in 1:replicationsNum){
    repMean[i] = mean(sample(vecNoMiss, size=n, replace=TRUE))
  }
  repMean = repMean[order(repMean)]
  res = rep(NA, 2)
  res[1] = 2*sampleMean - repMean[ceiling((replicationsNum+1)*(1-halfAlpha))]
  res[2] = 2*sampleMean - repMean[floor((replicationsNum+1)*(halfAlpha))]
  #res = sampleMean - c(repMean[ceiling((replicationsNum+1)*confLevel)], repMean[floor((replicationsNum+1)*(1-confLevel))])
  res
}

wilsonCIsForProportion = function(p, n, confLevel=.95){
  ### from book Statistica methods for rates and proportions by Fleiss, Levin, and Cho Paik
  z = qnorm((1+confLevel)/2)
  q = 1-p
  ciL = ((2*n*p + z^2 - 1) - z*sqrt(z^2 - (2 + 1/n) + 4*p*(n*q+1)))/(2*(n + z^2))
  ciH = ((2*n*p + z^2 + 1) + z*sqrt(z^2 + (2 - 1/n) + 4*p*(n*q-1)))/(2*(n + z^2))
  c(ciL, ciH)
}

CIsForSurvivalFunctionNoCensoring = function(initialNumberOfSubjects, times, confLevel=.95){
  # some lecture notes, the files is called CIForSurvivalCurv.pdf
  #initialNumberOfSubjects = 240
  #times = c(rep(1, 12), rep(2, 9), rep(3, 17), rep(4, 36), rep(5,6), rep(6, 18), rep(7, 13), rep(8, 11), rep(9, 14), 33, 35, 36)
  z = qnorm((1+confLevel)/2)
  data = data.frame(times = times, events = 1, atrisk = 0)
  #surv.all = survfit(Surv(times, event=events)~1, data=data)    
  data = data[order(data$times), ]
  eventsPerTime = tapply(data$times, data$times, length)
  data$events = eventsPerTime[as.character(data$times)]
  data = unique(data)
  data = data[order(data$times), ]
  data$atrisk = initialNumberOfSubjects
  data$surv[1] = (data$atrisk[1] - data$events[1])/data$atrisk[1]
  data$cumulevents = data$events
  data$CIUPPER = NA
  data$CILOWER = NA
  data$GreedwoodCILower = NA
  data$GreedwoodCIUpper = NA
  for (i in 2:nrow(data)){
    data$cumulevents[i] = data$cumulevents[i-1] + data$events[i]
    data$atrisk[i] = initialNumberOfSubjects - data$cumulevents[i-1]
    data$surv[i] = data$surv[i-1]*(data$atrisk[i] - data$events[i])/data$atrisk[i]
  }
  sdLt = rep(NA, nrow(data))
  GreedwoodSD = rep(NA, nrow(data))
  for (i in 1:nrow(data)){
    numerator  = sum(data$events[1:i]/(data$atrisk[1:i]*(data$atrisk[1:i] - data$events[1:i])))
    denominator = log(data$surv[i])^2
    GreedwoodSD[i] = sqrt(numerator) * data$surv[i]
    sdLt[i] = sqrt(numerator/denominator)
  }
  data$GreedwoodCILower = data$surv - z*GreedwoodSD
  data$GreedwoodCIUpper = data$surv + z*GreedwoodSD
  data$Lt = log(-log(data$surv))
  ltLower = data$Lt - z*sdLt
  ltUpper = data$Lt + z*sdLt
  data$CILOWER = exp(-exp(ltUpper))
  data$CIUPPER = exp(-exp(ltLower))
  data$GreedwoodCILower[data$GreedwoodCILower<0] = 0
  data$GreedwoodCIUpper[data$GreedwoodCIUpper>1] = 1
  data
}

regressionImpFromOLS = function(modelFormulaStr, idVar, dataToImputeFrom, howManyImp = 10){
  ### takes a model of varToImpute ~ dataToImpute
  ### and returns several predicted values of the varToImpute
  ### return object is a matrix with nrow of dataToImpute and ncol of howManyImp
  ### and rownames that are obtained from idVar
  ### so it would be possible to merge predicted data with the dataframe
  ### idVar has to contain unique values.
  ## number of observations
  if(length(dataToImputeFrom[[idVar]]) != length(unique(dataToImputeFrom[[idVar]]))) stop("idVar has to contain unique values.")
  varToImpute = gsub(" +", "", unlist(strsplit(as.character(modelFormulaStr), split="~"))[1])
  dataToRunTheModel = dataToImputeFrom[!is.na(dataToImputeFrom[[varToImpute]]),]
  dataToUseToPredict = dataToImputeFrom[is.na(dataToImputeFrom[[varToImpute]]),]
  modelOne = ols(as.formula(modelFormulaStr), data = dataToRunTheModel)
  res = data.frame(imputeFromOLS(modelOne, dataToImpute = dataToUseToPredict, howManyImp = howManyImp))
  res[[idVar]] = dataToUseToPredict[[idVar]]
  res
}

imputeFromOLS = function(linearModel, dataToImpute, howManyImp = 10){
  ### imputes from a model howManyImp times
  ### imputes from a model howManyImp times
  ### imputes from a model howManyImp times
  nj = length(linearModel$residuals) 
  ## number of variables in the model including non-linear terms
  k = length(linearModel$coefficients)-1 
  sigmaSquare = sum(linearModel$residuals^2)/(nj-k-1)
  g = rchisq(n=1, df=nj-k-1)
  sigmaSquareStar = sigmaSquare*(nj-k-1)/g
  ############ prepare matrix for Choleski decomposition, 1/(X'X) page 166-167 of Rubin
  X = (model.matrix(linearModel))
  InvOfXPrimeX = solve(t(X) %*% X)
  choleskiDecUppTrian = chol(InvOfXPrimeX)
  res = matrix(NA, nrow=nrow(dataToImpute), ncol=howManyImp)
  ############ draw betas
  for (i in 1:howManyImp){
    zVars = rnorm(k+1, 0, 1)
    sampledCoef = linearModel$coefficients + sqrt(sigmaSquareStar) * 
                  choleskiDecUppTrian %*% matrix(zVars, ncol=1)
    ############ put coefficients into the model object to get predicted values
    names(sampledCoef) = names(linearModel[["coefficients"]])
    sampledModel = linearModel
    sampledModel[["coefficients"]] = sampledCoef
    res[,i] = predict(sampledModel, newdata = dataToImpute) + sqrt(sigmaSquareStar) * rnorm(n=nrow(dataToImpute))
  }
  res
}


if(FALSE){
  initialNumberOfSubjects = 240
  times = c(rep(1, 12), rep(2, 9), rep(3, 17), rep(4, 36), rep(5,6), rep(6, 18), rep(7, 13), rep(8, 11), rep(9, 14), 39, 40, 44, 45, 47)
  #times = c(.5, rep(1, 2), 3)
  tmpdata = data.frame(times = times, events = 1, atrisk = 0)
  surv.all = survfit(Surv(time=times, event=events)~1, data=tmpdata)
  summary(surv.all) 
  tmp = CIsForSurvivalFunctionNoCensoring(initialNumberOfSubjects=length(times), times, confLevel=.95)
  tmp
  par(mfrow=c(2, 1))
  plot(tmp$times, tmp$surv, type="l", col="gray", ylim=c(0,1))
  lines(tmp$times, tmp$CILOWER, col="red")
  lines(tmp$times, tmp$CIUPPER, col="blue")
  plot(tmp$times, tmp$surv, type="l", col="gray", ylim=c(0,1))
  lines(tmp$times, tmp$GreedwoodCILower, col="red")
  lines(tmp$times, tmp$GreedwoodCIUpper, col="blue")
  d1 = read.csv("tmp.csv")
  system.time(tmp<-CIsForSurvivalFunctionNoCensoring(initialNumberOfSubjects=nrow(d1), d1[,1], confLevel=.95))
  tmp = CIsForSurvivalFunctionNoCensoring(initialNumberOfSubjects=nrow(d1), d1[,1], confLevel=.95)
  tmp1 = tmp[tmp$times<30 & tmp&times>29,]
  range(tmp1$times, na.rm=TRUE)
  range(tmp1$surv, na.rm=TRUE)
}
# set.seed(1)
# 
# replVec = c(19, 39, 99, 199, 499, 999, 1999, 5999)
# somenum = 10
# samplesize = 30000
# cil = matrix(NA, ncol=length(replVec), nrow=somenum)
# cih = matrix(NA, ncol=length(replVec), nrow=somenum)
# for (i in 1:length(replVec)){
#   for(j in 1:somenum){
#     vec = rnorm(samplesize)
#     #res = bootstrapCIsForMean(vec, replicationsNum=replVec[i])
#     res = gossettCIsForMean(vec)
#     cil[j, i] = res[1]
#     cih[j, i] = res[2]
#   }
# }

# plot(cil[1, ], ylim=c(min(cil, cih), max(cil, cih)), type="l")
# #plot(cil[1, ], ylim=c(-1,1), type="l")
# for (i in 2:somenum){
#   lines(cil[i, ])
# }
# for (i in 2:somenum){
#   lines(cih[i, ])
# }
# abline(h=c(qnorm(.975)/sqrt(samplesize), qnorm(.025)/sqrt(samplesize)), col="red", lwd=2)
# lines(apply(cil, 2, mean), col="blue", lwd=2)
# lines(apply(cih, 2, mean), col="blue", lwd=2)
# 
# normalCIsForMean(vec)
# gossettCIsForMean(vec)
# bootstrapCIsForMean(vec)

if(FALSE){ ### old code not in use
#   convertToCoxTimeVarientFormat = function(data, idVarName, startVarName, newStopVarName, newStartVarName, timeVarName, statusVarsName, eventVarName){
#     ### data = data.frame(id=c(1,1,1,1, 3, 4,4,5 ,6,6,6, 7,7,7, 8,8, 9,9,9), start = 0, time=c(1,2,3,4, 3, 1,4, 6, 3,7,10, 0,1,2, 0,4, 1,2,3), st1=c(1,1,3,3, 3, 1,1, 3, 3,1,1, 1,1,1, 3,3, 1, 3, 1), event=c(0,0,1,0, 0, 0,1, 1, 0,0,1, 1,0,0, 0,0, 0, 1, 0))
#     ### idVarName="id"; startVarName="start"; timeVarName="time"; statusVarsName = "st1"; eventVarName = "event"; newStartVarName="newstart"; newStopVarName="stop"
#     if (any(data[[startVarName]] > data[[timeVarName]])) stop (paste(startVarName, "cannot be larger than", timeVarName))
#     data = data[order(data[[idVarName]], data[[timeVarName]]), ]
#     data[[newStartVarName]] = data[[startVarName]]
#     data[[newStopVarName]] = c(NA, data[[timeVarName]][1:(nrow(data)-1)])
#     currentID = data[[idVarName]][1]
#     for (i in 1:nrow(data)){
#       if (currentID != data[[idVarName]][i]){
#         data[[newStopVarName]][i] = NA
#         currentID = data[[idVarName]][i]
#       }
#     }
#     data[[newStartVarName]][!is.na(data[[newStopVarName]])] = data[[newStopVarName]][!is.na(data[[newStopVarName]])]
#     data[[newStopVarName]] = data[[timeVarName]]
#     data
#   }
# 
#   completeStatusInCoxTimeVarientFormat = function(data, idVarName, startVarName, stopVarName, statusVarsNames){
#     ### data = data.frame(id=c(1,1,1,1, 3, 4,4,5 ,6,6,6, 7,7,7, 8,8, 9,9,9), start = 0, time=c(1,2,3,4, 3, 1,4, 6, 3,7,10, 0,1,2, 0,4, 1,2,3), st1=c(1,1,3,3, 3, 1,1, 3, 3,1,1, 1,1,1, 3,3, 1, 3, 1), event=c(0,0,1,0, 0, 0,1, 1, 0,0,1, 1,0,0, 0,0, 0, 1, 0))
#     ### idVarName="id"; startVarName="start"; timeVarName="time"; statusVarsName = "st1"; eventVarName = "event"; newStartVarName="newstart"; newStopVarName="stop"
#     ### data = convertToCoxTimeVarientFormat(data, idVarName="id", startVarName="start", timeVarName="time", statusVarsName = "st1", eventVarName = "event", newStartVarName="newstart", newStopVarName="stop")
#     ### data$st1[c(3, 5, 10, 11, 15)] = NA
#     data = data[order(data[[idVarName]], data[[stopVarName]]), ]
#     for (st in statusVarNames){
#       currentID = data[[idVarName]][1]
#       currentSt = data[[st]][1]
#       for (i in 2:nrow(data)){
#         if (currentID != data[[idVarName]][i]){
#           currentID = data[[idVarName]][i]
#           currentSt = data[[st]][i]
#         }else{
#           if (is.na(data[[st]][i])){
#             data[[st]][i] = currentSt
#           }else{
#             if (is.na(currentSt) | data[[st]][i] != currentSt){
#               currentSt = data[[st]][i]
#             }
#           }
#         }
#       }
#     }
#     data
#   }
}

