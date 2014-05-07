####################### author: Svetlana K. Eden, svetlana.eden@vanderbilt.edu
################################################

### defineDistance from the center of the group
distFromCenter = function(x, y, xSet=x, ySet=y, precision=12){
  if(FALSE){
    x = rnorm(500)
    y = rnorm(500)
    dist = distFromCenter(x, y, precision=36)
    plot(x, y, col=paste(gray(0), transparency(dist), sep=""), pch=19)
    plot(x, y, col=paste(gray(0), transparency(-dist^(1/4)), sep=""), pch=19)
  }
  ### precision - how many different lines you want to draw to separate the dot from the rest of the set.
  if (length(x)!=length(y) | length(xSet)!=length(ySet)){ stop("Lengths of x and y (xSet and ySet) should be the same.")}
  ### slopes
  Ks = seq(0, (2*pi - 2*pi/precision), 2*pi/precision)
  ### intercepts
  dist = rep(NA, length(x))
  xySetMatr = cbind(xSet, ySet, one=rep(1, length(xSet)))
  for (i in 1:length(x)){
    Cs = y[i] - Ks*x[i]
    KsCsMatr = rbind(Ks, minusOne = rep(-1, length(Ks)), Cs)
    location = xySetMatr %*% KsCsMatr
    distance = apply(location>0, 2, sum)
    dist[i] = min(distance, length(xSet) - distance)
  }
  dist
}

transparency = function(x, minimum = 255, maximum = 1){
  ### takes numbers and converts them to transparency col (00 - FF)
  as.hexmode(round( (x)*(minimum - maximum)/(diff(range(x))) - (minimum - maximum)/(diff(range(x))) *min(x) + maximum))
}

missingVariables = function(data, roundn=3){
  tmp = sapply(data, function(x){sum(is.na(x))})
  if (all(tmp==0)){
    res = "There are no variables with missing values"
    return(res)
  }
  tmp = tmp[tmp>0]
  res = matrix(tmp, nrow=length(tmp))
  rownames(res) = names(tmp)
  res = cbind(res, round(100*tmp/nrow(data), roundn))
  colnames(res) = c(paste("# of NA out of", nrow(data)), "  % of NA")
  return(res)
}


myformula = function(outcome, covariates, splineCovarsWithNumOfKnots=NULL, interactions = NULL){
  ### outcome is a character string
  ### when used in models (lrm, ols, ... ) needs to be converted to formula using function as.formula()
  ### covariates is a list of characters
  ### splineCovarsWithNumOfKnots is a list with names as covariates and number of knots for each name
  ### for example splineCovarsWithNumOfKnots = list(age=3, literacy=24)
  ### interactions is a list with character vectors
  ### for example: list(c("age", "tx"), c("visit", "a1c", "literacy"))
  if(FALSE){ ### example:
    outcome = "someScore"
    covariates = c("age", "tx", "visit", "a1c", "literacy")
    myformula(outcome, covariates)

    splineCovarsWithNumOfKnots = list()
    splineCovarsWithNumOfKnots[["age"]] = 3
    splineCovarsWithNumOfKnots[["literacy"]] = 4
    myformula(outcome, covariates, splineCovarsWithNumOfKnots = splineCovarsWithNumOfKnots)

    interactions = list(c("age", "tx"), c("visit", "a1c", "literacy"))
    myformula(outcome, covariates, interactions = interactions)

    myformula(outcome, covariates, splineCovarsWithNumOfKnots = splineCovarsWithNumOfKnots, interactions = interactions)
  }
  nonLinearCovarsStr = c()
  interactStr = c()
  covariatesStr = covariates
  if (!is.null(splineCovarsWithNumOfKnots)){
    rcsCovars = unlist(splineCovarsWithNumOfKnots)
    if (any(rcsCovars<3)) stop("splines cannot have less than three knots")
    if (!all(names(rcsCovars) %in% covariates)) stop(paste("The following spline covariates are not in the list of the main covariates:", paste(names(rcsCovars)[names(rcsCovars) %nin% covariates])))
    covariatesStr = setdiff(covariates, c(names(rcsCovars), outcome))
    nonLinearCovarsStr = paste("rcs(", names(rcsCovars), ",", rcsCovars, ")", sep="")
    names(nonLinearCovarsStr) = names(rcsCovars)
  }
  if (!is.null(interactions)){
    interCheck = unlist(interactions)
    if (any(interCheck %nin% covariates)) stop(paste("The following interaction covariates are not in the list of the main covariates:", paste(interCheck[interCheck %nin% covariates])))
    if (!is.null(splineCovarsWithNumOfKnots)){
      interWithNonLinear = sapply(interactions, function(x){if (any(x %in% names(nonLinearCovarsStr))) c(x[x %nin% names(nonLinearCovarsStr)], nonLinearCovarsStr[x[x %in% names(nonLinearCovarsStr)]]) else x})
    }else{
      interWithNonLinear = interactions
    }
    if (length(interactions) > 1){
      interactStr = sapply(interWithNonLinear, paste, collapse="*")
    }else{
      interactStr = paste(unlist(interWithNonLinear), collapse="*")
    }
  }
  formRes = paste(outcome, "~", paste(c(covariatesStr, nonLinearCovarsStr, interactStr), collapse=" + "))
  formRes 
}

### non-statistical and non-latex stuff

checkDerivativeVariables = function(data, componentsVars, derivativeVars, message, byVar=NULL, numToPick=2){
### componentsVars - what variables the score was created from
### derivativeVars - what were created from componentsVars
### byVar - if not NULL, the function will also give examples by different levels of this variable (must be the same length as "data")
### numToPick - how many rows to sample for each frequency of missing componentsVars
  # depressScaleVars = c("ces_d_1", "ces_d_2", "ces_d_3", "ces_d_4", "ces_d_5", "ces_d_6", "ces_d_7", "ces_d_8", "ces_d_9", "ces_d_10", "ces_d_11", "ces_d_12", "ces_d_13", "ces_d_14", "ces_d_15", "ces_d_16", "ces_d_17", "ces_d_18", "ces_d_19", "ces_d_20")
  # derivativeVars = c("cesd")
  if (length(componentsVars)<2){
    howManyMissing = as.numeric(is.na(data[, componentsVars]))
  }else{
    howManyMissing = apply(data[, componentsVars], 1, function(x)sum(is.na(x)))
  }
  howManyMisTable = table(howManyMissing)
  howManyMissMatr = matrix(c(as.numeric(names(howManyMisTable)), howManyMisTable), nrow=2, byrow = TRUE)
  rownames(howManyMissMatr) = c("# of missing items", "# of subjects with missing items")
  colnames(howManyMissMatr) = rep("", ncol(howManyMissMatr))
  set.seed(45)
  randomPickForEachMissNum = unlist(tapply(1:nrow(data), howManyMissing, function(x){if (length(x)==1) return(x) else sample(x=x, size=min(numToPick, length(x)), replace = FALSE)}))
  outputdata = cbind(data[randomPickForEachMissNum, c("id", componentsVars, derivativeVars)])
  rownames(outputdata) = 1:nrow(outputdata)
  cat("\n\n", message, "\n", paste(rep("-", nchar(message)), collapse=""), "\n", sep="")
  print(howManyMissMatr)
  cat("\nExamples of variable/s:", paste(derivativeVars, collapse=", "), "\nfor each case of missing:\n")
  print(outputdata)
  if (!is.null(byVar)){
    numToPick1 = 1
    if (nrow(data) != length(data[[byVar]])) stop("\nLength of variable byVar is not equal to number of rows in data")
    set.seed(4309)
    randomPickForEachVarVal = unlist(tapply(1:nrow(data), data[[byVar]], function(x){if (length(x)==1) return(x) else sample(x=x, size=min(numToPick1, length(x)), replace = FALSE)}))
    outputdata = cbind(data[randomPickForEachVarVal, c("id", componentsVars, derivativeVars)])
    rownames(outputdata) = 1:nrow(outputdata)
    cat("\nExamples of variable/s:", paste(derivativeVars, collapse=", "), "\nfor each value of variable", byVar, "\n")
    print(outputdata)
  }
}


lookUpInMonotonVec = function(vector, value){
  ### works on a numeric ordered vector with unique values
  ### returns a matrix of c(index, index) of value in vector if vector contains value,
  ### If there exist two values (v1, v2) in vector such that
  ### v1< value <v2, then returns indices of these two vector elements
  ### if value is greater than the max(vector) returns the index of the last (largest) element of vector and NA
  ### if value is less than the min(vector) returns NA and the index of the first (smallest) element of vector
  vector = as.numeric(vector)
  value = as.numeric(value)
  if (any(is.na(vector))) stop ("\"vector\" should not have missing values")
  if (any(vector!=vector[order(vector)])) stop ("\"vector\" must be ordered")
  if (length(vector)!=length(unique(vector))) stop ("All values of \"vector\" must be unique")
  retMatr = matrix(NA, nrow = length(value), ncol=2)
  for (i in 1:length(value)){
    if (is.na(value[i])){
      retMatr[i, ] = c(NA, NA)
    }else{
      matchIndex = match(value[i], vector)
      if (!is.na(matchIndex)){
        retMatr[i, ] = c(matchIndex, matchIndex)
      }else{
        if (value[i] < vector[1]){
          retMatr[i, ] = c(NA, 1)
        }else{
          if (value[i] > vector[length(vector)]){
            retMatr[i, ] = c(length(vector), NA)
          }else{
            secondIndex = match(TRUE, as.numeric(vector) > value[i])
            retMatr[i, ] = c(secondIndex-1, secondIndex)
          }
        }
      }
    }
  }
  retMatr
}
 
dumpToCSV = function(dataframe, fileName, sep=",", quotes="\""){
  bigsep = paste(quotes, sep, quotes, sep="")
  for (n in names(dataframe)){
    dataframe[[n]] = as.character(dataframe[[n]])
  }
  cat(quotes, paste(names(dataframe), collapse=bigsep), quotes, "\n", sep="", file=fileName)
  for (i in 1:nrow(dataframe)){
    str = paste(quotes, paste(dataframe[i,], collapse=bigsep), quotes, "\n", sep="")
    str = gsub("\"NA\"", "\"\"", str)
    cat(str, file=fileName, append=TRUE)
  }
}

dumpToCSVForSAS = function(dataframe, fileName){
  ### the difference b/w dumpToCSV and this function is that this function wrap in quotes only factor or character variales
  sep=","
  quotes=rep("\"", ncol(dataframe))
  varClass = sapply(dataframe, class)
  quotes[!(varClass == "character" | varClass == "factor")] = ""
  bigsep = paste(quotes, sep, quotes, sep="")
  for (n in names(dataframe)){
    dataframe[[n]] = as.character(dataframe[[n]])
  }
  #cat(quotes, paste(names(dataframe), collapse=bigsep), quotes, "\n", sep="", file=fileName)
  #cat(paste(names(dataframe), collapse=bigsep), "\n", sep="", file=fileName)
  cat(paste(names(dataframe), collapse=sep), "\n", sep="", file=fileName)
  for (i in 1:nrow(dataframe)){
    #str = paste(quotes, paste(dataframe[i,], collapse=bigsep), quotes, "\n", sep="")
    #str = paste(paste(quotes, dataframe[i,], quotes, sep=""), collapse=bigsep)
    str = paste(paste(quotes, dataframe[i,], quotes, sep=""), collapse=sep)
    str = gsub("\"NA\"", "\"\"", str)
    str = gsub("(?!=,)NA(?!>,)", "", str, perl=TRUE)
    #str = gsub(",NA,", ",,", str)
    cat(str, "\n", file=fileName, append=TRUE)
  }
}

substituteSKE = function(vector, what, towhat, therest=NULL){
### function substitutes the 'what' values in 'vector' to towhat 'values' in this vector.
### if other values in 'vector' are not in 'what' then put therest instead
### if 'therest' is set to NULL, makes sure all values in 'vector' are mentioned in 'what'
  if (is.null(therest)){
    #stop()
    #cat(paste(unique(vector), collapse=","), length(setdiff(vector, c(what, NA))), "=============\n")
    if (length(setdiff(vector, c(what, NA)))>0) stop("variable \"what\" has to contain all unique values of \"vector\"")
  }
  if (length(what)!=length(towhat)){stop("the lengths of \"what\" and \"towhat\" should be the same")}
  substvector = towhat
  names(substvector) = as.character(what)
  res = substvector[as.character(vector)]
  if (!is.null(therest)){
    res[is.na(res)] = therest
  }
  names(res) = NULL
  res
}

kama=function(vec){
  length(unique(vec))
}

fw = function(vector, str){
  grep(str, vector, ignore.case=TRUE, value=TRUE)
}

easyRectUniv = function(startXY, xLen, yLen, xDirection = "right", yDirection="up", color="gray",...){
### plots rectangle with given coordinates startXY=c(. , .)
### where startXY is its beginning point
### xDirection = "right" where to draw the rect. to the right or to the left or both (center alignment) from startXY
### yDirection = "up" where to draw the rect. up or down or both (center alignment) from startXY
  xDirectionVector = c(0, 0, 0, 0)
  if (xDirection == "right"){
    xDirectionVector = c(0, 1, 1, 0)
  }else{
    if (xDirection == "left"){
      xDirectionVector = c(0, -1, -1, 0)
    }else{
      if (xDirection == "both"){
        xDirectionVector = .5*c(-1, 1, 1, -1)
      }else{
        stop("xDirection can be 'right', 'left', or 'both'")
      }
    }
  }
  yDirectionVector = c(0, 0, 0, 0)
  if (yDirection == "up"){
    yDirectionVector = c(0, 0, 1, 1)
  }else{
    if (yDirection == "down"){
      yDirectionVector = c(0, 0, -1, -1)
    }else{
      if (yDirection == "both"){
        yDirectionVector = .5*c(1, 1, -1, -1)
      }else{
        stop("yDirection can be 'up', 'down', or 'both'")
      }
    }
  }
  polygon(x = startXY[1] + xDirectionVector*xLen, y = startXY[2] + yDirectionVector*yLen, col=color, ...)
}

easyRectUnivMany = function(startX, startY, xLen, yLen, xDirection = "right", yDirection="up", color=rep("gray", length(startX)), ...){
### plots rectangle with given coordinates startX, startY
### where startX, startY are vectors where the rectangles start
### xLen = vector of lengths in x dimensions
###    if constant then all rectangles are of the same x length
### yLen = vector of lengths in y dimensions
###    if constant then all rectangles are of the same y length
### xDirection = "right" where to draw the rect. to the right or to the left or both (center alignment)
###    if constant then all rectangles are drwan in one direction
### yDirection = "up" where to draw the rect. up or down or both (center alignment)
###    if constant then all rectangles are drwan in one direction
### EXAMPLE:
### plot(c(0, length(vect)+1), range(c(0, vect)), type="n")
### vect = c(1, 4, 5, 3, 3.5)
### easyRectUnivMany(1:length(vect), 0, .45, vect, xDirection="both", yDirection="up", col="gray", border=FALSE)
  howMany = max(length(startX), length(startY))
  if (length(startX)==1){startX = rep(startX, howMany)}
  if (length(startY)==1){startY = rep(startY, howMany)}
  if (length(xLen)==1){xLen = rep(xLen, howMany)}
  if (length(yLen)==1){yLen = rep(yLen, howMany)}
  if (length(xDirection)==1){xDirection = rep(xDirection, howMany)}
  if (length(yDirection)==1){yDirection = rep(yDirection, howMany)}
  for (i in 1:howMany){
    easyRectUniv(c(startX[i], startY[i]), xLen[i], yLen[i], xDirection[i], yDirection[i], color[i], ...)
  }
}

plotAbsValueMatrix = function(xCoor, yCoor, values, color, grid=FALSE, ...){
### absolute values are plotted at xCoor (from left to right starting from the left column)
###                            and yCoor (from bottom to top starting from the bottom row)
### bigger absolute values get bigger rectangles
### EXAMPLE:
###   set.seed(2)
###   values = matrix(round(rnorm(12)*2), nrow=3)
###   xCoor = 1:ncol(values)
###   yCoor = 1:nrow(values)
###   color = matrix("green", ncol=ncol(values), nrow=nrow(values))
###   color[values<0] = "red"
###   plot(range(xCoor)+c(-1,1), range(yCoor)+c(-1,1), type="n", ann=FALSE)
###   plotAbsValueMatrix(xCoor, yCoor, values, color, border=NA)

  if (length(xCoor) != ncol(values)) stop("error 1 in plotAbsValueMatrix\n")
  if (length(yCoor) != nrow(values)) stop("error 2 in plotAbsValueMatrix\n")
  if (any(dim(values) != dim(color))) stop("error 3 in plotAbsValueMatrix\n")
  startX = rep(xCoor, length(yCoor))
  startY = rep(yCoor[length(yCoor):1], each=length(xCoor))
  xLen = (abs(values)/max(abs(values))) * max(abs(diff(xCoor)))
  yLen = (abs(values)/max(abs(values))) * max(abs(diff(yCoor)))
  xDirection="both"
  yDirection="both"
  if (grid){
    abline(h = yCoor, lty=3, col="gray")
    abline(v = xCoor, lty=3, col="gray")
  }
  easyRectUnivMany(startX, startY, as.vector(t(xLen)), as.vector(t(yLen)), xDirection, yDirection, as.vector(t(color)), ...)
}
if(FALSE){
  set.seed(2)
  values = matrix(round(rnorm(120)*2), nrow=11)
  xCoor = 1:ncol(values)
  yCoor = 1:nrow(values)
  color = matrix("green", ncol=ncol(values), nrow=nrow(values))
  color[values<0] = "red"
  plot(range(xCoor)+c(-1,1), range(yCoor)+c(-1,1), type="n", ann=FALSE)
  plotAbsValueMatrix(xCoor, yCoor, values, color, border=NA, grid=TRUE)
}

printMissingPattern = function(datam, variables){
  ### example:
#   load(file = paste(currentVersonFolder, "/data/rda/datam.rda", sep=""))
#   names(datam)
#   covariates = c("patient_age", "patient_sex", "caregiverage", "caregiver_sex", "caregiver_type", "demo_d1", "demo_d1_binary", "demo_d2", "demo_d2_binary", "racelatino", "language", "wicstatus", "caregiveredu", "income", "stofhla", "stofhla_binary", "phlatpercentcorrect2mo", "site")
#   covariates %in% names(datam)
#   covar1 = covariates[covariates %in% names(datam)]
#   tmp = datam[,covar1]
#   printMissingPattern(tmp, covar1)
  datam = datam[, variables]
  missPerVar = sapply(datam, function(x){sum(is.na(x))})
  missmatr = is.na(as.matrix(datam)) + 0
  colnames(missmatr) = names(datam)
  missmatr = missmatr[, names(missPerVar[order(-missPerVar)])]
  ordervar = apply(missmatr, 1, paste, collapse="")
  ordervar
  missmatr = missmatr[order(ordervar), ]
  missmatr = missmatr[apply(missmatr, 1, paste, collapse="") != paste(rep(0, ncol(missmatr)), collapse=""), ]
  col = as.matrix(ifelse(missmatr==1, "gray", "white"), nrow=nrow(missmatr), ncol=ncol(missmatr))
  par(mar=c(5, 0, 5, 0))
  plot(range(1:ncol(missmatr)), range(1:nrow(missmatr)), type="n", axes=FALSE, xlab="", ylab="")
  axis(side=1, at=1:ncol(missmatr), labels=colnames(missmatr), las=2)
  axis(side=3, at=1:ncol(missmatr), labels=colnames(missmatr), las=2)
  plotAbsValueMatrix(1:ncol(missmatr), nrow(missmatr):1, missmatr, color=col, grid=FALSE)
}

boxWithWiskers = function(data, x, boxWidth, wiskWidth=boxWidth/2, jitterX=FALSE, jitterXPar=0, jitterY=FALSE, jitterYPar=0, boxCol="black", lineWidth=1, ...){
  ### data as a vector
  objectPar = quantile(data, c(.05, .25, .50, .75, .95))
  if (jitterX){
    xp = jitter(rep(x, length(data)), amount=jitterXPar)
  }else{
    xp = rep(x, length(data))
  }
  if (jitterX){
    yp = jitter(data, amount=jitterYPar)
  }else{
    yp = data
  }
  points(xp, yp, ...)
  rect(x - boxWidth/2, objectPar[2], x + boxWidth/2, objectPar[4], border=boxCol, lwd=lineWidth)
  segments(x - boxWidth/2, objectPar[3], x + boxWidth/2, objectPar[3], lwd=lineWidth*3, col=boxCol)
  segments(c(x - wiskWidth/2, x - wiskWidth/2), c(objectPar[1], objectPar[5]), c(x + wiskWidth/2, x + wiskWidth/2), c(objectPar[1], objectPar[5]), lwd=lineWidth, col=boxCol)
  segments(c(x, x), c(objectPar[1], objectPar[4]), c(x, x), c(objectPar[2], objectPar[5]), lty=2, col=boxCol, lwd=lineWidth)
}

categorizeIntervals = function(variable, intervals, categories, lessEqual){
  ### gives categories to different segments of the variable
  ### variable - a variable to be split into categories
  ### intervals - splitting points
  ### categories - what categories to assign (number of categories is number of splitting points plus one)
  ### lessEqual - if lessEqual=TRUE then <=, if lessEqual=FALSE then < (or >=)
  ### example: 
  ### variable =  seq(0, 1.7, .1);  intervals = c(.25, .7, 1.4); categories = c("less or equal .25", "less or equal .7", "less or equal 1.4", "greater than 1.4"); lessEqual=TRUE
  ### categorizeIntervals(variable, intervals, categories, lessEqual)
  ### variable =  seq(0, 1.7, .1);  intervals = c(.25, .7, 1.4); categories = c("less .25", "less .7", "less 1.4", "greater or equal than 1.4"); lessEqual=FALSE
  ### categorizeIntervals(variable, intervals, categories, lessEqual)

  if (length(categories) != length(intervals)+1) stop("Number of categories should be the number interval points (length(intervals)) plus one")
  tolerance = .Machine$double.eps ^ 0.5
  newVar = rep(categories[length(categories)], length(variable))
  for (i in length(intervals):1){
    if (lessEqual){
      newVar = ifelse(abs(variable-intervals[i]) < tolerance | variable < intervals[i], categories[i], newVar)
    }else{
      newVar = ifelse(variable < intervals[i] - tolerance, categories[i], newVar)
    }
  }
  newVar
}

sumquick = function(data, fileName, caption, varNames, strataName = "", overAll = FALSE, ...){
  if (strataName==""){
    overAll = FALSE
    test = FALSE
  }else{
    test = TRUE
    varNames = setdiff(varNames, strataName)
  }
  summ <- summary(as.formula(paste(strataName, "~", paste(varNames, collapse="+"))), data = data, method = "reverse", continuous = 1, overall = overAll, test = test, ...)
  trash = latex(summ, file = fileName, where = "!h", caption = caption, long = TRUE, exclude1 = FALSE, prmsd = TRUE, digits = 2, prtest = "P", npct = "both", size="tiny")
}


###-------------------------------------------------------------------------------------------
###-------------------------------R tricks --------------------------------------------------
###-------------------------------------------------------------------------------------------
myfactor = function(var, levels=NULL){
  if (is.null(levels)){
    return(factor(var))
  }
  if(!all(setdiff(unique(var), c(NA, NaN)) %in% levels)){
    cat(paste(charToRaw(setdiff(setdiff(unique(var), c(NA, NaN)), levels)), collapse="\",\""), "\n")
    stop("Values [\"", paste(setdiff(setdiff(unique(var), c(NA, NaN)), levels), collapse="\",\""), "\"] are not in levels [\"", paste(levels, collapse="\",\""), "\"]\n")
  }
  return(factor(var, levels=levels))
}

myasnumeric = function(var){
  ### safely converts to numeric
  ### (unlike as.numeric() it checks if some values cannot be converted to numeric
  ### and notifies about it). Also, when converting a factor variable to numeric
  ### it makes sure to convert it to character first.
  ### examples:
  if(FALSE){
    somestr = c("$3", "+5", "4.5", NA, "4,000")
    as.numeric(somestr)

    myasnumeric(somestr)

    somefactorstr = as.factor(c("0", "1", "2"))
    as.numeric(somefactorstr)

    myasnumeric(somefactorstr)

  }
  if (sum(is.na(var)) != sum(is.na(as.numeric(as.character(var))))){
    cat("Values [\"", paste(var[!is.na(var)&is.na(as.numeric(as.character(var)))], collapse="\",\""), "\"] cannot be converted to numeric\n\n")
    stop()
  }
  as.numeric(as.character(var))
}

compare = function(a, b){
###################### returns TRUE if A==B, and FALSE if A!=B. 
######################         if A=NA and B=1 will return FALSE
######################         if A=NA and B=NA will return TRUE
  c = a==b
  c[is.na(a) & !is.na(b)] = FALSE
  c[is.na(a) & is.na(b)] = TRUE
  c
}

positive = function(x){
### used by spline functions
  x[x <= 0]=0
  x
}

components = function(x, knots=NULL){
### gives sline components
  if (is.null(knots)) {knots = quantile(x, c(.1, .5, .9))}
  k = length(knots)
  res = matrix(NA, nrow=length(x), ncol=k-1)
  res[, 1] = x
  kd = (knots[k] - knots[1])^(2/3)
  for (j in 1:(k-2)){
    res[, j+1] = positive(((x-knots[j])/kd)^3) - positive(((x-knots[k-1])/kd)^3) * (knots[k]-knots[j])/(knots[k]-knots[k-1]) + positive(((x-knots[k])/kd)^3) * (knots[k-1]-knots[j])/(knots[k]-knots[k-1])
  }
  res
}

findLineParam = function(point1, point2){
  x1 = point1[1]; y1 = point1[2]
  x2 = point2[1]; y2 = point2[2]
  k = (y1-y2)/(x1-x2)
  b = y1 - x1*(y1-y2)/(x1-x2)
  c(k,b)
}

lineTrans = function(x, lineK_and_B){
  x*lineK_and_B[1]+lineK_and_B[2]
}

interp = function(x, point1, point2){
  kab = findLineParam(point1, point2)
  x*kab[1] + kab[2]
}

interpInASetOfPoints = function(pointA, Xs, Ys){
  ### takes a set of Xs and a corresponding set of Ys
  ### and interpolates a point pointA (returns its corresponding y)
  ### if pointA is out of the range of Xs returns NA
  ### Example:
  ###   Xs = c(1, 2, 4, 6, 7)
  ###   Ys = c(1.5, 1.0, 3.0, 3.0, 5.0)
  ###   pointA = c(3, 10,  1,  0)
  ###   > interpInASetOfPoints(pointA, Xs, Ys)
  ###   [1] 2.0  NA 1.5  NA
  if (length(Xs) != length(Ys)) stop("Xs and Ys should have the same length")
  lengthOfXs = length(Xs)
  Xs = Xs[order(Xs)]
  Ys = Ys[order(Xs)]
  res = rep(NA, length(pointA))
  for (i in 1:length(pointA)){
    eless = Xs <= pointA[i]
    emore = Xs >= pointA[i]
    index1 = sum(eless)
    index2 = lengthOfXs - sum(emore) + 1
    if (index1 != 0 & index1 != index2){
      x1 = Xs[index1]; y1 = Ys[index1]
      x2 = Xs[index2]; y2 = Ys[index2]
      k = (y1-y2)/(x1-x2)
      b = y1 - x1*(y1-y2)/(x1-x2)
      res[i] = b + k * pointA[i]
    }else{
      if (index1 == index2){
        res[i] = Ys[index1]
      }else{
        res[i] = NA ### return NA when pointA is less than the lowest Xs
      }
    }
  }
  res
}

# lookUpInMonotonVec = function(vector, value){
#   ### works on a numeric ordered vector with unique values
#   ### returns a matrix of c(index, index) of value in vector if vector contains value,
#   ### If there exist two values (v1, v2) in vector such that
#   ### v1< value <v2, then returns indices of these two vector elements
#   ### if value is greater than the max(vector) returns the index of the last (largest) element of vector and NA
#   ### if value is less than the min(vector) returns NA and the index of the first (smallest) element of vector
#   vector = as.numeric(vector)
#   value = as.numeric(value)
#   if (any(is.na(vector))) stop ("\"vector\" should not have missing values")
#   if (any(vector!=vector[order(vector)])) stop ("\"vector\" must be ordered")
#   if (length(vector)!=length(unique(vector))) stop ("All values of \"vector\" must be unique")
#   retMatr = matrix(NA, nrow = length(value), ncol=2)
#   for (i in 1:length(value)){
#     if (is.na(value[i])){
#       retMatr[i, ] = c(NA, NA)
#     }else{
#       matchIndex = match(value[i], vector)
#       if (!is.na(matchIndex)){
#         retMatr[i, ] = c(matchIndex, matchIndex)
#       }else{
#         if (value[i] < vector[1]){
#           retMatr[i, ] = c(NA, 1)
#         }else{
#           if (value[i] > vector[length(vector)]){
#             retMatr[i, ] = c(length(vector), NA)
#           }else{
#             secondIndex = match(TRUE, as.numeric(vector) > value[i])
#             retMatr[i, ] = c(secondIndex-1, secondIndex)
#           }
#         }
#       }
#     }
#   }
#   retMatr
# }
# splineValue = function(x, knots, betas){
# ### gives sline values
#   if (length(betas) != length(knots)-1) stop("length(betas) should be equal to number of knots minus 1")
#   k = length(knots)
#   kd = (knots[k] - knots[1])^(2/3)
#   allbetas = c(betas, rep(NA, 2))
#   allbetas[k] = sum(betas[2:(k-1)]*(knots[1:(k-2)] - knots[k]))/(knots[k] - knots[k-1])
#   allbetas[k+1] = sum(allbetas[2:(k-1)]*(knots[1:(k-2)] - knots[k-1]))/(knots[k-1] - knots[k])
#   cat(allbetas, "\n")
#   #allbetas[1]*x + sum(allbetas[2:(k+1)]*(x-knots)^3)
#   sumAll = allbetas[1]*x
#   for (i in 2:(k+1)){
#     sumAll = sumAll + positive(allbetas[i]*((x-knots[i-1])/kd)^3)
#   }
#   sumAll
# }
