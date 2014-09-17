
######################################## A set of functions allowing to run ANY model with splines and iteractions
######################################## This code is "tuned" to competing risk regression by Gray and Fine
######################################## but after some modification, it can be used with any regression
######################################## implemented in R as long as this regression has usual coefficient
######################################## and variance-covariance matrix
### Author: Svetlana Eden 

positive = function(x){
### used by spline functions
  x[x <= 0]=0
  x
}

splineKnots  = function(var, knotNumber, ...){
  if (knotNumber <= 2 | knotNumber >=8)stop("function splineKnots() is not implemented for knotNumber =", knotNumber)
  if (knotNumber == 3) varQuantiles = c(.1, .5, .9)
  if (knotNumber == 4) varQuantiles = c(.05, .35, .65, .95)
  if (knotNumber == 5) varQuantiles = c(.05, .275, .5, .725, .95)
  if (knotNumber == 6) varQuantiles = c(.05, .23, .41, .59, .77, .95)
  if (knotNumber == 7) varQuantiles = c(.025, .1833, .3417, .5, .6583, .8167, .975)
  return(quantile(var, probs=varQuantiles, na.rm=TRUE, ...))
}

components = function(x, knots=NULL){
### gives sline components
  if (is.null(knots)) {knots = splineKnots(x, 3)}
  if (length(knots) != length(unique(knots))){stop("To computed spline components need unique knots. Please change your knots: ", paste(knots, collapse=", "))}
  k = length(knots)
  res = matrix(NA, nrow=length(x), ncol=k-1)
  res[, 1] = x
  kd = (knots[k] - knots[1])^(2/3)
  for (j in 1:(k-2)){
    res[, j+1] = positive(((x-knots[j])/kd)^3) - positive(((x-knots[k-1])/kd)^3) * (knots[k]-knots[j])/(knots[k]-knots[k-1]) + positive(((x-knots[k])/kd)^3) * (knots[k-1]-knots[j])/(knots[k]-knots[k-1])
  }
  res
}

recursiveInt = function(interaction, rcsComponents, interTag = "*"){
  ### creates a list of variables needed to create an interaction term (whether it has rcs or not)
  ### example:
  ### interList = list(c("age", "lit", "tx"), c("dnt", "lit"))
  ### rcsComponents = list(age = c("age1", "age2"), lit = c("lit1", "lit2", "lit3"), tx = c("tx1", "tx2"), dnt = c("dnt"))
  ### recursiveInt(interaction = interList[[1]], rcsComponents)
  ### recursiveInt(interaction = interList[[2]], rcsComponents)
  intStrVec = c()
  currentRcsComponents = rcsComponents[[interaction[1]]]
  if (length(interaction) == 1){
    return(currentRcsComponents)
  }
  remainingInteractionVars = interaction[2:length(interaction)]
  for (n in currentRcsComponents){
    bunchOfSplineNames = recursiveInt(remainingInteractionVars, rcsComponents)
    intStrVec = c(intStrVec, paste(n, bunchOfSplineNames, sep="*"))
  }
  intStrVec
}

createOneVariable = function(sdata, variableFormula, splineKnotsStorage=NULL, rcsTag = "\'"){
  ### create One variable from a formula that can look in general like: "age'*lit''*tx" (interaction term)
  ### as you can see it has spline components of age and literacy, and also original treatment
  ### resulting variable is computed as a product of each interaction component 
  ###     Note: this function could be optimized in the following way:
  ###           it would create a component if it doesn't exist in sdata. But I chose not to do it
  ###           because I wanted the data to be updated just in case.
  variable = rep(1, nrow(sdata))
  ### check for interactions first:
  breakFormulaApart = unlist(strsplit(variableFormula, "\\*"))
  breakFormulaApart = gsub(" +", "", breakFormulaApart)
  for (vi in breakFormulaApart){
    ### now deal with splines if any
    originalVarName = gsub(rcsTag, "", vi)
    #cat("*******", variableFormula, vi, originalVarName, "\n")
    if (!(originalVarName %in% names(sdata))){
      stop("Variable ", originalVarName, " is needed to compute ", vi, " in function createOneComponent")
    }
    splineOrder = nchar(vi) - nchar(originalVarName) + 1
    if (splineOrder > 1){
      if (is.null(splineKnotsStorage) | !(originalVarName  %in% names(splineKnotsStorage))){stop("Need knot info for variable ", vi)}
      variable = variable * components(sdata[[originalVarName]], knots = splineKnotsStorage[[originalVarName]])[,splineOrder]
    }else{
      variable = variable * sdata[[vi]]
    }
  }
  variable
}

runCRRwithRCSandIteractions = function(crrdata, timeToEventVar, eventVar, failCodeValue,
                                       varNames, rcsVarNames = NULL, interactions = NULL,
                                       coefName = "coef", varCovarName = "var", rcsTag = "'", interTag = "\\*"){
  ### creates splines and interactions, runs crr, returns crr object with spline knots info,
  ###       and p-values for spline variables and for non linear components separately
  ### !!!!!!!! need to check: how it works with interactions (if it creates the correct variable,
  ###                          especially if the interaction is with spline variable)
  ### !!!!!!!! need to implement: p-values for interactions
  ### crrdata - data.frame with data
  ### timeToEventVar - name of the variable in crrdata that has time to event  
  ### eventVar - name of the variable in crrdata that has event
  ### failCodeValue - in crr this variable can have several codes
  ###    for example: 0 for censoring, 1 for recurrent AKI, and 2 for death
  ### varNames = variables that are included as linear terms
  ### rcsVarNames = variables that are modeled as splines
  ### interactions = interaction terms.
  if (FALSE){
    varNames = c("dnt", "age", "tx", "cam")
    rcsVarNames = list(age=3, lit=4)
    interactions = list(c("age", "lit", "tx"), c("dnt", "tx"))
    varNames = c("albbl", "age", "phos6", "realm")
    interactions = list(c("age", "phos6"))
    rcsVarNames = list(age=3, realm=3)
  }
  allVars = unique(c(varNames, names(rcsVarNames), unlist(interactions)))
  nonRCSvars = setdiff(allVars, names(rcsVarNames))
  ############## create spline variable names
  ############## create spline variable names
  repTag = function(x, tag = rcsTag){sapply(1:(x-1), function(x){paste(rep(tag, x-1), collapse="")})}
  RCSsuffix = sapply(rcsVarNames, repTag)
    if(class(RCSsuffix)!="list"){RCSsuffix = as.list(as.data.frame(RCSsuffix, stringsAsFactors=FALSE))}
  modelMatrVarForRCS = sapply(names(RCSsuffix), function(n){paste(n, RCSsuffix[[n]], sep="")})
    if(class(modelMatrVarForRCS)!="list"){modelMatrVarForRCS = as.list(as.data.frame(modelMatrVarForRCS, stringsAsFactors=FALSE))}
  modelMatrVarForNonRCS = sapply(nonRCSvars, function(n){n})
  #modelMatrVarForBoth = c(unlist(modelMatrVarForRCS), modelMatrVarForNonRCS)
  modelMatrVarForBoth = c(modelMatrVarForRCS, modelMatrVarForNonRCS)
  ############## create interaction variable names
  ############## create interaction variable names
  modelMatrVarForInteractions = sapply(interactions, function(x){recursiveInt(x, modelMatrVarForBoth )})
  if(class(modelMatrVarForInteractions) == "matrix"){
    modelMatrVarForInteractions = list(as.vector(modelMatrVarForInteractions))
  }
  if(length(modelMatrVarForInteractions)!=0){
    #names(modelMatrVarForInteractions) = paste("*", sapply(interactions, paste, collapse = " * "), "*", sep="")
    names(modelMatrVarForInteractions) = paste("*", sapply(interactions, paste, collapse = "*"), "*", sep="")
  }
  modelMatrVarForAll = unlist(c(modelMatrVarForBoth, modelMatrVarForInteractions))
  splineKnotsStorage = list()
  ############## create spline knots
  ############## create spline knots
  for (n in names(rcsVarNames)){
    splineKnotsStorage[[n]] = splineKnots(crrdata[[n]], rcsVarNames[[n]], TRUE)
  }
  ############## create a list of regression components to compute model effect
  ##############   for each original variable this list will contain
  ##############   all the spline components and all the interaction term components
  ##############   that relate to this variable, so all the coefficients are accounted for
  varsToGetToComputeEffect = list()
  for (n in allVars){
    varsToGetToComputeEffect[[n]] = c(n, as.character(modelMatrVarForRCS[[n]]))
    if(length(modelMatrVarForInteractions)!=0){
      isThereInteraction = grep(paste("*", n, "*", sep=""), names(modelMatrVarForInteractions), value=TRUE)  
      if (length(isThereInteraction)!=0){
        varsToGetToComputeEffect[[n]] = c(varsToGetToComputeEffect[[n]], unlist(sapply(isThereInteraction, function(n)modelMatrVarForInteractions[[n]])))
      }
    }
    varsToGetToComputeEffect[[n]] = unique(varsToGetToComputeEffect[[n]])
  }
  ############## create actual interactions and spline variables
  ############## create actual interactions and spline variables
  for (n in setdiff(modelMatrVarForAll, names(crrdata))){
    cat("creating variable", n, "\n")
    crrdata[[n]] = createOneVariable(crrdata, variableFormula = n, splineKnotsStorage, rcsTag = rcsTag)
  }
  cat("All vars", modelMatrVarForAll, "\n")
  if(TRUE){
    crrModel = crr(crrdata[[timeToEventVar]], crrdata[[eventVar]], crrdata[, modelMatrVarForAll], failcode = failCodeValue)
  }else{
    form = paste(timeToEventVar, "~", paste(setdiff(modelMatrVarForAll, timeToEventVar), collapse="+"))
    nonMissingRows = apply(crrdata[, c(modelMatrVarForAll, timeToEventVar)], 1, function(x)(all(!is.na(x))))
    tmp = crrdata[nonMissingRows, c(modelMatrVarForAll, timeToEventVar)]
    crrModel = ols(as.formula(form), data = tmp)
    names(crrModel[[coefName]]) = gsub(" +", "", names(crrModel[[coefName]]))
  }
  ############## name the model var-covar matrix:
  ############## name the model var-covar matrix:
  if(length(crrModel[[coefName]]) == 1){names(crrModel[[coefName]]) = modelMatrVarForAll[1]}
  colnames(crrModel[[varCovarName]]) = names(crrModel[[coefName]])
  rownames(crrModel[[varCovarName]]) = names(crrModel[[coefName]])
  ############## wrap the model in an object:
  ############## wrap the model in an object:
  returnObject = list(crrObject = crrModel)
  ############## p-values for each variable including splines:
  ############## p-values for each variable including splines:
  pValues = rep(NA, length(varsToGetToComputeEffect))
  names(pValues) = names(varsToGetToComputeEffect)
  for (n in names(varsToGetToComputeEffect)){
    pValues[n] = waldPvalueFromGray(crrModel[[coefName]], crrModel[[varCovarName]], kk = varsToGetToComputeEffect[[n]])["p"]
  }
  interPValues = rep(NA, length(modelMatrVarForInteractions))
  names(interPValues) = names(modelMatrVarForInteractions)
  for (n in names(modelMatrVarForInteractions)){
    interPValues[n] = waldPvalueFromGray(crrModel[[coefName]], crrModel[[varCovarName]], kk = modelMatrVarForInteractions[[n]])["p"]
  }
  ############## add all the necessary info to the crr object:
  ############## add all the necessary info to the crr object:
  returnObject[["varsToGetToComputeEffect"]] = varsToGetToComputeEffect
  returnObject[["modelMatrVarForInteractions"]] = modelMatrVarForInteractions
  returnObject[["splineKnotsStorage"]] = splineKnotsStorage
  returnObject[["pValues"]] = pValues
  returnObject[["interPValues"]] = interPValues
  returnObject[["modelResNames"]] = list(coefName = coefName, varCovarName = varCovarName, rcsTag = rcsTag)
  returnObject
}

  waldPvalueFromGray <- function (coefficients, varCovarMatrix, kk){  ## multiparameter Wald test
    w <- coefficients[kk]
    v <- varCovarMatrix[kk, kk]
    stat <- sum(w * solve(v, w)) ## solves the equation ‘v %*% x = w’ for ‘x’
    df <- length(w)
    c(stat = stat, df = df, p = 1 - pchisq(stat, df))
  }

predictCRR = function(crrWithSplinesKnots, crrdata){
  coefVarName = crrWithSplinesKnots$modelResNames$coefName
  allVars = names(crrWithSplinesKnots$crrObject[[coefVarName]])
  for (n in setdiff(allVars, names(crrdata))){
    cat("creating variable", n, "\n")
    crrdata[[n]] = createOneVariable(crrdata, variableFormula = n, crrWithSplinesKnots$splineKnotsStorage, rcsTag = crrWithSplinesKnots$modelResNames$rcsTag)
  }
  predict(crrWithSplinesKnots$crrObject, cov1 = as.matrix(crrdata[, allVars]))
}

effectsCRR = function(crrWithSplinesKnots, varName, effectRange, varValuesForInteractions=NULL, exponentiateEffects=FALSE){
  ####### Example:
  ###     covariates = c("age", "Female", "MostRecentPreadmitOutptCr", "MaxInptCr", "RaceBlack", "RaceWhite")
  ###     crrObjWithSplines = runCRRwithRCSandIteractions(crrdata = maindata, timeToEventVar = "timeToRecurAKIOrHospDeathFromDis365",
  ###          eventVar = "hadRecurAKIOrHospDeathFromDis365", failCodeValue=1, varNames = covariates,
  ###          rcsVarNames = list(age = 4, MostRecentPreadmitOutptCr=4, MaxInptCr=4))
  ###     effectsCRR(crrModel = crrObjWithSplines, varName="cr", effectRange=c(0, 1))
  ### returns the effect in the format of Hazard Ratio (HR)
  ### varName - a name of a covariate
  ### effectRange - is in general a matrix of 2 columns and as many rows as effects are needed: from - to
  ### varValuesForInteractions - is a vector for variables included in interaction terms of the varName
  ###   for example if interaction of age is with tx and lit and if varName = "age"
  ###     then varValuesForInteractions = list(tx=1, lit=17)
  ############## multidimentional wald test, like chunk test or exactly chunk test
  ############## multidimentional wald test, like chunk test or exactly chunk test
  ############## multidimentional wald test, like chunk test or exactly chunk test
  ############## check dimension of effectRange
  ############## check dimension of effectRange
  if (is.null(dim(effectRange))){
    if (length(effectRange)!=2) stop("variable effectRange has to be either a vector of length 2 or an Nx2 matrix")
      effectRange = matrix(effectRange, nrow = 1)
  }
  if (dim(effectRange)[2]!=2) stop("variable effectRange has to be either a vector of length 2 or an Nx2 matrix")
  ############## check variables and names variance-covariance matrix
  ############## check variables and names variance-covariance matrix
  crrModel = crrWithSplinesKnots[["crrObject"]]
  coefName = crrWithSplinesKnots[["modelResNames"]][["coefName"]]
  varCovarName = crrWithSplinesKnots[["modelResNames"]][["varCovarName"]]
  varNames = names(crrModel[[coefName]])
  if (!(varName %in% varNames)) stop("Variable ", varName, " was not found in the model -  cannot provide the effect.")
  relatedComponents = crrWithSplinesKnots[["varsToGetToComputeEffect"]][[varName]]
  ################# create a matrix of "from"
  ################# create a matrix of "from"
  ################# create a data.frame that can be used in createOneVariable()
  ################# create a data.frame that can be used in createOneVariable()
  compdata = data.frame(effectRange[,1])
  names(compdata) = varName
  for (n in names(varValuesForInteractions)){
    compdata[[n]] = varValuesForInteractions[[n]]
  }
  fromMatrix = matrix(NA, nrow=dim(effectRange)[1], ncol=length(relatedComponents))
  colnames(fromMatrix) = relatedComponents
  for (n in relatedComponents){
    fromMatrix[, n] = createOneVariable(compdata, n, crrWithSplinesKnots[["splineKnotsStorage"]], crrWithSplinesKnots[["modelResNames"]][["rcsTag"]])
  }
  ################# create a matrix of "to"
  ################# create a matrix of "to"
  ################# create a data.frame that can be used in createOneVariable()
  ################# create a data.frame that can be used in createOneVariable()
  compdata = data.frame(effectRange[,2])
  names(compdata) = varName
  for (n in names(varValuesForInteractions)){
    compdata[[n]] = varValuesForInteractions[[n]]
  }
  toMatrix = matrix(NA, nrow=dim(effectRange)[1], ncol=length(relatedComponents))
  colnames(toMatrix) = relatedComponents
  for (n in relatedComponents){
    toMatrix[, n] = createOneVariable(compdata, n, crrWithSplinesKnots[["splineKnotsStorage"]], crrWithSplinesKnots[["modelResNames"]][["rcsTag"]])
  }
  ################# compute total effect
  ################# compute total effect
  diffMatrix = as.matrix(as.data.frame(toMatrix)[, relatedComponents] - as.data.frame(fromMatrix)[,relatedComponents]) ### this painful conversion is necessary in case when nrow(toMatrix) is 1, and R CHOOSES (without you suspecting it) to convert matrices to vectors ...
  betas = matrix(crrModel[[coefName]][relatedComponents], ncol=1)
  rownames(betas) = relatedComponents
  effect = diffMatrix %*% betas
  vars = crrModel[[varCovarName]][relatedComponents, relatedComponents]
  ################# compute total SE
  ################# compute total SE
  totalSE = sqrt(diag(diffMatrix %*% vars %*% t(diffMatrix)))
  CIsLower = effect + -1*1.96*totalSE
  CIsUpper = effect + 1.96*totalSE
  pval = waldPvalueFromGray(crrModel[[coefName]], crrModel[[varCovarName]], relatedComponents)["p"]
  effectPval = 1 - pchisq((effect/totalSE)^2, 1) ### or: effectPval = 2*pnorm(abs(effect/totalSE)*(-1))
  res = data.frame(rep(varName, nrow(effectRange)), effectRange[,1], effectRange[,2], exp(effect), exp(CIsLower), exp(CIsUpper), effectPval, rep(pval, nrow(effectRange)))
  names(res) = c("Variable", "From", "To", "Effect", "Lower", "Upper", "EffectPvalue", "OverallPvalue")
  if(exponentiateEffects){
    for (n in c("Effect", "Lower", "Upper")){
      res[[n]] = exp(res[[n]])
    }
  }
  res
}

effectsForDiffVarInterValues = function(crrWithSplinesKnots, effectRangeFrom, effectRangeTo, exponentiateEffects = FALSE){
  ### computes effect difference for ANY covariates levels (one of the things that "contrast" can do)
  ####### Example:
  ###     covariates = c("age", "Female", "MostRecentPreadmitOutptCr", "MaxInptCr", "RaceBlack", "RaceWhite")
  ###     crrObjWithSplines = runCRRwithRCSandIteractions(crrdata = maindata, timeToEventVar = "timeToRecurAKIOrHospDeathFromDis365",
  ###          eventVar = "hadRecurAKIOrHospDeathFromDis365", failCodeValue=1, varNames = covariates,
  ###          rcsVarNames = list(age = 4, MostRecentPreadmitOutptCr=4, MaxInptCr=4))
  ###     effectsForDiffVarInterValues(crrModel = crrObjWithSplines, effectRangeFrom=list(MaxInptCr=1, Female=0), effectRangeTo=list(MaxInptCr=2, Female=1)))
  ### effectRangeFrom - is a data.frame of all "from" values
  ### effectRangeTo   - is a data.frame of all "to" values
  ### effectRangeFrom and effectRangeTo should have the same elements
  ###   for example MaxInptCr and Female
  ###   they also have to contain all the variable values that are needed to compute interaction with the variable specified in them
  ############## multidimentional wald test, like chunk test or exactly chunk test
  ############## multidimentional wald test, like chunk test or exactly chunk test
  ############## multidimentional wald test, like chunk test or exactly chunk test
  ############## check dimension of effectRange
  ############## check dimension of effectRange
  if (class(effectRangeFrom)!= "data.frame" | class(effectRangeTo)!= "data.frame"){stop("effectRangeFrom and effectRangeTo should be data.frame's")}
  if (!setequal(colnames(effectRangeFrom), colnames(effectRangeTo))) stop("effectRangeFrom and effectRangeTo should be have the same variables specified in them")
  ############## check variables and names variance-covariance matrix
  ############## check variables and names variance-covariance matrix
  crrModel = crrWithSplinesKnots[["crrObject"]]
  coefName = crrWithSplinesKnots[["modelResNames"]][["coefName"]]
  varCovarName = crrWithSplinesKnots[["modelResNames"]][["varCovarName"]]
  if (any(!(colnames(effectRangeFrom) %in% names(crrModel[[coefName]])))) stop("All variables specified in effectRangeFrom should be found in the model.")
  relatedComponents = unique(unlist(crrWithSplinesKnots[["varsToGetToComputeEffect"]][names(effectRangeFrom)]))
  ################# create a matrix of "from"
  ################# create a matrix of "from"
  fromMatrix = effectRangeFrom
  for (n in relatedComponents){
    fromMatrix[, n] = createOneVariable(effectRangeFrom, n, crrWithSplinesKnots[["splineKnotsStorage"]], crrWithSplinesKnots[["modelResNames"]][["rcsTag"]])
  }
  ################# create a matrix of "to"
  ################# create a matrix of "to"
  toMatrix = effectRangeTo
  for (n in relatedComponents){
    toMatrix[, n] = createOneVariable(effectRangeTo, n, crrWithSplinesKnots[["splineKnotsStorage"]], crrWithSplinesKnots[["modelResNames"]][["rcsTag"]])
  }
  ################# compute total effect
  ################# compute total effect
  diffMatrix = as.matrix((toMatrix - fromMatrix)[, relatedComponents], nrow=1)
  betas = matrix(crrModel[[coefName]][relatedComponents], ncol=1)
  rownames(betas) = relatedComponents
  effect = diffMatrix %*% betas
  vars = crrModel[[varCovarName]][relatedComponents, relatedComponents]
  ################# compute total SE
  ################# compute total SE
  totalSE = sqrt(diag(diffMatrix %*% vars %*% t(diffMatrix)))
  CIsLower = effect + -1*1.96*totalSE
  CIsUpper = effect + 1.96*totalSE
  pval = waldPvalueFromGray(crrModel[[coefName]], crrModel[[varCovarName]], relatedComponents)["p"]
  effectPval = 1 - pchisq((effect/totalSE)^2, 1) ### or: effectPval = 2*pnorm(abs(effect/totalSE)*(-1))
  res = data.frame(effect, CIsLower, CIsUpper, effectPval, rep(pval, nrow(effectRangeFrom)),
          stringsAsFactors = FALSE)
  names(res) = c("Effect", "Lower", "Upper", "EffectPvalue", "OverallPvalue")
  if(exponentiateEffects){
    for (n in c("Effect", "Lower", "Upper")){
      res[[n]] = exp(res[[n]])
    }
  }
  res
}

summaryForCRR = function(crrWithSplinesKnots, varRangeList, varValuesForInteractions=NULL){
  ### ATTENTION: computes effects for variable difference one at a time: for example:
  ###   age 20 vs 30, with the contstant interaction value of tx (say 0)
  ###      in other words it doesn't compute age=20, tx=0   vs   age=30, tx=1
  ### varRangeList- a vector of variables (spline components are considered as one variable)
  ###   if a variable was model with splines the function will need knots for each variable
  ### effectRange - a vector of length 2: from - to
  ### knots - if the variable is a spline, knots must be supplied to calculate the effect
  ### 
  ### summaryForCRR(crrWithSplinesKnots, varRangeList=list(age = c(50, 60)))
  res = data.frame(Variable = c(), From = c(), To = c(), Effect = c(), Lower = c(), Upper = c(), EffectPvalue = c(), OverallPvalue = c(),
                   stringsAsFactors = FALSE)
  for (v in names(varRangeList)){
    #cat(v, "____________", setdiff(names(varValuesForInteractions), v), "\n")
    res = rbind(res, effectsCRR(crrWithSplinesKnots, v, varRangeList[[v]], varValuesForInteractions=varValuesForInteractions[setdiff(names(varValuesForInteractions), v)]))
  }
  leftOverVars = setdiff(names(crrWithSplinesKnots[["varsToGetToComputeEffect"]]), names(varRangeList))
  if(!is.null(varValuesForInteractions)){
    cat("Effects with interactions were adjusted for:\n"); print(varValuesForInteractions)
  }
  if(length(leftOverVars)){
    cat("The following variables weren't summarized:", paste(leftOverVars, collapse=", "), "\n")
  }
  rownames(res) = 1:nrow(res)
  res$Variable = as.character(res$Variable)
  res
}

################################# example of how to run this ###############################
################################# example of how to run this ###############################
################################# example of how to run this ###############################
################################# example of how to run this ###############################
################################# example of how to run this ###############################
################################# example of how to run this ###############################
################################# example of how to run this ###############################

### example for linear regression
### example for linear regression
### example for linear regression
if(FALSE){
  source("skeModel.R")
  load("main.rda")
  varNames = c("age", "phos6", "realm")
  interactions = list(c("age", "phos6"))
  rcsVarNames = list(age=3, realm=3)
  allVarMissing = sapply(main, function(x){sum(!is.na(x))<2})
  dd = datadist(main[, setdiff(names(main), names(allVarMissing[allVarMissing]))])
  options(datadist = "dd")
  #crrdata = main; timeToEventVar = "albbl"; eventVar = "albbl"; failCodeValue=1; varNames = c("age", "phos6", "realm"); rcsVarNames = rcsVarNames; rcsTag = "9"; coefName = "coefficients"; varCovarName = "var"; interactions = interactions
  mymod = runCRRwithRCSandIteractions(crrdata = main, timeToEventVar = "albbl",
                              eventVar = "albbl", failCodeValue=1,
                              interactions = interactions,
                              varNames = varNames, rcsVarNames = rcsVarNames, rcsTag = "9",
                              coefName = "coefficients", varCovarName = "var")
  mymod$crrObject
  fmodel = ols(albbl ~ rcs(age, 3)*phos6 + rcs(realm, 3), data = crrdata)

  #crrWithSplinesKnots = mymod; varName="age"; effectRange=c(20, 25); varValuesForInteractions=NULL
  source("skeModel.R")
  effectsCRR(crrWithSplinesKnots = mymod, varName="age", effectRange=c(20, 25), varValuesForInteractions=list(phos6 = 3))
  summaryForCRR(crrWithSplinesKnots = mymod, varRangeList = list(age=matrix(c(20, 25, 40, 50), byrow=TRUE, ncol=2), phos6=matrix(c(3, 6, 6, 9), byrow=TRUE, ncol=2), realm=c(0, 17)), varValuesForInteractions=list(age=40, phos6=3))
  summary(fmodel, age = c(20, 25), phos=c(3, 6))

  effectsForDiffVarInterValues(crrWithSplinesKnots=mymod, effectRangeFrom=data.frame(age = c(25, 25), phos6 = c(3, 3)), effectRangeTo=data.frame(age = c(20, 20), phos6 = c(3, 7)))
  contrast(fmodel, list(age = 20, phos6 = 3), list(age = 25, phos6 = 3))
  contrast(fmodel, list(age = 20, phos6 = 7), list(age = 25, phos6 = 3))

}

collapseCRRModelsAfterImputation = function(models){
   ### take several models (the only difference b/w them is in imputed values of missing vairables)
   ### and effects and covariance matrix based on:
   ### The function implies that betas are in models[[i]]$coefficients, and covariance matrix isi in the models[[i]]$var
   ### REF: Statistical Analysis with Missing Data, Second edition, Roderick J.A. Little, Donald B. Rubin, pages 85-87, 211.
   numOfImp = length(models)
   meanBetas = rep(0, length(models[[1]]$crrObject$coef))
   for (i in 1:numOfImp){
     meanBetas = meanBetas + models[[i]]$crrObject$coef
   }
   meanBetas = meanBetas/numOfImp

   withinVar = matrix(0, nrow=nrow(models[[1]]$crrObject$var), ncol=ncol(models[[1]]$crrObject$var))
   betweenVar = matrix(0, nrow=nrow(models[[1]]$crrObject$var), ncol=ncol(models[[1]]$crrObject$var))
   for (i in 1:numOfImp){
     withinVar = withinVar + models[[i]]$crrObject$var
     betweenVar = betweenVar + matrix(models[[i]]$crrObject$coef-meanBetas, ncol=1) %*% matrix(models[[i]]$crrObject$coef-meanBetas, nrow=1)
   }
   withinVar = withinVar/numOfImp
   if(numOfImp - 1 != 0){
     betweenVar = betweenVar/(numOfImp - 1)
   }
   totalVar = withinVar + ((numOfImp + 1)/numOfImp)*betweenVar

   newObject = models[[1]]$crrObject
   coefName = models[[1]]$modelResNames$coefName
   varCovarName = models[[1]]$modelResNames$varCovarName
   newObject[[coefName]] = meanBetas
   newObject[[varCovarName]] = totalVar

   ############## recomputed p-values:
   ############## recomputed p-values:
   ############## recomputed p-values:
   varsToGetToComputeEffect = models[[1]]$varsToGetToComputeEffect
   modelMatrVarForInteractions = models[[1]]$modelMatrVarForInteractions
   pValues = rep(NA, length(varsToGetToComputeEffect))
   names(pValues) = names(varsToGetToComputeEffect)
   for (n in names(varsToGetToComputeEffect)){
     pValues[n] = waldPvalueFromGray(newObject[[coefName]], newObject[[varCovarName]], kk = varsToGetToComputeEffect[[n]])["p"]
   }
   interPValues = rep(NA, length(modelMatrVarForInteractions))
   names(interPValues) = names(modelMatrVarForInteractions)
   for (n in names(modelMatrVarForInteractions)){
     interPValues[n] = waldPvalueFromGray(newObject[[coefName]], newObject[[varCovarName]], kk = modelMatrVarForInteractions[[n]])["p"]
   }
   modifiedModel = models[[1]]
   modifiedModel$crrObject = newObject
   modifiedModel$pValues = pValues
   modifiedModel$interPValues = interPValues
   modifiedModel
}

