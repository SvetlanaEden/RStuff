
####################### author: Svetlana K. Eden, svetlana.eden@vanderbilt.edu
################################################

################################################
################################################
################################################
### 
### These set of functions built a latex table for a matrix.
### As an application, it is suitable for two-by-two tables for sensitivity and specificity.
###

dbbackslash = function(){
  "\\\\"
}

cline = function(from, to){
  paste("\\cline{", from, "-", to,"}", sep="")
}
 
hline = function(n=1){
  paste(rep("\\hline", n), collapse="")
}

firsthead = function(caption, name, width1, width2, align="c", colnames){
  cr="\n"
  mlt = "\\multicolumn"
  nameSep = "} & {"
  str = paste("\\caption{",caption,"}", dbbackslash(), cr, sep=" ")
  str = paste(str, mlt, "{", width1, "}", "{", align, "}{} & ", mlt, "{", width2, "}", "{", align, "}{",name,"}", dbbackslash(), cr, cline(width1+1, width1+width2), cr, sep="")
  str = paste(str, mlt, "{", width1, "}", "{", align, "}{} & {", paste(colnames, collapse=nameSep), "}", dbbackslash(), cr, cline(width1+1, width1+width2), cr, sep="")
  str = paste(str, "\\endfirsthead", cr, sep="")
  str
}

nexthead = function(name, width1, width2, align="c", colnames){
  cr="\n"
  mlt = "\\multicolumn"
  nameSep = "} & {"
  str = paste("\\caption[]{(continued)}", dbbackslash(), hline(), cr, sep=" ")
  str = paste(str, mlt, "{", width1, "}", "{", align, "}{} & ", mlt, "{", width2, "}", "{", align, "|}{",name,"}", dbbackslash(), cr, cline(width1+1, width1+width2), cr, sep="")
  str = paste(str, mlt, "{", width1, "}", "{", align, "}{} & {", paste(colnames, collapse=nameSep), "}", dbbackslash(), cr, cline(width1+1, width1+width2), cr, sep="")
  str = paste(str, "\\endhead", cr, sep="")
  str
}

footnote = function(nextfootnote="Continued on next page", lastfootnote="", width, align="r"){
  cr="\n"
  mlt = "\\multicolumn"
  italbold = "\\bfseries\\em"
  str = paste(mlt, "{", width, "}{", align, "}{",italbold," ",nextfootnote,"}", dbbackslash(), cr, sep="")
  str = paste(str, "\\endfoot", cr)
  str = paste(str, mlt, "{", width, "}{", align, "}{",italbold," ",lastfootnote,"}", dbbackslash(), cr, sep="")
  str = paste(str, "\\endlastfoot", cr)
  str
}

multirow = function(nrows, name, namewidth){
  cr = "\n"
  str = paste("\\multirow", "{", nrows, "}{", namewidth, "pt}{",name,"}", cr, sep="")
  str
}

processRowsNByN = function(rownames, datamatrix){
  if (length(rownames) != nrow(datamatrix)) stop("function processRows: the length of 'rownames' has to be equal to the number of rows of the 'matrix'")
  cr="\n"
  am= " & "
  width1 = 2
  str = ""
  for (i in 1:nrow(datamatrix)){
    if (i != nrow(datamatrix)){
      str = paste(str, am, rownames[i], am, paste(datamatrix[i,], collapse=am), dbbackslash(), cr, cline(width1, width1+ncol(datamatrix)), cr)
    }else{
      str = paste(str, am, rownames[i], am, paste(datamatrix[i,], collapse=am), dbbackslash(), cr)
    }
  }
  str
}
#cat(processRows(c("iooio","mmmm"), data))


matrixTable <- function(file="", longcaption = "", fontSize="small", placement=NULL,
                   data, vertName, vertNames=rownames(data), horiName, horiNames=colnames(data), vertNameWidth = 30,
                   toLatexChar=TRUE, marker=NULL, append=FALSE){
  ### file - file name, if file="", then ... (useful for sweave)
  ### longcaption - caption for the table in the document
  ### caption - not in use: caption for the table in the table of contents
  ### fontsize - font size
  ### placement - if you'd like to put your table in the very spot it is in the .tex file then placement ="h"
  ### data - "matrix" or "table"
  ### vertName - what is displayed in the vertical orientation (example: "Gender")
  ### vertNames - the names of the categories, displayed in the vertical orien. (ex.: "Men", "Women")
  ### horiName - what is displayed in the horizontal orientation (example: "Smoking status")
  ### horiNames - the names of the categories, displayed in the horizontal orien. (ex.: "Yes", "No")
  ### toLatexChar - convert the content of the table into latex-text characters ? (if "data" is numeric then only names are converted)
  ### marker - if not NULL gives the label for the table
  ### append - whether to append to the "file" or not.
  # not in use: rightVars=c(), centerVars=c(), fixedColVars=c(), fixedColWdths=c(), pallet=NULL, pattern=NULL, landscape=FALSE
  
  longtable=TRUE
  fontSizes <- c("tiny","scriptsize","footnotesize","small",
                "normalsize","large","Large","LARGE","huge","Huge")
  
  #data = matrix(c(1,2,3,400,500,600), ncol=3, nrow=2, byrow=TRUE)
  
  ###-------------------getting data parameters:
  transpose=FALSE  ### if transpose is TRUE, the table is represented sideways (as a transposed matrix)
  if (class(data)=="table"){
    data = as.matrix(data)
  }
  if (transpose) {data = t(data)}
  nrows=nrow(data)
  ncols=ncol(data)
  
  
  ###-------------------column formats
  verticalHeadWidth = 2
  horizontHeadWidth = ncols
  multirowWidth = nrows
  cr="\n"
  colSep="|"
  colFormats = c("c", "l", "r", "pt")
  names(colFormats) = c("central", "left", "right", "fixedwidth")
  majorColFormatStr = paste("\\begin{longtable}{", paste(rep(colFormats["left"], verticalHeadWidth+ncols), collapse=colSep),"}", cr, sep="")
  
  ###-------------------headers
  cat(majorColFormatStr, file=file, append=FALSE)
  cat(firsthead(caption=longcaption, name=horiName, width1=verticalHeadWidth, width2=horizontHeadWidth, colnames=horiNames), file=file, append=TRUE)
  cat(nexthead(name=horiName, width1=verticalHeadWidth, width2=horizontHeadWidth, colnames=horiNames), file=file, append=TRUE)
  cat(footnote(width = verticalHeadWidth + horizontHeadWidth), file=file, append=TRUE)
  cat(multirow(nrows, vertName, vertNameWidth), file=file, append=TRUE)
  #cat(paste(vertNames, collapse=","), "\n")
  cat(processRowsNByN(vertNames, data), file=file, append=TRUE)
  cat(paste("\\end{longtable}", cr, sep=""), file=file, append=TRUE)
}

#data = matrix(c(1,2,3,400,500,600), ncol=3, nrow=2, byrow=TRUE)
#rownames(data) = c("rname1", "rname2")
#colnames(data) = c("cname1", "cname2", "cname2")
#matrixTable(file="tmp.tex", longcaption = "Table Example", #data=data, vertName="TRUE", vertNames=rownames(data), #horiName="TEST", horiNames=colnames(data), vertNameWidth = 30)


#data = matrix(c(1,2,3,400,500,600, .1, .2, .3), ncol=3, #byrow=TRUE)
#rownames(data) = c("rname1", "rname2", "rname3")
#colnames(data) = c("cname1", "cname2", "cname2")
#matrixTable(file="tmp.tex", longcaption = "Table Example", #data=data, vertName="TRUE", vertNames=rownames(data), #horiName="TEST", horiNames=colnames(data), vertNameWidth = 30)

################################################
################################################
################################################

####################### function charToStr
####################### replaces a character with a word. Can operate on a character vector
####################### Arguments: character vector "words" where replacement should be done
#######################            character vector "charToWord" with names as characters
#######################            to be replaced and
#######################            elements as words to replace these characters with.
charToStr <- function(words, charToWord){
  rawCharRange=as.raw(32:126)
  wordNames=charToRaw(paste(names(charToWord), collapse=""))
  if (!is.vector(charToWord)){
    stop("Argument 'charToWord' has to be a vector.")
  }
  if (typeof(charToWord)!="character"){
    stop("Argument 'charToWord' has to be a character type.")
  }
  if (length(wordNames)!=length(charToWord) | !all(wordNames %in% rawCharRange)){
    stop("The names of 'charToWord' have to be non-empty single characters.")
  }
  indBase = as.numeric(rawCharRange[1])-1
  notChar = as.raw(30)
  charMatrT = matrix(rawCharRange, nrow=max(nchar(charToWord)), ncol=length(rawCharRange))
  for (i in 1:length(rawCharRange)){
    rawChar = rawCharRange[i]
    ind = as.numeric(rawChar)-indBase
    if (rawToChar(rawChar) %in% names(charToWord)){
      charMatrT[,ind] = c(charToRaw(charToWord[rawToChar(rawChar)]),
                          rep(notChar, nrow(charMatrT)-nchar(charToWord[rawToChar(rawChar)])))
    }else{
      charMatrT[,ind] = c(rawChar, rep(notChar, nrow(charMatrT)-1))
    }
  }
  modifWords = rep(NA, length(words))
  for (i in 1:length(words)){
    ind = as.numeric(charToRaw(words[i]))-indBase
    modifWords[i] = rawToChar(charMatrT[,ind][charMatrT[,ind]!= notChar])
  }
  modifWords
}

if(FALSE){
####################### function latexTextModeNoRegExpr
####################### replaces special latex characters with words
#######################          that look like regular characters in Latex
####################### Arguments: character vector "words" where replacement should be done
#######################
latexTextModeNoRegExpr <- function(words){
  toLatexChar = c("#"="\\#", "$"="\\$", "%"="\\%", "^"="\\verb*+^+", "_"="\\_",
        "{"="\\{", "}"="\\}", "~"="\\verb*+~+", "&"="\\&", "\\"="$\\backslash$",
                                                "<"="$<$", ">"="$>$", "|"="$|$")
  charToStr(words, toLatexChar)
}
}

####################### function latexTextMode
####################### replaces special latex characters with words
#######################          that look like regular characters in Latex using regular expr.
####################### Arguments: character vector "words" where replacement should be done
#######################
latexTextMode <- function(words){
  ### put double back slash in frount of these ${}#%_&
  newWords = gsub("([${}#%_&])", "\\\\\\1", words)
  ### put dollar sign around < or >
  newWords = gsub("([<|>])", "$\\1$", newWords)
  ### if you see a double backslash not followed by one of these ${}#%_& , replace it with $\\\\backslash$
  newWords = gsub("\\\\(?!{|}|#|%|_|&|\\$)", "$\\\\backslash$", newWords, perl=TRUE)
  ### replace ~ or ^ with \\verb*+~+ or \\verb*+^+
  newWords = gsub("([~^])", "\\\\verb*+\\1+", newWords)
  newWords
}

zebraTable <- function(file="",
                      longtable=TRUE, landscape=FALSE,
                      caption = "", fontSize="small", placement=NULL,
                      dataframe, zebraPattern="none", by=names(dataframe)[1], orderGroups=FALSE,
                      colNames = names(dataframe),
                      vars =names(dataframe),
                      rightVars=c(), centerVars=c(), fixedColVars=c(), fixedColWdths=c(),
                      toLatexChar=TRUE, pallet=NULL, pattern=NULL,
                      marker=NULL, append=FALSE) {
  ### file: if file=="" the output is redirected to standard output
  ###       if file is character vector of length 1 the output is flushed into
  ###       a file with this character string as a name
  ###       file can be an open connection.
  ### append: logical, if TRUE all the out put will be appended to a file,
  ###         if the file is given as character string
  ### marker: a character string, if you want to reference this table in your latex document
  ### toLatexChar: logical, if TRUE all latex special symbols will be converted
  ###              to latex text symbols in the entire table.

  #internal constants and functions definitions
  fontSizes <- c("tiny","scriptsize","footnotesize","small",
                 "normalsize","large","Large","LARGE","huge","Huge")
  zebraPatterns <- c("none", "plain", "group", "plaingroup")
  markVar=""
  markVarVal=""
  #NOTE: pattern "plaingroup" is recommended only for large groups (more than 4 objects in a group)

  defaultPallet <- list(lightwhite=c(0,0,0,0), darkwhite=c(0,0,0,0.07),
                        lightgray =c(0,0,0,0.2), darkgray=c(0,0,0,0.27),
                        middlegray=c(0,0,0,0.18), red=c(0,0.9,0.3,0))
  
  if (nrow(dataframe)==0) stop("Argument dataframe has zero number of rows. Nothing to work with.")
  if (length(setdiff(c(rightVars, centerVars, fixedColVars), vars))>0) stop("Please make sure that all variables from rightVars, centerVars, fixedColVars are in vars")
  if (length(c(unique(rightVars), unique(centerVars), unique(fixedColVars)))>length(unique(c(rightVars, centerVars, fixedColVars)))) stop("Please make sure that variables sets rightVars, centerVars, fixedColVars do not intersect")
  if (length(vars)!=length(colNames)) stop("Make sure the length of 'vars' is the same as the length of 'colNames'")
  
  processBeginCommand <- function(beginCom=c(), outFile,
                                  caption = "", fontSize="small", colNames = c(),
                                  dataframe, zebraPattern,
                                  fixedColVars=c(), fixedColWdths=c(),
                                  markVar="", markVarVal="", placementStr, marker=NULL){
    
    #internal constants and functions definitions
    commandBegin <- function(command, outFile){
      if (command %in% fontSizes){
        cat("\n{\\",command,"\n",sep="", file=outFile)
      }else{
        cat("\n\\begin{",command,"}",sep="", file=outFile)
      }
    }
    commandEnd <- function(command, outFile){
      if (command %in% fontSizes){
        cat("}\n", file=outFile)
      }else{
        if (command == "table" | command == "tabular" | command == "longtable"){
          if (!is.null(marker)){cat("\\label{",marker,"}\n", sep="", file=outFile)}
        }
        cat("\\end{",command,"}\n",sep="", file=outFile)
      }
    }
    latexCaption <- function(captionFill, longtable, outFile){
      cap <- paste("\n\\caption{",captionFill,"}",sep="")
      if (longtable){
        cat(cap,"\\\\\n", sep="", file=outFile)
      }else{
        cat(cap,"\n", sep="", file=outFile)
      }
    }
    processColsFormat <- function(data, fixedColVars=c(), fixedColWdths=c(), outFile, placementStr){
      colFormat <- c()
      for (n in names(data)){
        if (n %in% fixedColVars){
          colFormat <- c(colFormat,paste("p{",fixedColWdths[fixedColVars==n],"pt}", sep=""))
        }else{
          if (n %in% rightVars){
            colFormat <- c(colFormat,"r")
          }else{
            if (n %in% centerVars){
              colFormat <- c(colFormat,"c")
            }else{
              colFormat <- c(colFormat,"l")
            }
          }
        }
      }
      cat(placementStr," {",paste(colFormat, collapse=""),"}", sep="", file=outFile)
    }
    hline <- function(outFile, number=1){
      cat(paste(rep("\\hline",number),collapse=""),"\n", file=outFile)
    }

    processColsHead <- function(colNames, longtable, caption="", outFile){
      processColNames <- function(colNames, style, outFile){
        hline(outFile,2)
        cat(style, colNames[1], file=outFile)
        for (i in 2:length(colNames)){
          cat("&", style, colNames[i], file=outFile)
        }
        cat("\\\\\n", file=outFile)
        hline(outFile)
      }
      headStyle <- "\\bfseries"
      otherStyle <- "\\bfseries\\em"
      colNum <- length(colNames)
      processColNames(colNames, headStyle, outFile)
      hline(outFile)
      if (longtable){
        cat("\\endfirsthead\n", file=outFile)
        cat("\\caption[]{",caption,"{",otherStyle," (continued)}} \\\\\n", file=outFile)
        processColNames(colNames, headStyle, outFile)
        cat("\\endhead\n", file=outFile)
        hline(outFile)
        cat("\\multicolumn{",colNum,"}{r}{",otherStyle," Continued on next page}\\\\\n",
            sep="", file=outFile)
        cat("\\endfoot\n", file=outFile)
        hline(outFile)
        cat("\\multicolumn{",colNum,"}{r}{",otherStyle," End}\\\\\n",
            sep="", file=outFile)
        cat("\\endlastfoot\n", file=outFile)
      }
    }
    processRows <- function(data, zebraPattern, outFile, markVar, markVarVal){
      #internal constants and functions definitions

      processRow <- function(row, color, outFile, markVarIndex){
        rowStr <- ""
        if (!is.na(color)){
          rowStr <- paste(rowStr,"\\rowcolor{",color,"}\n",sep="")
        }
        rowStr <- paste(rowStr,row[[1]],sep="")
        for (i in c(2:length(row))){
          if (i==markVarIndex){
            rowStr <- paste(rowStr," &", "\\color{red}{\\bfseries\\em ", row[[i]],"}",sep="")
          }else{
            rowStr <- paste(rowStr," &", row[[i]],sep="")
          }
        }
        rowStr <- paste(rowStr,"\\\\\n", sep="")
        cat(rowStr, file=outFile)
      }
      
    #beginning of the function processRows
      markVarIndex <- match(markVar, names(data))
      if (markVar!="" & is.na(markVarIndex)) stop("Please make sure that markVar is in the dataframe names\n")
      if (length(data[[1]]) != 0){
        for (i in c(1:length(data[[1]]))){
          if (!is.na(markVarIndex)){
            if (data[i,markVarIndex]==markVarVal){
              processRow(data[i,], zebraPattern[i], outFile, markVarIndex)
            }else{
              processRow(data[i,], zebraPattern[i], outFile, -1)
            }
          }else{
            processRow(data[i,], zebraPattern[i], outFile, -1)
          }
        }
      }
      hline(outFile)
    }

        
  #beginning of the function processBeginCommand
    if (length(beginCom[!is.na(beginCom)])>0){
      commandBegin(beginCom[1],outFile)
      if (beginCom[1] == "table"){
        latexCaption(caption, beginCom[1] == "longtable",outFile)
      }
      if (beginCom[1] == "tabular"){
        processColsFormat(data=dataframe, fixedColVars, fixedColWdths, outFile, placementStr)
        processColsHead(colNames = colNames, beginCom[1] == "longtable", outFile=outFile)
        processRows(data=dataframe, zebraPattern=zebraPattern, outFile,
                    markVar, markVarVal)
      }
      if (beginCom[1] == "longtable"){
        processColsFormat(data=dataframe, fixedColVars, fixedColWdths, outFile, placementStr)
        latexCaption(caption, beginCom[1] == "longtable",outFile)
        processColsHead(colNames = colNames, beginCom[1] == "longtable", caption, outFile)
        processRows(data=dataframe, zebraPattern=zebraPattern, outFile,
                    markVar, markVarVal)
      }
      processBeginCommand(beginCom[2:(length(beginCom)+1)], outFile,
                          caption, fontSize, colNames,
                          dataframe, zebraPattern,
                          fixedColVars, fixedColWdths,
                          markVar, markVarVal, placementStr, marker)
#      if (beginCom[1] == "table" | beginCom[1] == "tabular" | beginCom[1] == "longtable"){
#        if (!is.null(marker)){cat("\\label{",marker,"}\n", sep="", file=outFile)}
#      }
      commandEnd(beginCom[1], outFile)
    }
  }

  makePattern <- function(zebraPattern, by){
    #internal constants and functions definitions

    plain <- function(col1, col2, len){
      pattern <- rep(c(col1,col2), len)
      pattern <- pattern[1:len]  
      pattern
    }
    group <- function(by){
      col <- 1
      pattern <- rep(col, length(by))
      if (length(pattern)>1){
        current <- by[1]
        for (i in 2:length(by)){
          nextg <- by[i]
          if (current==nextg) {
            pattern[i] <- col
          }else{
            col <- col*(-1)
            pattern[i] <- col
          }
          current <- nextg
        }
      }
      pattern
    }
    
  #beginning of the function makePattern
    if (zebraPattern != "none"){
      if (!(zebraPattern %in% zebraPatterns)) stop("Error: Zebra Pattern can be one of the following: none, group, plain, plaingroup\n")
      lightwhite = "lightwhite"
      darkwhite = "darkwhite"
      lightgray = "lightgray"
      middlegray = "middlegray"
      darkgray = "darkgray"
      if (zebraPattern=="plain"){
        pattern <- plain(middlegray, lightwhite, length(by))
      }
      if (zebraPattern=="group"){
        pattern <- group(by)
        pattern <- ifelse(pattern>0, middlegray, lightwhite)
      }
      if (zebraPattern=="plaingroup"){
        pattern <- group(by)
        light <- plain(darkwhite, lightwhite, length(by))
        dark  <- plain(lightgray, darkgray, length(by))
        pattern <- ifelse(pattern>0,dark,light)
      }
      pattern
    }else{
      rep(NA, length(by))
    }
  }

  defineColors <- function(pallet, outFile){
    for (n in names(pallet)){
      cat("\\definecolor{",n,"}{cmyk}{",paste(pallet[[n]], collapse=","),"}\n",sep="",file=outFile)
    }
  }
    
#beginning of the function listTable

###-----------------------------open connection
  if (append){
    openMode <- "at"
  }else{
    openMode <- "wt"
  }
  if (typeof(file) == "character" && length(file) == 1){
    if (file=="") {
      outFile = stdout()
    }else{
      outFile = file(file, openMode)
      on.exit(close(outFile))
    }
  } else if (inherits(file,"connection") && (summary(file)$opened == "opened")){
    outFile = file
  }

###-----------------------------define pallet
  if (is.null(pallet)){
    pallet = defaultPallet
  }else{
    if(class(pallet)!="list") stop("Supply the pallet as a list")
    if(!all(sapply(pallet, length)==4)) stop("Make each element of pallet a vector of length 4")
    if(!all(sapply(pallet, max)<=1 & sapply(pallet, min)>=0)) stop("Make each element of pallet a vector with values from 0 to 1")
  }

###-----------------------------order dataframe
  if (orderGroups){
    dataframe <- dataframe[order(dataframe[[by]]),]
  }
  patternVar = dataframe[[by]]
  dataframe <- dataframe[,vars]
  if (zebraPattern %in% c("group","plaingroup")){
    ordered <- patternVar[order(patternVar)]
    if (!all(patternVar==ordered) && !all(patternVar==rev(ordered))){
      cat("\nWARNING: It is recommended to order the data by", by, "\n",
          "         when argument 'zebraPattern' is set to 'group' or 'plaingroup'.\n",
          "         It can be done by setting argument 'orderGroups' to TRUE.\n")
    }
  }

###-----------------------------create pattern
  cat(unique((pattern)), "\n")
  if (!is.null(pattern)){
    cat(unique(names(pattern)), "\n")
    if (!all(unique(pattern) %in% names(pallet))){
      stop("Please make sure that the color names of the pattern are in the pallet color names")
    }
    if (length(pattern)!=nrow(dataframe)) stop("The length of the pattern should be equal to the number of rows of the dataframe")
  }else{
    if (!all(names(defaultPallet) %in% names(pallet))) stop(paste("Please supply a pallet with names:", names(defaulePallet)))
    pattern <- makePattern(zebraPattern, patternVar)
  }

###-----------------------------take care of placement
  placementStr = ""
  if (!is.null(placement))
    if (!(placement %in% c("h","t","b","p"))){
      stop("Argument 'placement' should be one of the following 'h','t','b','p'\n")
    }else{
      placementStr = paste("[",placement,"]", sep="")
    }

###-----------------------------convert characters to latex text characters
  for (n in names(dataframe)){
    dataframe[[n]] <- as.character(dataframe[[n]])
    if (toLatexChar){
      dataframe[[n]] <- sapply(X = dataframe[[n]], FUN=latexTextMode)
      #dataframe[[n]] <- sapply(X = dataframe[[n]], FUN=latexTextModeNoRegExpr)
    }
  }
  
###-----------------------------process data
  beginCommands <- c()
  if (landscape) beginCommands <- c(beginCommands,"landscape")
  beginCommands <- c(beginCommands, fontSize)
  if (!landscape) beginCommands <- c(beginCommands,"center")
  
  if (longtable) beginCommands <- c(beginCommands,"longtable")
  else beginCommands <- c(beginCommands,c("table","tabular"))
  defineColors(pallet, outFile)
  processBeginCommand(beginCommands, outFile,
                      caption, fontSize, colNames,
                      dataframe, pattern,
                      fixedColVars, fixedColWdths,
                      markVar, markVarVal, placementStr, marker)
}#end of the function listTable





if (FALSE){

###---------------------------------------------
### MANAGING LATEX REFERENCES
###---------------------------------------------

getReferenceObject <- function(optionName){
  refD <- options(optionName)[[1]]
  if (is.null(refD)){
    refD <- data.frame(marker=c(), keyword=c(), label=c())
  }
  refD
}

putReferenceObject <- function(refD){
  options(rreport.reference.list=refD)
}

#print.latexReference <- function(refD){
#  if (is.null(refD)){
#    cat("The list of markers has not been created. Use function getReferenceObject() to create it\n")
#  }else{
#    print(refD)
#  }
#}

updateMarkers <- function(newMarker, keyword="", label="", optionName){
  ### puts a new marker into the dataframe of the existing ones
  ### checks if it is different from the existing ones
  ### returns updated latexReference
  refD = getReferenceObject(optionName)
  if (newMarker %in% refD$marker){
    stop(paste("Duplicated marker", newMarker))
  }
  newM <-data.frame(marker=newMarker, keyword=keyword, label=label)
  newM$marker <- as.character(newM$marker)
  newM$keyword <- as.character(newM$keyword)
  newM$label <- as.character(newM$label)
  refD = rbind(refD, newM)
  for (n in names(refD)) refD[[n]] <- as.character(refD[[n]])
  putReferenceObject(refD)
}

generateRef <- function(){
  generate <- function(){paste("marker",abs(round(rnorm(1)*(10^8))), sep="")}
  existingMarkers <- getRefsByKey()
  newMarker <- generate()
  while (newMarker %in% existingMarkers){
    newMarker <- generate()
  }
  newMarker
}

getReference <- function(keyword="", label=""){
  optionName = "manage.Latex.Reference"
  newMarker <- generateRef()
  updateMarkers(newMarker = newMarker, keyword=keyword, label=label, optionName)
  newMarker
}

getRefsByKey <- function(keyword=NULL){
  ### returns all markers with a given keyword 
  ### if keyword==NULL returns all markers
  refD = getReferenceObject()
  if (!is.null(keyword)){
    refD$marker[refD$keyword==keyword]
  }else{
    refD$marker
  }
}

getLabelsByKey <- function(keyword=NULL){
  ### returns all labels with a given keyword 
  ### if keyword==NULL returns all labels
  refD = getReferenceObject()
  if (!is.null(keyword)){
    refD$label[refD$keyword==keyword]
  }else{
    refD$label
  }
}

getReferenceString <- function(keyword){
  ### returns a vector of strings "see section \\ref{m1} (page\\pageref{m1})"
  ### for all markers with a given keyword 
  markers <- getRefsByKey(keyword)
  labels <- getLabelsByKey(keyword)
  keys <- paste(labels," in section ", "\\ref{", markers, "}", " (page ", "\\pageref{",markers,"}",")", sep="")
  paste("See", paste(keys, collapse=", "))
}




optionName = "manage.Latex.Reference"

getReferenceObject <- function(optionName){
  refD <- options(optionName)[[1]]
  if (is.null(refD)){
    refD <- data.frame(marker=c(), keyword=c(), label=c())
  }
  refD
}


getReference1 <- function(keyword="", label=""){
  updateMarkers <- function(newMarker, keyword="", label="", optionName){
    ### puts a new marker into the dataframe of the existing ones
    ### checks if it is different from the existing ones
    ### returns updated latexReference
    putReferenceObject <- function(refD){
      options(rreport.reference.list=refD)
    }
    refD = getReferenceObject(optionName)
    if (newMarker %in% refD$marker){
      stop(paste("Duplicated marker", newMarker))
    }
    newM <-data.frame(marker=newMarker, keyword=keyword, label=label)
    newM$marker <- as.character(newM$marker)
    newM$keyword <- as.character(newM$keyword)
    newM$label <- as.character(newM$label)
    refD = rbind(refD, newM)
    for (n in names(refD)) refD[[n]] <- as.character(refD[[n]])
    putReferenceObject(refD)
  }
  generateRef <- function(){
    generate <- function(){paste("marker",abs(round(rnorm(1)*(10^8))), sep="")}
    existingMarkers <- getRefsByKey()
    newMarker <- generate()
    while (newMarker %in% existingMarkers){
      newMarker <- generate()
    }
    newMarker
  }
  newMarker <- generateRef()
  updateMarkers(newMarker = newMarker, keyword=keyword, label=label, optionName)
  newMarker
}

getRefsByKey1 <- function(keyword=NULL){
  ### returns all markers with a given keyword 
  ### if keyword==NULL returns all markers
  refD = getReferenceObject()
  if (!is.null(keyword)){
    refD$marker[refD$keyword==keyword]
  }else{
    refD$marker
  }
}

getLabelsByKey1 <- function(keyword=NULL){
  ### returns all labels with a given keyword
  ### if keyword==NULL returns all labels
  refD = getReferenceObject()
  if (!is.null(keyword)){
    refD$label[refD$keyword==keyword]
  }else{
    refD$label
  }
}

}

pdftex = function(pathToPdf="pdf", pathToTex="gentex", fileName, caption, longCaption=caption, label=NULL, width=6, height=6, fontsize="normalsize"){
### opens pdf connection and writes a pdf wrapper in a tex file
### fileName - name of the file without the pdf extension, and NOT a PATH
### caption - caption of the pdf figure
### label - label that you can reference later in latex
  if (!any(fontsize %in% c("tiny","scriptsize","footnotesize","small", "normalsize","large","Large","LARGE","huge","Huge")))
    stop("Please pick argument fontsize out of: tiny,scriptsize,footnotesize,small, normalsize,large,Large,LARGE,huge,Huge" )
  append=FALSE
  cleanFileName = gsub("\\.pdf$", "", fileName)  ### strips off pdf extension
  pdfFileName = paste(cleanFileName, ".pdf", sep="")
  texFileName = paste(cleanFileName, ".tex", sep="")
  fullPdfFileName = paste(pathToPdf, pdfFileName, sep="/")
  fullTexFileName = paste(pathToTex, texFileName, sep="/")
  pdf(fullPdfFileName, width=width, height=height)
  cat("\\begin{figure}[hbp!]\\leavevmode\\centerline{\\includegraphics{", pdfFileName, "}}\n", sep="", file=fullTexFileName, append=append)
  append=TRUE
  cat("\\caption[",caption,"]{","\\",fontsize,"{",longCaption,"}","}\n", sep="", file=fullTexFileName, append=append)
  if (!is.null(label)){
    cat("\\label{",label,"}\n", sep="", file=fullTexFileName, append=append)
  }
  cat("\\end{figure}\n", sep="", file=fullTexFileName, append=append)
}

sinktex = function(fileName=NULL, verbatim = TRUE){
  if (is.null(fileName)){
    if(verbatim){
      cat("\\end{verbatim}\n")
    }
    sink()
  }else{
    sink(fileName)
    if(verbatim){
      cat("\\begin{verbatim}\n")
    }
  }
}
