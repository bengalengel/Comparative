ExpandPhos <- function(phospho){
  ##now to ensure each multiplicity has its own row and variables are condensed. I am converting from wide to long format to ensure that each
  # observation is uniquely represented in the data. Then I cast back and remove empty rows
  # This function accepts a dataframe 'phospho' loaded from MQ output
  require(plyr)
  require(reshape2)
  
  expression <- phospho[,grep("Intensity.[0-9]*", colnames(phospho))]
  
  melted <- melt(phospho, measure.vars = names(expression))
  
  # Here I split (that is add a variable identifier) the melted 'variable' column so that the sample, replicate, 
  # and multiplicity are now explicit
  melted <- cbind(melted, colsplit(melted$variable, "___", c("sample", "multiplicity")))
  melted$sample <- gsub(melted$sample, pattern = "Intensity.", replacement = "") ##remove redundant information next 3 lines
  
  #cast s.t. each unique sample has its own column
  casted <- dcast(melted, ... ~ sample, value.var = "value")
  
  #make the multiplicity explicity table
  
  ##gives index of experiment and replicate
  data <- grep("[0-9]", colnames(casted))
  
  ##gives string of experiment and replicate
  data2 <- colnames(casted)[data]
  
  #produces a new string with proper alpha leading R variable names (otherwise below won't work)
  newnames <- paste0("Int",data2)
  
  # rename the experiment variables within the dataframe
  colnames(casted)[data] <- newnames
  
  ## columnwise application of mean to condense the dataframe. PRODUCES CORRECT NUMBERS!
  out <- ddply(casted, .(id, multiplicity), colwise(mean,newnames,na.rm=T))
  
  #merge with identifying information by id to produce the multiplicity expanded table (each obs has a row)
  other_data <- phospho[, -which(names(phospho) %in% names(expression))]         
  
  multExpanded <- merge(other_data, out, by="id")
  
  ## remove rows with only NAs in expression columns
  
  expCol <- grep("Int[0-9]", colnames(multExpanded))

  multExpanded <- multExpanded[rowSums(multExpanded[,expCol])>0,]##removes rows containing all zeros
  
  return(multExpanded)
  }


