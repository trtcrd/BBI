
## behaviour: sum: colSums, mean: colMeans
mergeOTUTable = function(otutable, metadata, by="Group", behaviour = "sum"){
  if (dim(otutable)[1] != dim(metadata)[1]){
    print("Error: number of rows in OTU table is not equal to metadata!")
    return()
  }
  if (!identical(row.names(otutable),row.names(metadata))){
    print("Warning: row names differ between OTU table and metadata")
  }
  if (!(behaviour %in% c("sum","mean"))){
    print("behaviour must be either sum or mean!")
    return
  }
  groups = unique(metadata[[by]])
  ng = length(groups)
  mergedOTUs = matrix(nrow = ng, ncol = dim(otutable)[2])
  for (i in c(1:ng)){
    toMerge = otutable[metadata[[by]] == groups[i],]
    if (behaviour == "sum") mergedOTUs[i,] = colSums(toMerge)
    else if (behaviour == "mean") mergedOTUs[i,] = colMeans(toMerge)
  }
  mergedOTUs = data.frame(mergedOTUs)
  row.names(mergedOTUs) = groups
  names(mergedOTUs) = names(otutable)
  
  return(mergedOTUs)
}



mergeMetadata = function(metadata, by="Group"){
  groups = unique(metadata[[by]])
  ng = length(groups)
  np = dim(metadata)[2]
  merged = matrix(nrow = ng, ncol = np)
  merged = as.data.frame(merged)
  for (i in c(1:ng)){
    toMerge = metadata[metadata[[by]] == groups[i],]
    if (dim(toMerge)[2]>0) {
      for (j in c(1:np)){
        if (is.numeric(toMerge[1,j])) {
          #print(paste("DEBUG: ",names(toMerge)[j], "is numeric"))
          merged[i,j] = mean(toMerge[,j],na.rm=T)
        } else {
          #print(paste("DEBUG: ",names(toMerge)[j], "is text"))
          merged[i,j] = as.character(toMerge[1,j])
        }
      }
    }
  }
  #merged = as.data.frame(merged)
  row.names(merged) = groups
  names(merged) = names(metadata)
  
  return(merged)
}

