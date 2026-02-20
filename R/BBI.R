#' BBI : Benthic Biotic Indices calculation function
#'
#' @description The \code{BBI} package is meant to calculate Benthic Biotic Indices from
#' composition data, obtained whether from morphotaxonomic inventories or
#' seBI_metab <- BBI(metab)quencing data. Based on reference ecological weights publicly available for
#' a set of commonly used marine biotic indices, such as AMBI (Borja et al., 2000)
#' NSI and ISI indices (Rygg 2013).
#'
#' New since v1.0:
#' Requires stringdist package
#' Does not take median of genus match but searches regressively if taxonomic path
#' is given until a hit is found. Genera for species that agree can be put in
#' index instead (as it should be done) 
#'
#' @param data A data frame containing samples as columns and taxa as rows by default, with
#'   taxonomic path (or best taxon) in the first column of \code{data}. Optionally the row.names 
#'   if the first column is numeric
#'
#' @return Function \code{BBI} returns a list of containing:
#'   \enumerate{
#'   \item \code{found} - the amount of taxa that matched an entry in the
#'   database and the amount that did not
#'   \item \code{BBI} - the BBI values per sample
#'   \item \code{table} - the subset of composition data that contains only taxa with
#'   at least a match in one of the BBI
#'   \item \code{taxa} - the list of taxa that matched an entry and the correspondant OTU, if from NGS data.
#'
#' @example BBI(my_table, log = FALSE)
#' @author Tristan Cordier
#' @author Anders Lanzen
#' Re-implemented by Anders 2025
#' @export
#' BBI()

require("stringdist")

# sourcing the others of functions
source("R/nEQR.R")
source("R/status.R")
source("R/mergeOTUTable.R")

# Match function
# returns index of best matching taxon in reference list or NA if not found
matchTaxon = function(taxText, bi_taxa, maxDist=1){
  
  # Text to substitute - trim to only taxon (after semicolon, comma or slash)
  # and remove any occurence of sp. or spp.
  tf = tolower(as.character(taxText))
  
  if(grepl("sp",tf)) {
    tf = gsub(" sp.$","",tf)
    tf = gsub(" sp$","",tf)
    tf = gsub(" spp.$","",tf)
    tf = gsub(" spp$","",tf)
  }
  
  tf = gsub("uncultured  ","",tf)
  tf = gsub("unknown  ","",tf)
  tf = gsub("cf.$","",tf)
  tf = gsub("indet.$","",tf)
  
  # Make list of taxonomic path if given, otherwise just the name
  if(grepl(";",tf)) {
    taxpath = unlist(strsplit(tf, split=";", fixed=TRUE))
  }
  else if(grepl("/",tf)) {
    taxpath = unlist(strsplit(tf, split="/", fixed=TRUE))
  }
  else if(grepl(",",tf)) {
    taxpath = unlist(strsplit(tf, split=",", fixed=TRUE))
  }
  else taxpath = tf
  
  
  # Search from lowest to highest rank until something is found, or not
  bestMatch = NA
  i=0
  while(i<length(taxpath) & is.na(bestMatch)){
    i = i+1
    searchTaxon = tail(taxpath,i)[1]
    if(nchar(searchTaxon)>1) bestMatch = amatch(searchTaxon, bi_taxa, maxDist=maxDist)
  }
  
  # Return the index of the best hit
  return(bestMatch)
}

# BBI main function
BBI <- function(data, log=FALSE, maxMatchErrors=1) {
  # if data is not a data.frame
  data <- as.data.frame(data)
  
  # If the first column is not numeric use it as names and remove
  if(!is.numeric(data[,1])) {
    queries=data[,1]
    data = data[,-1]
  }
  else queries = row.names(data)
  
  # Logging the search?
  if(log) log_file <- file(paste("Log_BBI_", format(Sys.time(), "%Y_%a_%b_%dth-%H.%M.%S"), ".txt", sep=""),open="a")
  
  # Keep only OTUs with non-zero abundances
  message("Removing OTUs with zero total abundance across samples:")
  if(log) cat("Removing OTUs with zero total abundance across samples:", 
              file=log_file, fill=T, append=T)
  otus_positive = data[rowSums(data)>0,]
  queries = queries[rowSums(data)>0]
  message(paste(dim(otus_positive)[1],"out of",dim(data)[1],"OTUs kept."))
  if(log) cat(paste(dim(otus_positive)[1],"out of",dim(data)[1],"OTUs kept."), 
                      file=log_file, fill=T, append=T)
  
  ## import the reference BI table
  message("Reading reference table")
  if(log) cat("Reading reference table", file=log_file, fill=T, append=T)
  
  #print(paste(system.file(package="BBI"), "/TABLE_REF.Rd", sep="")) #DEBUG
  eco_index <- read.table(paste(system.file(package="BBI"), "/TABLE_REF.Rd", sep=""), header=TRUE, sep="\t", dec=".")
  
  # Make data frame with all quieries with columns for best hit found and its index values
  query_results = data.frame(row.names=row.names(otus_positive), query=queries, found = FALSE, 
                             index=NA, bestHit=NA, AMBI=NA, ITI_GROUP=NA, ISI_value=NA,
                             NSI=NA, NSI.group=NA, Bentix=NA)
  
  # for each taxa
  for (i in 1:length(queries)){
    
    # Try to match to biotic indices table
    y = matchTaxon(queries[i], eco_index$query, maxDist = maxMatchErrors)
    
    if(!is.na(y)){
      query_results[i,"found"] = TRUE
      
      # Fill in query_results
      query_results[i,"index"] = y
      query_results[i,"bestHit"] = eco_index$TAXA[y]
      query_results[i,"AMBI"] = eco_index$AMBI[y]
      query_results[i,"ITI_GROUP"] = eco_index$Iti_group[y]
      query_results[i,"ISI_value"] = eco_index$ISI.2018[y]
      query_results[i,"NSI"] = eco_index$NSI.value[y]
      query_results[i,"NSI.group"] = eco_index$NSI.group[y]
      query_results[i,"Bentix"] = eco_index$bentix[y]
      
      message(paste("Matched",queries[i],"to",query_results$bestHit[i]))
      #message(query_results[i,]) #DEBUG
      if(log) cat(paste("Matched",queries[i],"to",query_results$bestHit[i]), file=log_file, fill=T, append=T)
    }
    else{
      message(paste("No hit found for ",queries[i]))
      if(log) cat(paste("No hit found for ",queries[i]), file=log_file, fill=T, append=T)
    }
  }
  # Summarise found / not found stats
  foundN = sum(query_results$found)
  unfound = length(queries)-foundN
  message(paste("=== In total found",foundN,"matches.===="))
  message(paste("=== No match for",unfound,"queries.===="))
          
  if(log) {
    cat(paste("=== In total found",foundN,"matches.===="),
        file=log_file, fill=T, append=T)
    cat(paste("=== No match for",unfound,"queries.===="),
        file=log_file, fill=T, append=T)
  }
  
  # Exit here if no taxa were identified
  if(foundN==0) return()

  # bind all taxa, eco-weights, and composition data
  #output <- cbind(tax_n, out, dat)
  found_otus = otus_positive[query_results$found,]
  found_results = query_results[query_results$found,]
  
  ## Create taxon table summing OTU abundances and corresponding results table
  found_taxa_table = mergeOTUTable(found_otus, by = "query", metadata = found_results)
  found_results_u = mergeMetadata(found_results, by = "query")
  totalTaxa = dim(found_taxa_table)[1]
  
  ## compute shannon (base 2) on taxa that got at least one value (not the best, but it is computed like that in Norway)
  tmp_sha = as.data.frame(t(found_taxa_table))
  output_shannon <- vegan::diversity(tmp_sha, index="shannon", base = 2)
 
  ############################################################
  ### now compute all indices.
  ############################################################

  indexNames = c("AMBI", "ISI", "NSI", "NQI1", "Shannon", "ITI", "Bentix")
  # indices is a matrix with these indices for each sample
  indices = matrix(nrow=length(indexNames),ncol=dim(found_taxa_table)[2],)
    # <- as.data.frame(array(NA, c(7,dim(data_)[2]-6)))
  # dimnames(indices)[[1]] <- c("AMBI", "ISI", "NSI", "NQI1", "Shannon", "ITI", "Bentix")
  # dimnames(indices)[[2]] <- dimnames(data_)[[2]][7:dim(data_)[2]]

  indices[5,] <- output_shannon
 
  # AMBI
  ambi_eg <- found_results_u$AMBI
  partial_ambi = matrix(nrow = dim(found_taxa_table)[1], ncol=dim(found_taxa_table)[2])
  
  # total abundance assigned
  ab_w_ambi = colSums(found_taxa_table[!is.na(ambi_eg),])
  
  # for each taxon in the sample with abundances
  for (j in 1:totalTaxa) {
    # if there is a value for AMBI
    if (!is.na(ambi_eg[j]))  {
      # for AMBI value 2 (weights starts at 2)
      if (ambi_eg[j] == 1)  partial_ambi[j,] <- 0
      else if (ambi_eg[j] == 2)  partial_ambi[j,] <- t(1.5 * found_taxa_table[j,])
      else if (ambi_eg[j] == 3)  partial_ambi[j,] <- t(3 * found_taxa_table[j,])
      else if (ambi_eg[j] == 4)  partial_ambi[j,] <- t(4.5 * found_taxa_table[j,])
      else if (ambi_eg[j] == 5)  partial_ambi[j,] <- t(6 * found_taxa_table[j,])
    }
  }
  ambi = colSums(partial_ambi, na.rm=T) / ab_w_ambi
  indices[1,] = ambi

  # NSI
  ab_w_nsi = colSums(found_taxa_table[!is.na(found_results_u$NSI),])
  partial_nsi = matrix(nrow = dim(found_taxa_table)[1], ncol=dim(found_taxa_table)[2])
  
  for (j in 1:totalTaxa) {
    # if there is a sensitivity value for NSI
    if (!is.na(found_results_u$NSI[j])) 
      partial_nsi[j,] <- t(found_results_u$NSI[j] * found_taxa_table[j,])
  }
  nsi = colSums(partial_nsi, na.rm=T) / ab_w_nsi
  indices[3,] = nsi
  
  # ISI 2018 (ignores counts and looks at p/a)
  ftt_pa = vegan::decostand(found_taxa_table, method="pa")
  isiN = colSums(ftt_pa[!is.na(found_results_u$ISI_value),])
  #print(isiN)
  partial_isi = matrix(nrow = dim(found_taxa_table)[1], ncol=dim(found_taxa_table)[2])
  
  for (j in 1:totalTaxa) {
    # if there is a sensitivity value for NSI
    if (!is.na(found_results_u$ISI_value[j])) 
      partial_isi[j,] <- t(found_results_u$ISI_value[j] * ftt_pa[j,])
  }
  isi = colSums(partial_isi, na.rm=T) / isiN
  indices[2,] = isi
  
  # ITI index
  iti_eg <- found_results_u$ITI_GROUP
  partial_iti = matrix(nrow = dim(found_taxa_table)[1], ncol=dim(found_taxa_table)[2])
  
  # total abundance assigned
  ab_w_iti = colSums(found_taxa_table[!is.na(iti_eg),])
  
  # for each taxon in the sample with abundances
  for (j in 1:totalTaxa) {
    # if there is a value for AMBI
    if (!is.na(iti_eg[j]))  {
      # for AMBI value 2 (weights starts at 2)
      if (iti_eg[j] == 1)  partial_iti[j,] <- 0
      else if (iti_eg[j] == 2)  partial_iti[j,] <- t(found_taxa_table[j,])
      else if (iti_eg[j] == 3)  partial_iti[j,] <- t(2 * found_taxa_table[j,])
      else if (iti_eg[j] == 4)  partial_iti[j,] <- t(3 * found_taxa_table[j,])
    }
  }
  iti = 100-(100/3*(colSums(partial_iti, na.rm=T) / ab_w_iti))
  indices[6,] = iti
  
  ## NQI1 (only works if we have OTU or individual counts, not relative abundaneces)
  if(sum(found_taxa_table<1 & found_taxa_table>0)==0){
    totalCounts = colSums(found_taxa_table)
    totalSp = colSums(vegan::decostand(found_taxa_table, method="pa"))
    nqi1 = 0.5*(1-ambi/7) +0.5*(log(totalSp)/log(log(totalCounts))/2.7)*(totalCounts/(totalCounts+5))
    indices[4,] = nqi1
  }
  
  
  ### BENTIX
  
  partial_b = matrix(nrow = dim(found_taxa_table)[1], ncol=dim(found_taxa_table)[2])
  
  # total abundance assigned
  ab_w_b = colSums(found_taxa_table[!is.na(found_results_u$Bentix),])
  
  # for each taxon in the sample with abundances
  for (j in 1:totalTaxa) {
    # if there is a value for AMBI
    if (!is.na(found_results_u$Bentix[j]))  {
      if (found_results_u$Bentix[j] == 1)  partial_b[j,] <- t(6 * found_taxa_table[j,])
      else if (found_results_u$Bentix[j] == 2 | found_results_u$Bentix[j] == 3)  partial_b[j,] <- t(2*found_taxa_table[j,])
    }
  }
  bentix = 100*colSums(partial_b, na.rm=T)/ab_w_b
  
  indices[7,] = bentix
 
  ## Convert indices to deataframe
  idf = as.data.frame(t(indices))
  names(idf) = indexNames
  row.names(idf) = names(found_taxa_table)
  #summary(idf)
  
  ## now we can return the dicrete assessment for each BBI
  
  idf$AMBI_group = NA
  idf$ISI_group = NA
  idf$NSI_group = NA
  idf$NQ1_group = NA
  idf$Shannon_group = NA
  
  for (i in 1:nrow(idf)){
    idf[i,"AMBI_group"] = e$status.ambi(idf$AMBI[i])
    idf[i, "ISI_group"] = e$status.isi(idf$ISI[i])
    idf[i,"NSI_group"] = e$status.nsi(idf$NSI[i])
    idf[i,"NQ1_group"] = e$status.nqi1(idf$NQI1[i])
    idf[i,"Shannon_group"] = e$status.shannon(idf$Shannon[i])
  }
  idf$AMBI_group=as.factor(idf$AMBI_group)
  idf$ISI_group=as.factor(idf$ISI_group)
  idf$NSI_group=as.factor(idf$NSI_group)
  idf$NQ1_group=as.factor(idf$NQ1_group)
  idf$Shannon_group=as.factor(idf$Shannon_group)
  
  # preparing the output
  output <- list("found" = c("Found match:", foundN, " Not found:", unfound ),
                 "queries" = queries,
                 "BBI" = idf,
                 "BBIclass" = idf[,c(8:12)],
                 "table" = query_results,
                 "taxa" = queries)
  return(output)
}
