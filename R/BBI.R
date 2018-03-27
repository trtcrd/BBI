################################################################
# BBI : Benthic Biotic Indices calculation function
# author : tristan cordier
# contact : tristan.cordier [a] gmail.com
################################################################

# BBI function
BBI <- function(data) {
  # package requirement
  require(vegan)
  # if data is not a data.frame
  data <- as.data.frame(data)
  ## import the reference BI table
  eco_index <- read.table("TABLE_REF/TABLE_REF_BBI.txt", header=TRUE, sep="\t", dec=",")
  # table morpho
  table_test <- data
  tax_n <- table_test[,1]
  tax_n <- cbind(as.character(tax_n), as.character(tax_n))

  otu_id <- rownames(table_test)
  
  ## counting
  cpt_found <- 0 
  cpt_not <- 0
  found_index <- c()
  
  # need to convert each comumns into numeric (WEIRD behaviour : check http://stackoverflow.com/questions/3418128/how-to-convert-a-factor-to-an-integer-numeric-without-a-loss-of-information)
  if (is.factor(table_test[,2]) == T) for (i in 2:dim(table_test)[2]) table_test[,i] <- as.numeric(levels(table_test[,i])[table_test[,i]])
  
  # # if duplicate row names, the aggregate by taxa name
  # if (dim(table_test)[1] != length(unique(tax_n[,1]))) 
  # {
  #   dat_dna <- aggregate(table_test[,2:dim(table_test)[2]], by = list(tax_n[,1]), FUN = "sum")
  #   tax_n <- dat_dna[,1]
  #   tax_n <- cbind(as.character(tax_n), as.character(tax_n))
  #   otu_list <- grep("OTU", tax_n[,1], ignore.case = F)
  #   tax_n <- tax_n[-otu_list,]
  #   dat_dna <- dat_dna[-otu_list,2:dim(dat_dna)[2]]
  #   
  # } else {
  # in case of correlation screening, removing OTU without correlation to speed up process
  otu_list <- grep("OTU", tax_n[,1], ignore.case = F)
  if (length(otu_list) > 0) tax_n <- tax_n[-otu_list,]
  if (length(otu_list) > 0) dat_dna <-  table_test[-otu_list,2:dim(table_test)[2]] else dat_dna <-  table_test[,2:dim(table_test)[2]]
  # }
  
  # creating the OTU table with indices values binded
  out <- as.data.frame(array(NA, c(dim(tax_n)[1],6)))
  dimnames(out)[[2]] <- c("query","AMBI", "ITI_GROUP", "ISI_value", "NSI", "Bentix")
  
  for (i in 1:dim(tax_n)[1])
  {
    # keep the last value in taxonomy assignment
    sp <- tail(unlist(strsplit(as.character(tax_n[i,2]), split="|", fixed=TRUE)),1)
    # if we got a species or sp. as assignement, replace "+" by " "
    sp <- gsub("+", " ", sp, fixed=TRUE)
    sp <- gsub("_", " ", sp, fixed=TRUE)
    sp <- gsub(" Cmplx.", "", sp, fixed=TRUE)
    
    #sp <- gsub("Capitella sp. LFR-2016", "Capitella capitata", sp, fixed=TRUE)
    
    sp <- gsub(" environmental sample", "", sp, fixed=TRUE)
    sp <- gsub(" indet.", "", sp, fixed=TRUE)
    sp <- gsub(" indet", "", sp, fixed=TRUE)
    ## check if "sp" instead of "sp."
    if (length(strsplit(sp, split=" sp", fixed=TRUE)[[1]]) == 1) sp <- strsplit(sp, split=" sp", fixed=TRUE)[[1]][1]
    # check if something after sp.
    if (length(strsplit(sp, split="sp.", fixed=TRUE)[[1]]) > 1)
    {
      sp <- strsplit(sp, split=" sp.", fixed=TRUE)[[1]][1]
    } else if (length(strsplit(sp, split=" cf.", fixed=TRUE)[[1]]) > 1)
    {
      sp <- strsplit(sp, split=" cf.", fixed=TRUE)[[1]][1]
    } else if (length(strsplit(sp, split="(", fixed=TRUE)[[1]]) > 1)
    {
      sp <- strsplit(sp, split=" (", fixed=TRUE)[[1]][1]
    }
    # check if there a value in reference eco values
    y <- grep(sp, eco_index[,"species"], ignore.case = TRUE)
    print(paste("Processing : ", sp, " - ", length(y), " match in database so far", sep=""))
    # if the "exact" species is not matching and there a value for genus level
    if (length(y) == 0)
    {
      sp <- strsplit(sp, split=" ", fixed=TRUE)[[1]][1]
      print(paste("   No match, trying exact match of the genus", sp))
      #sp <- paste(sp, "sp.", sep=" ")
      # and then exact match
      y <- grep(paste("^", sp, " sp.", "$", sep=""), eco_index[,"species"], ignore.case = TRUE)
      # if still not, there is no "genus sp." in the table, grab the species of the genus and median value
      if (length(y) == 0)
      {
        print(paste("   No exact match, trying to match other species of the genus", sp))
        y <- grep(paste(sp, " ",sep=""), eco_index[,"species"], ignore.case = F)
        if (length(y) >= 1)
        {
          print(paste("   Found ", length(y), " match for ", sp, ". Taking the median values of ", length(y), " multiple values", sep = ""))
          ambi <- eco_index[y,"AMBI"]
          iti  <- eco_index[y,"Iti_group"]
          isi  <- eco_index[y,"ISI.2012"]
          nsi  <- eco_index[y,"NSI.value"]
          ben  <- eco_index[y,"bentix"]
          out[i,"AMBI"] <- ceiling(median(as.numeric(as.vector(ambi))))  # round to the ceiling value if median is not integer
          out[i,"ITI_GROUP"] <- median(as.numeric(as.vector(iti)))
          out[i,"ISI_value"] <- median(as.numeric(as.vector(isi)))
          out[i,"NSI"] <- median(as.numeric(as.vector(nsi)))
          out[i,"Bentix"] <- ceiling(median(as.numeric(as.vector(ben))))
          print(paste("   Done - ", sp, " AMBI value ", out[i,"AMBI"], sep=""))
          cpt_found <- cpt_found + 1
          found_index <- c(found_index, 1)
        } else {
          print("   Not found.")
          cpt_not <- cpt_not + 1
          found_index <- c(found_index, 0)
        }
      } else if (length(y) >= 1)
      {
        print(paste("   Found ", length(y), " match for ", sp, ". Taking the median values of ", length(y), " multiple values", sep = ""))
        ambi <- eco_index[y,"AMBI"]
        iti  <- eco_index[y,"Iti_group"]
        isi  <- eco_index[y,"ISI.2012"]
        nsi  <- eco_index[y,"NSI.value"]
        ben  <- eco_index[y,"bentix"]
        out[i,"AMBI"] <- ceiling(median(as.numeric(as.vector(ambi))))  # round to the ceiling value if median is not integer
        out[i,"ITI_GROUP"] <- median(as.numeric(as.vector(iti)))
        out[i,"ISI_value"] <- median(as.numeric(as.vector(isi)))
        out[i,"NSI"] <- median(as.numeric(as.vector(nsi)))
        out[i,"Bentix"] <- ceiling(median(as.numeric(as.vector(ben))))
        print(paste("   Done - ", sp, " AMBI value ", out[i,"AMBI"], sep=""))
        cpt_found <- cpt_found + 1
        found_index <- c(found_index, 1)
      }
    }
    # if there is at least a match
    else if (length(y) > 0)
    {
      #print(sp)
      tmp <- eco_index[y,]
      if (length(y) > 1) 
      {
        print(paste("   ", sp, " matched ", length(y), " entries in database. Trying exact match", sep=""))
        # if more than one match, try exact match
        y <- grep(paste("^", sp, "$", sep=""), eco_index[,"species"], ignore.case = TRUE)
        # if no perfect match, try adding a sp.
        if (length(y) == 0) 
        {
          sp <- paste(sp, "sp.", sep=" ")
          print(paste("   No exact match. Trying exact match with ", sp, sep=""))
          y <- grep(paste("^", sp, "$", sep=""), eco_index[,"species"], ignore.case = TRUE)
        }
        print(paste("   Found ", length(y), " match for ", sp, sep=""))
        if (length(y) == 1)
        {
          ambi <- eco_index[y,"AMBI"]
          iti  <- eco_index[y,"Iti_group"]
          isi  <- eco_index[y,"ISI.2012"]
          nsi  <- eco_index[y,"NSI.value"]
          ben  <- eco_index[y,"bentix"]
          out[i,"AMBI"] <- as.numeric(as.vector(ambi))
          out[i,"ITI_GROUP"] <- as.numeric(as.vector(iti))
          out[i,"ISI_value"] <- as.numeric(as.vector(isi))
          out[i,"NSI"] <- as.numeric(as.vector(nsi))
          out[i,"Bentix"] <- as.numeric(as.vector(ben))
          print(paste("   Done - ", sp, " AMBI value ", eco_index[y,"AMBI"], sep=""))
          cpt_found <- cpt_found + 1
          found_index <- c(found_index, 1)
        } else {
          print(paste("   Taking median of the ", length(y), " match for ", sp, sep=""))
          ambi <- eco_index[y,"AMBI"]
          iti  <- eco_index[y,"Iti_group"]
          isi  <- eco_index[y,"ISI.2012"]
          nsi  <- eco_index[y,"NSI.value"]
          ben  <- eco_index[y,"bentix"]
          out[i,"AMBI"] <- ceiling(median(as.numeric(as.vector(ambi)))) # ceiling value if not integer
          out[i,"ITI_GROUP"] <- median(as.numeric(as.vector(iti)))
          out[i,"ISI_value"] <- median(as.numeric(as.vector(isi)))
          out[i,"NSI"] <- median(as.numeric(as.vector(nsi)))
          out[i,"Bentix"] <- ceiling(median(as.numeric(as.vector(ben))))
          print(paste("   Done - ", sp, " AMBI value ", out[i,"AMBI"], sep=""))
          cpt_found <- cpt_found + 1
          found_index <- c(found_index, 1)
        }
        
        # if it is a genus only (deeper assignments might or not have a value)
        if (length(unlist(strsplit(sp, split=" "))) == 1 & is.na(out[i,"AMBI"]) == T)
        {
          print(paste("   ", sp, " is is a genus query. Making the query : ", sp, " sp.", sep=""))
          # make the query "query sp."
          sp <- paste(sp,"sp.", sep=" ")
          y  <- grep(sp, eco_index[,"species"], fixed=TRUE)
          # if there was not only species of the genus (no value for the genus itself)
          if (length(y) > 0)
          {
            print(paste("   Found - ", length(y), " for ", sp, sep=""))
            tmp <- eco_index[y,]
            ambi <- eco_index[y,"AMBI"]
            iti  <- eco_index[y,"Iti_group"]
            isi  <- eco_index[y,"ISI.2012"]
            nsi  <- eco_index[y,"NSI.value"]
            ben  <- eco_index[y,"bentix"]
            out[i,"AMBI"] <- as.numeric(as.vector(ambi))
            out[i,"ITI_GROUP"] <- as.numeric(as.vector(iti))
            out[i,"ISI_value"] <- as.numeric(as.vector(isi))
            out[i,"NSI"] <- as.numeric(as.vector(nsi))
            out[i,"Bentix"] <- as.numeric(as.vector(ben))
            print(paste("   Done - ", sp, sep=""))
            cpt_found <- cpt_found + 1
            found_index <- c(found_index, 1)
          }
        }
      } else {
        print(paste("   Found ", length(y), " matche for ", sp, sep=""))
        ambi <- eco_index[y,"AMBI"]
        iti  <- eco_index[y,"Iti_group"]
        isi  <- eco_index[y,"ISI.2012"]
        nsi  <- eco_index[y,"NSI.value"]
        ben  <- eco_index[y,"bentix"]
        out[i,"AMBI"] <- as.numeric(as.vector(ambi))
        out[i,"ITI_GROUP"] <- as.numeric(as.vector(iti))
        out[i,"ISI_value"] <- as.numeric(as.vector(isi))
        out[i,"NSI"] <- as.numeric(as.vector(nsi))
        out[i,"Bentix"] <- as.numeric(as.vector(ben))
        print(paste("   Done - ", sp, " AMBI: ", eco_index[y,"AMBI"], sep=""))
        cpt_found <- cpt_found + 1
        found_index <- c(found_index, 1)
        # print(paste("   AMBI: ", eco_index[y,"AMBI"], " ITI: ", eco_index[y,"Iti_group"],
      }
    }
    out[i,"query"] <- sp
  }
  
  print(paste("==== Found match :", cpt_found, " Not found :", cpt_not, "===="))
  
  output <- cbind(tax_n, out, dat_dna)
  
  otu_id_list <- cbind(tax_n, otu_id)
  
  ## compute shannon on taxa that got at least one value
  tmp_sha <- t(dat_dna[,1:dim(dat_dna)[2]])
  tmp_sha[is.na(tmp_sha)] <- 0
  output_shannon <- diversity(tmp_sha, index="shannon", base = 2)
  
  # to make the rownames having unique names when same species is identified... 
  # because dimnames(output)[[1]] <- tax_n[,1] does not want it
  dimnames(output)[[1]] <- make.unique(tax_n[,1])
  
  output <- output[,4:dim(output)[2]]
  
  print(output[,1:5])
  
  output <- subset(output, found_index==1)
  otu_id_list <- subset(otu_id_list, found_index==1)
  
  ### 
  
  ### now compute all indices.
  
  data_ <- output
  
  indices <- as.data.frame(array(NA, c(7,dim(data_)[2]-5)))
  dimnames(indices)[[1]] <- c("AMBI", "ITI", "ISI", "NSI", "NQI1", "Bentix", "Shannon")
  dimnames(indices)[[2]] <- dimnames(data_)[[2]][6:dim(data_)[2]]
  
  indices["Shannon",] <- output_shannon
  
  # for each SAMPLE
  for (i in dimnames(indices)[[2]])
  {
    # AMBI
    # keep only OTU with abundances
    ambi <- subset(data_, data_[,i] !=0)
    ambi <- subset(ambi, is.na(ambi[,i]) == FALSE ) # is the same as above... added the one below (we want only to compute with species with an index value
    ambi <- subset(ambi, is.na(ambi$AMBI) == FALSE )
    
    # normalize into %
    ambi[,i] <- as.numeric(ambi[,i])*100/sum(as.numeric(ambi[,i]))
    # creating the values and 0 by default
    val2 <-  val3  <-  val4  <-  val5 <- 0
    gp2 <-  gp3  <-  gp4  <-  gp5 <- 0
    # for each OTU in the sample with abundances
    for (j in dimnames(ambi)[[1]])
    {
      # if there is a value for AMBI
      if (is.na(ambi[j,"AMBI"]) == FALSE)
      {
        # for AMBIE value 2 (weights starts at 2)
        if (ambi[j,"AMBI"] == 2)
        {
          val2 <- 1.5 * ambi[j,i]
          gp2 <- gp2 + val2
        } else { val2 <- 0}
        if (ambi[j,"AMBI"] == 3)
        {
          val3 <- 3 * ambi[j,i]
          gp3 <- gp3 + val3
        } else { val3 <- 0}
        if (ambi[j,"AMBI"] == 4)
        {
          val4 <- 4.5 * ambi[j,i]
          gp4 <- gp4 + val4
        } else { val4 <- 0}
        if (ambi[j,"AMBI"] == 5)
        {
          val5 <- 6 * ambi[j,i]
          gp5 <- gp5 + val5
        } else { val5 <- 0}
      }
    }
    
    ambi_value <- (gp2 + gp3 + gp4 + gp5) /100
    
    # NSI
    nsi <- subset(data_, data_[,i] !=0)
    nsi <- subset(nsi, nsi$NSI > 0)
    nsi <- sum(as.numeric(nsi[,i]) * nsi[,"NSI"])/sum(as.numeric(nsi[,i]))
    
    # ISI 2012
    isi <- subset(data_, data_[,i] !=0)
    isi <- subset(isi, isi$ISI_value > 0)
    isi <- sum(isi[,"ISI_value"])/dim(isi)[1]
    
    # ITI index
    iti <- subset(data_, data_[,i] !=0)
    iti <- subset(iti, iti$ITI_GROUP > 0)
    n1 <- subset(iti, iti$ITI_GROUP ==1)
    n2 <- subset(iti, iti$ITI_GROUP ==2)
    n3 <- subset(iti, iti$ITI_GROUP ==3)
    n4 <- subset(iti, iti$ITI_GROUP ==4)
    iti <- 100-100/3*(sum(as.numeric(n2[,i])) + 2*sum(as.numeric(n3[,i])) + 3*sum(as.numeric(n4[,i])))/(sum(as.numeric(n1[,i])) +sum(as.numeric(n2[,i]))+sum(as.numeric(n3[,i]))+sum(as.numeric(n4[,i])))
    
    ## NQI1
    ambi_raw <- subset(data_, data_[,i] !=0)
    nqi1 <- 0.5*(1-ambi_value/7)+0.5*((log(dim(ambi_raw)[1])/log(log(sum(as.numeric(ambi_raw[,i])))))/2.7)*sum(as.numeric(ambi_raw[,i]))/(sum(as.numeric(ambi_raw[,i]))+5)
    
    ### BENTIX
    
    # keep only OTU with abundances
    bentix <- subset(data_, data_[,i] !=0)
    bentix <- subset(bentix, is.na(bentix[,i]) == FALSE ) # is the same as above... added the one below (we want only to compute with species with an index value
    bentix <- subset(bentix, is.na(bentix$Bentix) == FALSE )
    
    # normalize into %
    bentix[,i] <- as.numeric(bentix[,i])*100/sum(as.numeric(bentix[,i]))
    # creating the values and 0 by default
    val1 <-  val2  <-  val3  <-  0
    gp1 <-  gp2  <-  gp3  <- 0
    # for each OTU in the sample with abundances
    for (j in dimnames(bentix)[[1]])
    {
      # if there is a value for AMBI
      if (is.na(bentix[j,"AMBI"]) == FALSE)
      {
        # for AMBIE value 2 (weights starts at 2)
        if (bentix[j,"Bentix"] == 1)
        {
          val1 <- 6 * bentix[j,i]
          gp1 <- gp1 + val1
        } else { val2 <- 0}
        if (bentix[j,"Bentix"] == 2)
        {
          val2 <- 2 * bentix[j,i]
          gp2 <- gp2 + val2
        } else { val3 <- 0}
        if (bentix[j,"Bentix"] == 3)
        {
          val3 <- 2 * bentix[j,i]
          gp3 <- gp3 + val3
        } else { val3 <- 0}
      }
    }
    
    bentix_value <- (gp1 + gp2 + gp3) /100
    
    # then paste values in the indices table for output
    indices["AMBI", i] <- ambi_value
    indices["NSI", i] <- nsi
    indices["ISI", i] <- isi
    indices["ITI", i] <- iti
    indices["NQI1", i] <- nqi1
    indices["Bentix", i] <- bentix_value
    
  }
  return(list("taxa" = c("Found match:", cpt_found, " Not found:", cpt_not), "bi_values" = t(indices), "table_indices" = output, "OTU_list" = otu_id_list))
}
