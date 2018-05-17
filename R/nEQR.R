#' nEQR function
#' @param data A data frame containing samples as row and BBI values as column.
#' @example nEQR(my_BBI$BBI)
#' @author Tristan Cordier
#' @export
#' nEQR()


#### nEQR calculation from BBI values
## the boundaries between classes are the normalized one from each of the indices

## create a "e" environment for storing the functions
e <-new.env()

## just a normalizing function used by nEQR
e$linMap <- function(x, from, to) {
  # Shifting the vector so that min(x) == 0
  x <- x - min(x)
  # Scaling to the range of [0, 1]
  x <- x / max(x)
  # Scaling to the needed amplitude
  x <- x * (to - from)
  # Shifting to the needed level
  x + from
}

## EQR functions for each BBI
e$eqr.ambi <- function(z){
  t <- e$linMap(c(0,z,6),0,1)
  return(t[2:(length(t)-1)])
}

e$eqr.nsi <- function(z){
  t <- e$linMap(c(0,z,31),0,1)
  return(t[2:(length(t)-1)])
}

e$eqr.isi <- function(z){
  t <- e$linMap(c(0,z,13),0,1)
  return(t[2:(length(t)-1)])
}

e$eqr.nqi1 <- function(z){
  t <- e$linMap(c(0,z,0.9),0,1)
  return(t[2:(length(t)-1)])
}

e$eqr.shannon <- function(z){
  t <- e$linMap(c(0,z,5.7),0,1)
  return(t[2:(length(t)-1)])
}

## normalized EQR functions for each BBI
e$nEQR.ambi <- function(z){
  out <- c(1:length(z))
  for (i in 1:length(z))
  {
    if (z[i] < 0.2) out[i] <- (z[i]-0)*(0.2-0)/(0.2-0) + 0
    else if (z[i] >= 0.2 && z[i] < 0.55) out[i] <- (z[i]-0.2)*(0.4-0.2)/(0.55-0.2)+0.2
    else if (z[i] >= 0.55 && z[i] < 0.7166667) out[i] <- (z[i]-0.55)*(0.6-0.4)/(0.7166667-0.55)+0.4
    else if (z[i] >= 0.7166667 && z[i] < 0.9166667) out[i] <- (z[i]-0.7166667)*(0.8-0.6)/(0.9166667-0.7166667)+0.6
    else if (z[i] >= 0.9166667 && z[i] <= 1) out[i] <- (z[i]-0.9166667)*(1-0.8)/(1-0.9166667)+0.8
  }
  return(1-out)
}

e$nEQR.nsi <- function(z){
  out <- c(1:length(z))
  for (i in 1:length(z))
  {
    if (z[i] < 0.3225806) out[i] <- (z[i]-0)*(0.2-0)/(0.2-0) + 0
    else if (z[i] >= 0.3225806 && z[i] < 0.4838710) out[i] <- (z[i]-0.3225806)*(0.4-0.2)/(0.4838710-0.3225806)+0.2
    else if (z[i] >= 0.4838710 && z[i] < 0.6451613) out[i] <- (z[i]-0.4838710)*(0.6-0.4)/(0.6451613-0.4838710)+0.4
    else if (z[i] >= 0.6451613 && z[i] < 0.8064516) out[i] <- (z[i]-0.6451613)*(0.8-0.6)/(0.8064516-0.6451613)+0.6
    else if (z[i] >= 0.8064516 && z[i] <= 1) out[i] <- (z[i]-0.8064516)*(1-0.8)/(1-0.8064516)+0.8
  }
  return(out)
}


e$nEQR.nqi1 <- function(z){
  out <- c(1:length(z))
  for (i in 1:length(z))
  {
    if (z[i] < 0.3444444) out[i] <- (z[i]-0)*(0.2-0)/(0.2-0) + 0
    else if (z[i] >= 0.3444444 && z[i] < 0.5444444) out[i] <- (z[i]-0.3444444)*(0.4-0.2)/(0.5444444-0.3444444)+0.2
    else if (z[i] >= 0.5444444 && z[i] < 0.7) out[i] <- (z[i]-0.5444444)*(0.6-0.4)/(0.7-0.5444444)+0.4
    else if (z[i] >= 0.7 && z[i] < 0.9111111) out[i] <- (z[i]-0.7)*(0.8-0.6)/(0.9111111-0.7)+0.6
    else if (z[i] >= 0.9111111 && z[i] <= 1) out[i] <- (z[i]-0.9111111)*(1-0.8)/(1-0.9111111)+0.8
  }
  return(out)
}

e$nEQR.isi <- function(z){
  out <- c(1:length(z))
  for (i in 1:length(z))
  {
    if (z[i] < 0.3461538) out[i] <- (z[i]-0)*(0.2-0)/(0.2-0) + 0
    else if (z[i] >= 0.3461538 && z[i] < 0.4692308) out[i] <- (z[i]-0.3461538)*(0.4-0.2)/(0.4692308-0.3461538)+0.2
    else if (z[i] >= 0.4692308 && z[i] < 0.5769231) out[i] <- (z[i]-0.4692308)*(0.6-0.4)/(0.5769231-0.4692308)+0.4
    else if (z[i] >= 0.5769231 && z[i] < 0.7384615) out[i] <- (z[i]-0.5769231)*(0.8-0.6)/(0.7384615-0.5769231)+0.6
    else if (z[i] >= 0.7384615 && z[i] <= 1) out[i] <- (z[i]-0.7384615)*(1-0.8)/(1-0.7384615)+0.8
  }
  return(out)
}

e$nEQR.shannon <- function(z){
  out <- c(1:length(z))
  for (i in 1:length(z))
  {
    if (z[i] < 0.1578947) out[i] <- (z[i]-0)*(0.2-0)/(0.2-0) + 0
    else if (z[i] >= 0.1578947 && z[i] < 0.3333333) out[i] <- (z[i]-0.1578947)*(0.4-0.2)/(0.3333333-0.1578947)+0.2
    else if (z[i] >= 0.3333333 && z[i] < 0.5263158) out[i] <- (z[i]-0.3333333)*(0.6-0.4)/(0.5263158-0.3333333)+0.4
    else if (z[i] >= 0.5263158 && z[i] < 0.8421053) out[i] <- (z[i]-0.5263158)*(0.8-0.6)/(0.8421053-0.5263158)+0.6
    else if (z[i] >= 0.8421053 && z[i] <= 1) out[i] <- (z[i]-0.8421053)*(1-0.8)/(1-0.8421053)+0.8
  }
  return(out)
}

## and now the main function to compute nEQR across several BBI

nEQR <- function (data) {
  ## if bentix indices is given (it is not included in the nEQR assessment..)
  message("Bentix and ITI are being removed (if any), because not included in the nEQR assessment regulations")
  BI_val <- data[,!(colnames(data)) %in% c("Bentix", "ITI")]
  ## now compute nEQR
  # preparing the output
  out_n <- cbind(BI_val, EQR = rep(0, nrow(BI_val)))
  for (j in colnames(BI_val))
  {
    if (j == "AMBI") out_n[,j] <- e$nEQR.ambi(e$eqr.ambi(BI_val[,j]))
    if (j == "ISI") out_n[,j] <- e$nEQR.isi(e$eqr.isi(BI_val[,j]))
    if (j == "NSI") out_n[,j] <- e$nEQR.nsi(e$eqr.nsi(BI_val[,j]))
    if (j == "NQI1") out_n[,j] <- e$nEQR.nqi1(e$eqr.nqi1(BI_val[,j]))
    if (j == "Shannon") out_n[,j] <- e$nEQR.shannon(e$eqr.shannon(BI_val[,j]))
  }
  # renaming the columns
  colnames(out_n) <- paste0("n", colnames(out_n))
  ## then average over the row for nEQR
  out_n[,"nEQR"] <- rowMeans(out_n[,!(colnames(out_n)) %in% c("nEQR")])

  ## now we can return the dicrete assessment for nEQR
  neqr_class  <- out_n[,"nEQR"]
  # getting the status
  for (i in 1:length(neqr_class)) neqr_class[i] <- e$status.nEQR(out_n[i,"nEQR"])

  ## preparing the output
  output <- list("nEQR" = out_n,
                 "nEQRclass" = cbind(nEQR = as.numeric(out_n[,"nEQR"]), nEQR_class = neqr_class))

  return(output)
}
