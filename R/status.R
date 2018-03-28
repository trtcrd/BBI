# function for converting the continuous values into discrete assessement

e$status.ambi <- function(z){
  if(is.na(z)) return("NA")
  if (z < 1.2) return("very good")
  else if (z >= 1.2 && z < 3.3)  return("good")
  else if (z >= 3.3 && z < 4.3)  return("moderate")
  else if (z >= 4.3 && z < 5.5)  return("bad")
  else if (z >= 5.5)  return("very bad")
}

e$status.nsi <- function(z){
  if(is.na(z)) return("NA")
  if (z < 10) return("very bad")
  else if (z >= 10 && z < 15) return("bad")
  else if (z >= 15 && z < 20) return("moderate")
  else if (z >= 20 && z < 25) return("good")
  else if (z >= 25) return("very good")
}

e$status.nqi1 <- function(z){
  if(is.na(z)) return("NA")
  if (z < 0.31) return("very bad")
  else if (z >= 0.31 && z < 0.49) return("bad")
  else if (z >= 0.49 && z < 0.63) return("moderate")
  else if (z >= 0.63 && z < 0.82) return("good")
  else if (z >= 0.82) return("very good")
}


e$status.isi <- function(z){
  if(is.na(z)) return("NA")
  if (z < 4.5) return("very bad")
  else if (z >= 4.5 && z < 6.1) return("bad")
  else if (z >= 6.1 && z < 7.5) return("moderate")
  else if (z >= 7.5 && z < 9.6) return("good")
  else if (z >= 9.6) return("very good")
}

e$status.shannon <- function(z){
  if(is.na(z)) return("NA")
  if (z < 0.9) return("very bad")
  else if (z >= 0.9 && z < 1.9) return("bad")
  else if (z >= 1.9 && z < 3) return("moderate")
  else if (z >= 3 && z < 4.8) return("good")
  else if (z >= 4.8) return("very good")
}


e$status.nEQR <- function(z){
  if(is.na(z)) return("NA")
  if (z < 0.2) return("very bad")
  else if (z >= 0.2 && z < 0.4) return("bad")
  else if (z >= 0.4 && z < 0.6) return("moderate")
  else if (z >= 0.6 && z < 0.8) return("good")
  else if (z >= 0.8) return("very good")
}
