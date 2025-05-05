# function for converting the continuous values into discrete assessement. 

# Added an extra layer of security. Wrong values larger than possible or reasonable for an given index 
# will appear now as 'NA' (instead to be dangerously classified as 'excelent quality').

# FIXME. Adjust better the maximum values allowed for each index if needed.

quality <- c("terrible","bad","moderate","good","excellent")  
# deliberately avoiding names of levels with more than one word

e$status.ambi <-function(z) {
   if(is.na(z)) {return ("NA")}
   else {cut(z, breaks = c(0,1.2,3.3,4.3,5.5,100), include.lowest=TRUE,labels=quality)}
   }  # check upper limit 

e$status.nsi <-function(z) {
   if(is.na(z)) {return ("NA")}
   else {cut(z, breaks = c(0,10,15,20,25,100), include.lowest=TRUE,labels=quality)}
   }  # check max 

e$status.nqi1 <-function(z) {
   if(is.na(z)) {return ("NA")}
   else {cut(z, breaks = c(0,0.31,0.49,0.63,0.82,1), include.lowest=TRUE,labels=quality)}
   }  # check max

e$status.isi <-function(z) {
   if(is.na(z)) {return ("NA")}
   else {cut(z, breaks = c(0,4.5,6.1,7.5,9.6,10), include.lowest=TRUE,labels=quality)}
   }  # check max

e$status.shannon <-function(z) {
   if(is.na(z)) {return ("NA")}
   else {cut(z, breaks = c(0,0.9,1.9,3,4.8,1000000), include.lowest=TRUE,labels=quality)}
   }  # shannon index does not have a defined upper limit. We assume than is lower than one million in any case.

e$status.nEQR <-function(z) {
   if(is.na(z)) {return ("NA")}
   else {cut(z, breaks = c(0,0.2,0.4,0.6,0.8,1), include.lowest=TRUE,labels=quality)}
   }  # check than max allowed = 1
