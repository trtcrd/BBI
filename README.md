# BBI

BBI is a R package for Benthic Biotic Indices calculation from composition data

It takes composition data with associated taxonomic assignments as input and output biotic indices values. 
Composition data can be infered unsing morpho-taxonomy or DNA based assignements.

## Requirements

To be able to use the BBI package, the following is required:
1. R or RStudio
2. vegan package (install.packages("vegan"))
3. BBI package (see the installation section below)
4. Composition data (see below for format)

## Installation

The BBI package can be installed in R or RStudio using the devtools package, by typing these commands in R :

```
install.packages("devtools")
library(devtools)
install_github("trtcrd/BBI")
```


## An example of composition data adequate for BBI 

taxa|Sample1|Sample2|Sample3|Sample4|...  
--- | --- | --- | --- | ---  | ---
Capitella capitata|11|204|100|299|...   
Paramphinome jeffreysii|3|2201|100|388|...   
Nematoda sp.|0|20|130|10|...   
Glycera alba|147|0|0|9|... 
Nemertea indet.|17|0|110|15|... 
...|...|...|...|...


## An example of output from BBI 

SampleID|AMBI|ITI|ISI|NSI|NQI1|Bentix|Shannon
--- | --- | --- | --- | ---  | --- | --- | ---
Sample1|3.179310|15.1489913|9.301250|20.90030|0.6147135|0|2.304191
Sample2|3.466875|3.5714286|9.571429|19.93133|0.5450725|0|2.755729
Sample3|3.269499|24.4289970|8.661667|18.00327|0.5796122|0|2.376514
Sample4|3.810226|18.5291309|9.056000|20.73982|0.5178981|0|2.069489
...|...|...|...|...




