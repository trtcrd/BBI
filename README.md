# BBI

Set of functions to calculate Benthic Biotic Indices from composition data, obtained whether from morphotaxonomic inventories or sequencing data.

It takes composition data with associated taxonomic assignments as input and output biotic indices values.
It also return the ecological quality status for each pair of sample-BBI.
These BBI values can be used to calculate the normalized Ecological Quality Ratio (nEQR)

The composition can be derived from morpho-taxonomy or DNA based assignements.

## Requirements

To be able to use the BBI package, the following is required:
1. R or RStudio
2. vegan package
3. BBI package (see the installation section below)
4. Composition data (see below for format)

## Installation

The BBI package, in its lastest stable release, is available at:
https://cran.r-project.org/web/packages/BBI/index.html
It can be installed using the following command within R :
```
install.packages("BBI")
```

Alternatively, the developpement versino can be installed by typing these commands in R :

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


## Usage example

There is two functions:
BBI calculates and return for each sample the biotic indices values and their ecological quality status.
nEQR calculates and return for each sample the nEQR values and the associated normalized ecological quality status.

```
library(BBI)
# calculating BBI values and ecological quality status
my_BBI <- BBI(my_composition_data)
# calculating nEQR values and ecological quality status
my_nEQR <- nEQR(my_BBI$BBI)

```

## An example of output from BBI function

```
> my_BBI$found
[1] "Found match:" "101"          " Not found:"  "31"

> my_BBI$BBI
            AMBI       ISI      NSI      NQI1  Shannon        ITI    Bentix
Sample1 3.178476  9.571429 20.81954 0.5830994 2.093312 15.0084794 2.0249856
Sample2 3.471889  9.571429 19.93133 0.5284062 2.465243 21.4910281 2.4569632
Sample3 3.272532  9.180000 18.35880 0.5415376 2.181589 23.0795610 2.8546545
Sample4 4.777285  9.042857 20.92236 0.4173243 2.562996 23.3007691 2.9201504

> my_BBI$BBIclass
          AMBI        ISI         NSI        NQI1       Shannon
Sample1   "good"      "good"      "good"     "moderate" "moderate"
Sample2   "moderate"  "good"      "moderate" "moderate" "moderate"
Sample3   "good"      "good"      "moderate" "moderate" "moderate"
Sample4   "bad"       "good"      "good"     "bad"      "moderate"

```

## An example of output from nEQR function

```
> my_nEQR$nEQR
            nAMBI      nISI      nNSI     nNQI1  nShannon      nEQR
Sample1 0.6115737 0.7342886 0.6327817 0.5329992 0.4351477 0.5893582
Sample2 0.5656222 0.7342886 0.5972531 0.4548661 0.5027715 0.5709603
Sample3 0.6026160 0.6995857 0.5343520 0.4736251 0.4511980 0.5522754
Sample4 0.3204526 0.6874270 0.6368945 0.3192493 0.5205447 0.4969136

> my_nEQR$nEQR_class
        nEQR                nEQR_class
Sample1 "0.58935817792461"  "moderate"
Sample2 "0.570960301407608" "moderate"
Sample3 "0.552275367599034" "moderate"
Sample4 "0.496913603640647" "moderate"

```

## Paper and citation

Cordier T., Pawlowski J. BBI: an R package for the computation of Benthic Biotic Indices from composition data. Metabarcoding and Metagenomics 2: e25649, doi: 10.3897/mbmg.2.25649

## Version history

### version 0.3.0 ###

Updated the BBI reference table

### version 0.2.0 ###

Various code cleaning and optimizations, added help pages for package functions

### version 0.1.0 ###

First version
