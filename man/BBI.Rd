\encoding{UTF-8}
\name{BBI}
\alias{BBI}
\title{BBI}
\usage{
BBI(data, log = FALSE)
}

\arguments{
  \item{data}{A data frame containing samples as columns and taxa as rows, with species (or last taxonomic rank) in the first column}
  \item{log}{Whether or not keeping the history of match searching in a separate file in the workinf directory. Default = F}
}

\description{
The BBI function searches the taxa of the composition data in the reference table for match.
The then compute benthic biotic indices. The function BBI returns a list of containing.
}

\value{

  Function \code{BBI} returns a list containing :

    \item{found}{The amount of taxa that matched an entry in the database and the amount that did not. }
    \item{BBI}{The BBI values per sample. }
    \item{table}{The subset of composition data that contains only taxa with at least a match in one of the BBI. }
    \item{taxa}{The list of taxa that matched an entry and the correspondant OTU, if from metabarcoding data data. }
}


\examples{
## Loading the exemple data
data("metab")
data("morpho")
## Computing BBI indices
BI_metab <- BBI(metab)
BI_metab <- BBI(morpho)
## print values
BI_metab$BBI
## print classes
BI_metab$BBIclass
}
