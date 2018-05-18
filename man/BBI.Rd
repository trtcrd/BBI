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
The \code{BBI} function searches the taxa of the composition data in the reference table for match.
The taxonomic assignement usually includes previous taxonomic ranks, e.g. Kingdom;Phylum;Class;Order;Family;Genus;species.
The taxonomic path must be separated by a semicolon ';' in order properly parsed by the function.
The \code{BBI} function then compute benthic biotic indices. The function BBI returns a list of containing.
}

\value{

  Function \code{BBI} returns a list containing :

    \item{found}{The amount of taxa that matched an entry in the database and the amount that did not. }
    \item{BBI}{The BBI values per sample. }
    \item{table}{The subset of composition data that contains only taxa with at least a match in one of the BBI. }
    \item{taxa}{The list of taxa that matched an entry and the correspondant OTU, if from metabarcoding data data. }
}

\references{
Borja, A., Franco, J., Pérez, V., 2000. A Marine Biotic Index to Establish the Ecological Quality of Soft-Bottom Benthos Within European Estuarine and Coastal Environments. \emph{Mar. Pollut. Bull.} 40, 1100–1114. doi:10.1016/S0025-326X(00)00061-8

Rygg, B., 2013. Norwegian Sensitivity Index (NSI) for marine macroinvertebrates, and an update of Indicator Species Index (ISI). \emph{Norwegian Institute for Water Research}

Rygg, B., 2002. Indicator species index for assessing benthic ecological quality in marine waters of Norway. \emph{Norwegian Institute for Water Research}

}


\examples{
## Loading the exemple data
data("metab")
data("morpho")
## Computing BBI indices
BI_metab <- BBI(metab)
BI_morpho <- BBI(morpho)
## print values
BI_metab$BBI
## print classes
BI_metab$BBIclass
}
