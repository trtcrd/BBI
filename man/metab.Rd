\encoding{UTF-8}
\name{metab}
\alias{metab}
\docType{data}
\title{Subset of metabarcoding data}
\usage{
data("metab")
}


\description{
The \code{metab} dataset is an example of OTU table produced from metabarcoding data.
OTUs are in rows and samples are in columns. It contains taxonomic assignements of OTUs in the first column.
These assignement usually includes previous taxonomic ranks, e.g. Kingdom;Phylum;Class;Order;Family;Genus;species
The taxonomic path must be separated by a semicolon ';' in order the properly parse it.

}


\examples{
## Loading the exemple data
data("metab")
## Computing BBI indices
BI_metab <- BBI(metab)
## And then computing nEQR
eqr <- nEQR(BI_metab$BBI)
## print nEQR values
eqr$nEQR
## print discrete status
eqr$nEQRclass
}
\keyword{datasets}
