\encoding{UTF-8}
\name{morpho}
\alias{morpho}
\docType{data}
\title{Subset of morphologic inventories data}
\usage{
data("morpho")
}


\description{
The \code{morpho} dataset is an example of morpho-taxonomic inventories table produced from microscopical analysis.
Taxa are in rows and samples are in columns. It contains taxonomic assignements of taxa in the first column.

}


\examples{
## Loading the exemple data
data("morpho")
## Computing BBI indices
BI_morpho <- BBI(morpho)
## And then computing nEQR
eqr <- nEQR(BI_morpho$BBI)
## print nEQR values
eqr$nEQR
## print discrete status
eqr$nEQRclass
}
\keyword{datasets}
