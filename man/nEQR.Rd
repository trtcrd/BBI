\encoding{UTF-8}
\name{nEQR}
\alias{nEQR}
\title{nEQR}
\usage{
nEQR(data)
}

\arguments{
  \item{data}{A data frame containing samples as rows and BBI values as columns}
}

\description{
The nEQR function compute the normalized Ecological Quality Ratio (nEQR).
}

\value{

  Function \code{nEQR} returns a list containing :

    \item{nEQR}{The nEQR value. }
    \item{nEQR_class}{The dicrete ecological quality assessment. }

}

\examples{
## Loading the exemple data
data("metab")
data("morpho")
## Computing BBI indices
BI_metab <- BBI(metab)
BI_morpho <- BBI(morpho)
## And then computing nEQR
eqr <- nEQR(BI_metab$BBI)
## print nEQR values
eqr$nEQR
## print discrete status
eqr$nEQRclass
}
