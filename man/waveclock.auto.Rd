\name{waveclock.auto}
\alias{waveclock.auto}
\title{
Run the waveclock function in another instance of R
}
\description{
Performs the waveclock function in another instance of R.
This may be desirable when running waveclock over many datasets in a loop.
}
\usage{
waveclock.auto( ... )
}
\arguments{
\item{...}{
Parameters passed to \code{\link{waveclock}}
}
}
\details{
The function returns the output from \code{\link{waveclock}}
}
\references{
"Practical Time-Frequency Analysis: Gabor and Wavelet Transforms with
an Implementation in S", by Rene Carmona, Wen L. Hwang and Bruno
Torresani, Academic Press, 1998.
\url{http://sgdp.iop.kcl.ac.uk/tprice/software.html}
}
\author{T.S.Price}
\seealso{
\code{\link{waveclock}}
}
\keyword{ts}
