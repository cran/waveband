\name{power.sum}
\alias{power.sum}
\title{Sums of wavelets raised to integer powers}
\description{
Computes a sum expressed in terms of mother wavelets raised to the power two, three, or four. Either the exact solution or a faster approximation can be computed. 
}
\usage{
power.sum(alphas.wd, pow = 2, verbose = TRUE, type = "approx", plotfn = FALSE)
}
\arguments{
\item{alphas.wd}{A \code{wd.object}, the \code{D} component of which contains the coefficients of the powers of wavelets. The entry which would normalLy be the coefficient of the wavelet at scale j and location k is the coefficient of the same wavelet raised to the power pow.

If \code{pow=2}, then the overall scaling function coefficient is included in the sum, otherwise the C component is ignored completely. 

The \code{filter.number} and \code{/link{family}} components of \code{alphas.wd} are used to determine which wavelet is used.}
\item{pow}{The power to which the wavelets are raised; it can take values 2, 3, or 4.}
\item{verbose}{If \code{verbose=TRUE}, progress reports are printed while the sum is being evaluated.}
\item{type}{If \code{type="approx"}, the approximation is computed, if \code{type="exact"}, the exact solution is computed, and if \code{type="both"} both the exact and approximate solutions are found.}
\item{plotfn}{If \code{plotfn=TRUE}, the solution(s) found are plotted.}
}
\details{
For the approximate method, the powers of mother wavelets are represented by scaling functions (father wavelets) at a finer level. This is discussed in Barber, Nason, & Silverman (2001). 

Sums of powers of wavelets are used in the computation of posterior credible intervals for wavelet regression estimators; see the documentation for the function \code{\link{wave.band}} for more details. 
}
\value{
A vector containing the solution (either exact or approximate), or a list containing both solutions, depending on the value of "type". 
}
\section{SIDE EFFECTS}{
If \code{plotfn=TRUE}, the solution(s) found are plotted.
}
\seealso{
\code{\link{wave.band}}
}
\keyword{manip}
 
