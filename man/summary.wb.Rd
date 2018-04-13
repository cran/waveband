\name{summary.wb}
\alias{summary.wb}
\title{Print information about a wb object.
}
\description{Summarize information about a wb object.
}
\usage{
\method{summary}{wb}(object, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{wb object to summarize
}
  \item{\dots}{Other arguments (not used)
}
}
\details{Summarizes information about a wb object. 
}
\value{No specific value
}
\author{G. P. Nason
}

\seealso{\code{\link{plot.wb}}, \code{\link{wave.band}}
}
\examples{
#
set.seed(1)
tmp <- wave.band(rnorm(32))
summary(tmp)
#Wave.band credible bands object
#Bands produced for object in data component of length:  32
#Credible intervals are in the bands component
#Wave.band Bayesian hyperparameter alpha was:  0.5
#Wave.band Bayesian hyperparameter beta was:  1
#Wave.band Wavelet filter number was:  8
#Wave.band Wavelet family was:  DaubLeAsymm
#Type of input (data or test signal): data
#Rsnr (if applicable):  3
}
\keyword{manip}
