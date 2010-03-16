\name{plot.wb}
\alias{plot.wb}
\title{Plots output from wave.band. }
\description{
Plots the output from \code{\link{wave.band}}.
}
\usage{
\method{plot}{wb}(x, col=FALSE, ...)
}
\arguments{
\item{x}{Output list from the function \code{\link{wave.band}}.}
\item{col}{Specifies whether the figure plotted is in colour or black and white.}
\item{...}{Any other arguments.}
}
\details{
The function \code{\link{wave.band}} offers a plotting option. This function will either reproduce the plot made by \code{\link{wave.band}} from that function's output, or produce a different plot which reproduces better in black and white.  
}
\value{
A plot is produced on the current graphics device 
}
\seealso{
\code{\link{wave.band}}
}
\keyword{manip}


