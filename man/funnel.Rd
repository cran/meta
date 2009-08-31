\name{funnel}
\alias{funnel}
\alias{radial}
\title{Plot to assess funnel plot asymmetry}
\description{
  Draw a funnel or radial plot to assess funnel plot asymmetry in the
  active graphics window.
}
\usage{
funnel(x, y, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL,
       comb.fixed=FALSE, axes=TRUE,
       pch=1, text=NULL, cex=1, col=NULL,
       log="", yaxis="se", sm=NULL,
       level=NULL, ...)

radial(x, y, xlim=NULL, ylim=NULL,
       xlab="Inverse of standard error",
       ylab="Standardised treatment effect (z-score)",
       comb.fixed=TRUE, axes=TRUE,
       pch=1, text=NULL, cex=1, col=NULL,
       level=NULL, ...)
}
\arguments{
  \item{x}{An object of class \code{meta}, or estimated treatment
    effect in individual studies.}
  \item{y}{Standard error of estimated treatment effect (mandatory if
    \code{x} not of class \code{meta}).}
  \item{xlim}{The x limits (min,max) of the plot.}
  \item{ylim}{The y limits (min,max) of the plot.}
  \item{xlab}{A label for the x axis.}
  \item{ylab}{A label for the y axis.}
  \item{comb.fixed}{A logical indicating whether the pooled fixed effect
    estimate should be plotted.}
  \item{axes}{A logical indicating whether axes should be drawn on the
    plot.}
  \item{pch}{The plotting symbol used for individual studies.}
  \item{text}{A character vector specifying the text to be used instead
    of plotting symbol.}
  \item{cex}{The magnification to be used for plotting symbol.}
  \item{col}{A vector with colour of plotting symbols.}
  \item{log}{A character string which contains \code{"x"} if the x axis
    is to be logarithmic, \code{"y"} if the y axis is to be logarithmic
    and \code{"xy"} or \code{"yx"} if both axes are to be logarithmic
    (applies only to function \code{funnel}).}
  \item{yaxis}{A character string indicating which type of weights are
    to be used. Either \code{"se"}, \code{"invvar"}, \code{"invse"}, or \code{"size"}
    (applies only to function \code{funnel}). }
  \item{sm}{A character string indicating underlying summary measure,
    e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"AS"}, \code{"MD"},
    \code{"SMD"} (applies only to function \code{funnel}).}
  \item{level}{The confidence level utilised in the plot. For the funnel
  plot, confidence limits are not drawn if \code{yaxis="size"}.}
  \item{\dots}{Graphical parameters as in \code{par} may also be
    passed as arguments.}
}
\details{
  A funnel plot or radial plot, also called Galbraith plot, is drawn in
  the active graphics window. If \code{comb.fixed} is TRUE, the pooled
  estimate of the fixed effect model is plotted. If \code{level} is not
  NULL, the corresponding confidence limits are drawn.

  In the funnel plot, if \code{yaxis} is \code{"se"}, the standard error
  of the treatment estimates is plotted on the y axis which is likely to
  be the best choice (Sterne & Egger, 2001). Other possible choices for
  \code{yaxis} are \code{"invvar"} (inverse of the variance),
  \code{"invse"} (inverse of the standard error), and \code{"size"}
  (study size).
}
\references{
  Galbraith RF (1988a),
  Graphical display of estimates having differing standard
  errors. \emph{Technometrics}, \bold{30}, 271--281.

  Galbraith RF (1988b),
  A note on graphical presentation of estimated odds ratios from several
  clinical trials. \emph{Statistics in Medicine}, \bold{7}, 889--894.
  
  Light RJ & Pillemer DB (1984),
  \emph{Summing Up. The Science of Reviewing Research}.
  Cambridge: Harvard University Press.

  Sterne JAC & Egger M (2001),
  Funnel plots for detecting bias in meta-analysis: Guidelines on choice
  of axis. \emph{Journal of Clinical Epidemiology}, \bold{54},
  1046--1055.
}
\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}
\seealso{\code{\link{metabias}}, \code{\link{metabin}}, \code{\link{metagen}}}

\examples{
data(Olkin95)
meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Olkin95, subset=c(41,47,51,59),
                 sm="RR", meth="I")

oldpar <- par(mfrow=c(2, 2))

funnel(meta1)
##
## Same result as code above:
##
funnel(meta1$TE, meta1$seTE, sm="RR")

funnel(meta1, comb.fixed=TRUE, level=0.95)

radial(meta1, level=0.95)

par(oldpar)
}
\keyword{hplot}
