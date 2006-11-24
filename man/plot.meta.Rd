\name{plot.meta}
\alias{plot.meta}
\title{Plot function for objects of class meta}
\description{
  Draws a forest plot in the active graphics window.
}
\usage{
plot.meta(x, byvar, bylab, sortvar, studlab = TRUE, level = 0.95,
          level.comb = level, comb.f = FALSE, comb.r = FALSE,
          text.f = "Fixed effect model", text.r = "Random effects model",
          lty.f = 2, lty.r = 3, xlab=NULL, xlim, ylim, lwd = 1, cex = 1,
          cex.comb = 1.2 * cex, cex.axis=cex, cex.lab=cex,
          log = ifelse(x$sm == "RR" | x$sm == "OR" | x$sm == "HR", "x", ""),
          axes = TRUE, allstudies = TRUE, weight = "fixed", scale.diamond = 1,
          scale.square =  1, col.i = "black", ...)
}
\arguments{
  \item{x}{An object of class \code{meta}).}
  \item{byvar}{An optional vector containing grouping information (must
    be of same length as \code{x$TE}).}
  \item{bylab}{A character string with a label for the grouping variable.}
  \item{sortvar}{An optional vector used to sort the individual studies
    (must be of same length as \code{x$TE}).}
  \item{studlab}{A logical indicating whether study labels should be
    printed in the graph.}
  \item{level}{The level used to calculate confidence intervals for
    individual studies.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{comb.f}{A logical indicating whether the pooled fixed effect
    estimate should be plotted.}
  \item{comb.r}{A logical indicating whether the pooled random effects
    estimate should be plotted.}
  \item{text.f}{A character string used in the plot to label the pooled fixed
    effects estimate.}
  \item{text.r}{A character string used in the plot to label the pooled random
    effects estimate.}
  \item{lty.f}{Line type (pooled fixed effect estimate).}
  \item{lty.r}{Line type (pooled random effects estimate).}
  \item{xlab}{A label for the x axis.}
  \item{xlim}{The x limits (min,max) of the plot.}
  \item{ylim}{The y limits (min,max) of the plot.}
  \item{lwd}{The line width.}
  \item{cex}{A numerical value giving the amount by which plotting text
    and symbols should be scaled relative to the default.}
  \item{cex.comb}{A numerical value giving the amount by which plotting
    text and symbols for pooled fixed and random effects estimates should be scaled.}
  \item{cex.axis}{The magnification to be used for axis annotation
    relative to the current setting of \code{cex}.}
  \item{cex.lab}{The magnification to be used for x and y labels
    relative to the current setting of \code{cex}.}
  \item{log}{A character string which contains \code{"x"} if the x axis
    is to be logarithmic (other values for \code{log} are not
    reasonable).}
  \item{axes}{A logical indicating whether the x axis should be drawn on
    the plot.}
  \item{allstudies}{A logical indicating whether studies with inestimable
    treatment effects should be plotted.}
  \item{weight}{A character string indicating which type of plotting
    symbols is to be used for individual treatment estimates. One of
    \code{"same"} and \code{"fixed"}, i.e., plot symbols have
    the same size for all studies or represent the weight of the
    study in the fixed effect model.}
  \item{scale.diamond}{A numerical value giving the amount by which the
    diamond representing pooled treatment effects should be scaled
    relative to the default.}
  \item{scale.square}{A numerical value giving the amount by which the
    square representing treatment effects in individual studies should
    be scaled relative to the default.}
  \item{\dots}{Graphical parameters as in \code{par} may also be
    passed as arguments.}
  \item{col.i}{The color for individual study results and confidence limits.}
}
\details{
  A forest plot, also called confidence interval plot, is drawn in the
  active graphics window. Sub-group analyses are conducted and displayed
  in the plot if \code{byvar} is not missing.   
}
\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}
\seealso{\code{\link{metabin}}, \code{\link{metacont}}, \code{\link{metagen}}}

\examples{
data(Olkin95)
meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Olkin95, subset=c(41,47,51,59),
                 sm="RR", meth="I")

plot(meta1)
plot(meta1, byvar=c(1,2,1,2), bylab="label")
plot(meta1, byvar=1:4, xlim=c(0.02, 10))
}
\keyword{hplot}
