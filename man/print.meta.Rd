\name{print.meta}
\alias{print.meta}
\alias{summary.meta}
\alias{print.summary.meta}
\title{Print and summary method for objects of class meta}
\description{
  Print and summary method for objects of class \code{meta}.
}
\usage{
print.meta(x, sortvar, level=0.95, level.comb=level,
           details=FALSE, ma=TRUE, digits=max(4, .Options$digits - 3), ...)

summary.meta(object, byvar, bylab, bystud=FALSE, level.comb=0.95, ...)

print.summary.meta(x, digits = max(3, .Options$digits - 3),
                   print.byvar = TRUE, ...)
}
\arguments{
  \item{x}{An object of class \code{meta} or \code{summary.meta}.}
  \item{object}{An object of class \code{meta}.}
  \item{sortvar}{An optional vector used to sort the individual studies
    (must be of same length as \code{x$TE}).}
  \item{level}{The level used to calculate confidence intervals for
    individual studies.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{details}{A logical indicating whether further details of
    individual studies should be printed.}
  \item{ma}{A logical indicating whether the summary results of the
    meta-analysis should be printed.}
  \item{byvar}{An optional vector containing grouping information (must
    be of same length as \code{x$TE}).}
  \item{bylab}{A character string with a label for the grouping variable.}
  \item{bystud}{A logical indicating whether results of individual
    studies should be printed by grouping variable.}
  \item{digits}{Minimal number of significant digits, see \code{print.default}.}
  \item{print.byvar}{A logical indicating whether the name of the grouping
    variable should be printed in front of the group labels.}
  \item{...}{other arguments}
}
\value{
  A list is returned by the function \code{summary.meta} with the
  following elements:
  \item{fixed}{Results for fixed effect model (a list with elements TE,
    seTE, lower, upper, z, p, level).}
  \item{random}{Results for random effects model (a list with elements
    TE, seTE, lower, upper, z, p, level).}
  \item{random}{A list with elements TE, lower, upper, seTE, z, p, level,
    scale.}
  \item{k}{Number of studies combined in meta-analysis.}
  \item{Q}{Heterogeneity statistic Q.}
  \item{tau}{Square-root of between-study variance (moment estimator of
    DerSimonian-Laird).}
  \item{H}{Heterogeneity statistic H  (a list with elements TE, lower, upper).}
  \item{I2}{Heterogeneity statistic I2 (a list with elements TE, lower, upper).}
  \item{k.all}{Total number of trials.}
  \item{Q.CMH}{Cochrane-Mantel-Haenszel heterogeneity statistic.}
  \item{sm}{A character string indicating underlying summary measure.}
  \item{method}{A character string with the pooling method.}
  \item{call}{Function call.}
  \item{ci.lab}{Label for confidence interval.}
  \item{within}{Results within groups (a list with elements TE, seTE,
    lower, upper, z, p, level) - if \code{byvar} is not missing.}
  \item{k.w}{Number of studies combined within groups - if \code{byvar} is not missing.}
  \item{Q.w}{Heterogeneity statistic Q within groups - if \code{byvar} is not missing.}
  \item{bylab}{Label for grouping variable - if \code{byvar} is not missing.}
  \item{by.levs}{Levels of grouping variable - if \code{byvar} is not missing.}
}
\references{
  Cooper H & Hedges LV (1994),
  \emph{The Handbook of Research Synthesis}.
  Newbury Park, CA: Russell Sage Foundation.
}
\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}
\seealso{\code{\link{metabin}}, \code{\link{metacont}}, \code{\link{metagen}}}

\examples{
data(Fleiss93cont)
meta1 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, data=Fleiss93cont, sm="SMD")
summary(meta1)
summary(meta1, byvar=c(1,2,1,1,2), bylab="label")
}
\keyword{print}
