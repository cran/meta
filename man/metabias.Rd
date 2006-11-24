\name{metabias}
\alias{metabias}
\title{Test for funnel plot asymmetry}
\description{
  Test for funnel plot asymmetry, based on rank correlation or linear
  regression method.
}
\usage{
metabias(x, seTE, TE.fixed, seTE.fixed,
         method = "rank",
         plotit = FALSE, correct = FALSE)
}
\arguments{
  \item{x}{An object of class \code{meta}, or estimated treatment
    effect in individual studies.}
  \item{seTE}{Standard error of estimated treatment effect (mandatory if
    \code{x} not of class \code{meta}).}
  \item{TE.fixed}{Overall treatment estimate (mandatory if \code{x} not
    of class \code{meta} and \code{method} = \code{"rank"}).}
  \item{seTE.fixed}{Standard error of overall treatment estimate
    (mandatory if \code{x} not of class \code{meta} and \code{method} =
    \code{"rank"}).}
  \item{method}{A character string indicating which test is to be
    used. Either \code{"rank"}, \code{"linreg"}, \code{"mm"} or
    \code{"count"}, can be abbreviated.}
  \item{plotit}{A logical indicating whether a plot should be produced
    for method \code{"rank"}, \code{"linreg"} or \code{"mm"}.}
  \item{correct}{A logical indicating whether a continuity corrected
    statistic is used for rank correlation methods \code{"rank"} and \code{"count"}.}
}
\details{
  If \code{method} is \code{"rank"}, the test statistic is based on the
  rank correlation between standardised treatment estimates and variance
  estimates of estimated treatment effects; Kendall's tau is used as
  correlation measure (Begg & Mazumdar, 1994). The test statistic
  follows a standard normal distribution. By default (if \code{correct}
  is FALSE), no continuity correction is utilised (Kendall & Gibbons,
  1990).
  
  If \code{method} is \code{"linreg"}, the test statistic is based on a
  linear regression of the standardised treatment effect (standard
  normal deviate) on the inverse of the standard error of the treatment
  estimate (Egger et al., 1997). The test statistic follows a t
  distribution with \code{number of studies - 2} degrees of freedom.
  
  If \code{method} is \code{"mm"}, the test statistic is based on a
  weighted linear regression using the method of moments estimator of
  the additive between-study variance component (method 3a in Thompson,
  Sharp, 1999). The test statistic follows a t distribution with
  \code{number of studies - 2} degrees of freedom.

  If \code{method} is \code{"count"}, the test statistic is based on the
  rank correlation between a standardised cell frequency and the inverse
  of the variance of the cell frequency; Kendall's tau is used as
  correlation measure (Schwarzer, 2006). The test statistic
  follows a standard normal distribution. By default (if \code{correct}
  is FALSE), no continuity correction is utilised (Kendall & Gibbons,
  1990).
}
\value{
  A list with class \code{"htest"} containing the following components:
  \item{estimate}{the estimated degree of funnel plot asymmetry, with
    name \code{"ks"} or \code{"bias"} corresponding to the method
    employed, i.e., rank correlation or regression method.}
  \item{statistic}{The value of the test statistic.}
  \item{parameter}{The degrees of freedom of the test statistic in the
    case that it follows a t distribution.}
  \item{p.value}{The p-value for the test.}
  \item{null.value}{The value of test statistic under the
    null hypothesis, always \code{0}.}
  \item{alternative}{A character string describing the alternative
    hypothesis.}
  \item{method}{A character string indicating what type of test was used.}
  \item{data.name}{A character string giving the names of the data.}
}
\references{
  Begg CB & Berlin JA (1994), Operating characteristics of a
  rank correlation test for publication bias. \emph{Biometrics},
  \bold{50}, 1088--1101.

  Kendall M & Gibbons JD (1990),
  \emph{Rank Correlation Methods}.
  London: Edward Arnold.

  Egger M, Smith GD, Schneider M & Minder C (1997),
  Bias in meta-analysis detected by a simple,
  graphical test. \emph{British Medical Journal}, \bold{315}, 629--634.

  Schwarzer G, Antes G, Schumacher M (2006),
  A test for publication bias in meta-analysis with sparse binary data.
  \emph{Statistics in Medicine},
  DOI: 10.1002/sim.2588, http://dx.doi.org/10.1002/sim.2588

  Thompson SG, Sharp, SJ (1999),
  Explaining heterogeneity in meta-analysis: A comparison of methods,
  \emph{Statistics in Medicine}, \bold{18}, 2693--2708.
}
\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{funnel}}, \code{\link{metabin}}, \code{\link{metacont}}, \code{\link{metagen}}}

\examples{
data(Olkin95)
meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Olkin95, subset=c(41,47,51,59),
                 sm="RR", meth="I")

metabias(meta1)
metabias(meta1, correct=TRUE)

metabias(meta1, method="linreg")
metabias(meta1, method="linreg", plotit=TRUE)


metabias(meta1, method="count")

##
## Same result:
##
metabias(meta1, method="linreg")$p.value
metabias(meta1$TE, meta1$seTE, method="linreg")$p.value
}
\keyword{htest}
