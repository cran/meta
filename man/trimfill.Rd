\name{trimfill}
\alias{trimfill}
\title{Trim and fill method for meta-analysis}
\description{
  Trim and fill method for estimating and adjusting for the number and
  outcomes of missing studies in a meta-analysis.
}
\usage{
trimfill(x, seTE, left=NULL, ma.fixed=TRUE,
         type="L", n.iter.max=50,
         sm=NULL, studlab=NULL,
         silent=TRUE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{meta}, or estimated treatment
    effect in individual studies.}
  \item{seTE}{Standard error of estimated treatment effect (mandatory if
    \code{x} not of class \code{meta}).}
  \item{left}{A logical indicating whether studies are suppossed to be
    missing on the left or right side of the funnel plot. If NULL, the
    linear regression test for funnel plot symmetry (i.e., function
    metabias(..., meth="linreg")) is used to determine
    whether trials are missing on the left or right side.}
  \item{ma.fixed}{A logical indicating whether a fixed or random effects
    model is used to estimate the number of missing studies.}
  \item{type}{A character indicating which method is used to estimate
    the number of missing studies. Either \code{"L"} or \code{"R"}.}
  \item{n.iter.max}{Maximum number of iterations to estimate number of
    missing studies.}
  \item{sm}{An optional character string indicating underlying summary
    measure, e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"WMD"},
    \code{"SMD"}; ignored if \code{x} is of class \code{meta}.}
  \item{studlab}{An optional vector with study labels; ignored if
    \code{x} is of class \code{meta}.}
  \item{silent}{A logical indicating whether bassic information on iterations shown.}
}
\details{
  The trim and fill method can be used for estimating and adjusting for
  the number and outcomes of missing studies in a meta-analysis. The
  method relies on scrutiny of one side of a funnel plot for asymmetry
  assumed due to publication bias.

  The function \code{metagen} is called internally.
}
\value{
  An object of class \code{c("metagen", "meta", "trimfill")}. The object is a
  list containing the following components:
  \item{studlab, sm}{}
  \item{left, ma.fixed}{As defined above.}
  \item{type, n.iter.max}{}
  \item{TE, seTE}{Estimated treatment effect and standard error of individual studies.}
  \item{w.fixed, w.random}{Weight of individual studies (in fixed and
    random effects model).}
  \item{TE.fixed, seTE.fixed}{Estimated overall treatment effect and
    standard error (fixed effect model).}
  \item{TE.random, seTE.random}{Estimated overall treatment effect and
    standard error (random effects model).}
  \item{k}{Number of studies combined in meta-analysis.}
  \item{Q}{Heterogeneity statistic Q.}
  \item{tau}{Square-root of between-study variance (moment estimator of
    DerSimonian-Laird).}
  \item{method}{Pooling method: \code{"Inverse"}.}
  \item{call}{Function call.}
  \item{n.iter}{Actual number of iterations to estimate number of
    missing studies.}
  \item{trimfill}{A logical vector indicating studies that have been
    added by trim and fill method.}
  \item{k0}{Number of studies added by trim and fill.}
}
\references{
  Duval S & Tweedie R (2000), A nonparametric "Trim and Fill" method of
  accounting for publication bias in meta-analysis. \emph{Journal of the
    American Statistical Association}, \bold{95}, 89--98.

  Duval S & Tweedie R (2000), Trim and Fill: A simple funnel-plot-based
  method of testing and adjusting for publication bias in
  meta-analysis. \emph{Biometrics}, \bold{56}, 455--463.
}
\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}
\seealso{\code{\link{metagen}}, \code{\link{metabias}}, \code{\link{funnel}}}
\examples{
data(Fleiss93)
meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Fleiss93, sm="OR")
tf1 <- trimfill(meta1)
summary(tf1)
funnel(tf1, pch=ifelse(tf1$trimfill, 1, 16))
}
\keyword{htest}
