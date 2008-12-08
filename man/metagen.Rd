\name{metagen}
\alias{metagen}
\title{Generic inverse variance meta-analysis}
\description{
  Fixed and random effects meta-analysis based on estimates and their
  standard errors; inverse variance weighting is used for pooling.
}
\usage{
metagen(TE, seTE, studlab, data=NULL, subset=NULL, sm="",
        level = 0.95, level.comb = level)
}
\arguments{
  \item{TE}{Estimate of treatment effect.}
  \item{seTE}{Standard error of treatment estimate.}
  \item{studlab}{An optional vector with study labels.}
  \item{data}{An optional data frame containing the study information.}
  \item{subset}{An optional vector specifying a subset of studies to be used.}
  \item{sm}{A character string indicating underlying summary measure,
    e.g., \code{"RD"}, \code{"RR"}, \code{"OR"}, \code{"AS"},
    \code{"MD"}, \code{"SMD"}.}
  \item{level}{The level used to calculate confidence intervals for
    individual studies.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
}
\details{
  Generic method for meta-analysis, only treatment estimates and their
  standard error are needed. The method is useful, e.g., for pooling of
  log hazard ratios. The inverse variance method is used for
  pooling. Random effects estimate is based on the DerSimonian-Laird
  method.
}
\value{
  An object of class \code{c("metagen", "meta")} with corresponding
  \code{print}, \code{summary}, \code{plot} function. The object is a
  list containing the following components:
  \item{TE, seTE, studlab,}{As defined above.}
  \item{sm}{}
  \item{w.fixed, w.random}{Weight of individual studies (in fixed and
    random effects model).}
  \item{TE.fixed, seTE.fixed}{Estimated overall treatment effect and
    standard error (fixed effect model).}
  \item{TE.random, seTE.random}{Estimated overall treatment effect and
    standard error (random effects model).}
  \item{k}{Number of studies combined in meta-analysis.}
  \item{Q}{Heterogeneity statistic.}
  \item{tau}{Square-root of between-study variance (moment estimator of
    DerSimonian-Laird).}
  \item{method}{Pooling method: \code{"Inverse"}.}
  \item{level}{The level used to calculate confidence intervals for
    individual studies.}
  \item{level.comb}{The level used to calculate confidence intervals for
    pooled estimates.}
  \item{call}{Function call.}
}
\references{
  Cooper H & Hedges LV (1994),
  \emph{The Handbook of Research Synthesis}.
  Newbury Park, CA: Russell Sage Foundation.
}
\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metabin}}, \code{\link{metacont}}, \code{\link{print.meta}}}

\examples{
data(Fleiss93)
meta1 <- metabin(event.e, n.e, event.c, n.c, data=Fleiss93, sm="RR", meth="I")
meta1

##
## Identical results by using the following commands:
##
meta1
metagen(meta1$TE, meta1$seTE, sm="RR")
}
\keyword{htest}
