\name{metacont}
\alias{metacont}
\title{Meta-analysis of continuous outcome data}
\description{
  Calculation of fixed and random effects estimates for meta-analyses
  with continuous outcome data; inverse variance weighting is used for
  pooling.
}
\usage{
metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, studlab,
         data=NULL, subset=NULL, sm="WMD")
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{n.e}{Number of observations in experimental group.}
  \item{mean.e}{Estimated mean in experimental group.}
  \item{sd.e}{Standard deviation in experimental group.}
  \item{n.c}{Number of observations in control group.}
  \item{mean.c}{Estimated mean in control group.}
  \item{sd.c}{Standard deviation in control group.}
  \item{studlab}{An optional vector with study labels.}
  \item{data}{An optional data frame containing the study information,
    i.e., n.e, mean.e, sd.e, n.c, mean.c, and n.c.}
  \item{subset}{An optional vector specifying a subset of studies to be used.}
  \item{sm}{A character string indicating which summary measure
    (\code{"WMD"} or \code{"SMD"}) is to be used for pooling of
    studies.}
}
\details{
  Calculation of fixed and random effects estimates for meta-analyses
  with continuous outcome data; inverse variance weighting is used for
  pooling. The DerSimonian-Laird estimate is used in the
  random effects model. For the summary measure \code{"SMD"},
  Hedges' adjusted g is utilised for pooling.

  The function \code{metagen} is called internally to calculate
  individual and overall treatment estimates and standard errors.  
}
\value{
  An object of class \code{c("metacont", "meta")} with corresponding
  \code{print}, \code{summary}, \code{plot} function. The object is a
  list containing the following components:
  \item{n.e, mean.e, sd.e,}{}
  \item{n.c, mean.c, sd.c,}{As defined above.}
  \item{studlab, sm}{}
  \item{TE, seTE}{Estimated treatment effect and standard error of individual studies.}
  \item{w.fixed, w.random}{Weight of indiviudal studies (in fixed and
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
  \item{call}{Function call.}
}
\references{
  Cooper H & Hedges LV (1994),
  \emph{The Handbook of Research Synthesis}.
  Newbury Park, CA: Russell Sage Foundation.
}
\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}
\seealso{\code{\link{metabin}}, \code{\link{metagen}}}
\examples{
data(Fleiss93cont)
meta1 <- metacont(n.e, mean.e, sd.e, n.c, mean.c, sd.c, data=Fleiss93cont, sm="SMD")
meta1

meta2 <- metacont(Fleiss93cont$n.e, Fleiss93cont$mean.e,
                  Fleiss93cont$sd.e,
                  Fleiss93cont$n.c, Fleiss93cont$mean.c,
                  Fleiss93cont$sd.c,
                  sm="SMD")
meta2
}
\keyword{htest}

