\name{metacum}
\alias{metacum}
\title{Cumulative meta-analysis}
\description{
  Performs a cumulative meta-analysis.
}
\usage{
metacum(x, pooled="fixed", sortvar)
}
\arguments{
  \item{x}{An object of class \code{meta}.}
  \item{pooled}{A character string indicating whether a fixed or random
    effects model is used for pooling. Either \code{"fixed"} or
    \code{"random"}, can be abbreviated.}
  \item{sortvar}{An optional vector used to sort the individual studies
    (must be of same length as \code{x$TE}).}
}
\details{
  A cumulative meta-analysis is performed. Studies are included
  sequentially as defined by \code{sortvar}.
}
\value{
  An object of class \code{c("metacum", "meta")} with corresponding
  \code{print}, \code{plot} function. The object is a
  list containing the following components:
  \item{TE, seTE}{Estimated treatment effect and standard error of
    pooled estimate in cumulative meta-analyses.}
  \item{studlab}{Study label describing addition of studies.}
  \item{sm}{Summary measure.}
  \item{method}{Method used for pooling.}
  \item{k}{Number of studies combined in meta-analysis.}
  \item{pooled}{As defined above.}

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
meta1 <- metabin(event.e, n.e, event.c, n.c,
                 data=Fleiss93, studlab=study,
                 sm="RR", meth="I")
meta1

metacum(meta1)
metacum(meta1, pooled="random")
plot(metacum(meta1, pooled="random"))
}
\keyword{htest}

