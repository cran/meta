\name{read.mtv}
\alias{read.mtv}
\title{Import RevMan data files (.mtv)}
\description{
  Reads a file created with RevMan 4.1 and creates a data frame from it.
}
\usage{
read.mtv(file, sep="$")
}
\arguments{
  \item{file}{The name of the file which the data are to be read from.}
  \item{sep}{The field separator character. Values on each line of the
          file are separated by this character.}
}
\details{
  Reads a file created with RevMan 4.1 (Menu: "File" - "Export" -
  "Analysis data file...") and creates a data frame from it.
}
\value{
  A data frame containing the following components:
  \item{comp.no}{Comparison number.}
  \item{outcome.no}{Outcome number.}
  \item{group.no}{Group number.}
  \item{studlab}{Study label.}
  \item{year}{Year of publication.}
  \item{event.e}{Number of events in experimental group.}
  \item{n.e}{Number of observations in experimental group.}
  \item{event.c}{Number of events in control group.}
  \item{n.c}{Number of observations in control group.}
  \item{mean.e}{Estimated mean in experimental group.}
  \item{sd.e}{Standard deviation in experimental group.}
  \item{mean.c}{Estimated mean in control group.}
  \item{sd.c}{Standard deviation in control group.}
  \item{O.E}{Observed minus expected (IPD analysis).}
  \item{V}{Variance of \code{O.E} (IPD analysis).}
  \item{order}{Ordering of studies.}
  \item{conceal}{Concealment of treatment allocation.}
  \item{grplab}{Group label.}
  \item{type}{Type of outcome. D = dichotomous, C = continuous, P = IPD.}
  \item{outclab}{Outcome label.}
  \item{graph.exp}{Graph label for experimental group.}
  \item{graph.cont}{Graph label for control group.}
  \item{label.exp}{Label for experimental group.}
  \item{label.cont}{Label for control group.}
  \item{complab}{Comparison label.}
}
\author{Guido Schwarzer \email{sc@imbi.uni-freiburg.de}}

\seealso{\code{\link{metabin}}, \code{\link{metacont}}, \code{\link{metagen}}}

\examples{
## Locate MTV-data file "FLEISS93.MTV" in sub-directory of package "meta"
##
filename <- paste(searchpaths()[seq(along=search())[search()==
                  "package:meta"]], "/data/FLEISS93.MTV", sep="")
##
fleiss93.cc <- read.mtv(filename)
}
\keyword{datagen}
