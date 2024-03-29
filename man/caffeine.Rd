% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/caffeine.R
\docType{data}
\name{caffeine}
\alias{caffeine}
\title{Caffeine for daytime drowsiness}
\format{
A data frame with the following columns:
\tabular{rl}{
\bold{\emph{study}}\tab study label \cr
\bold{\emph{year}}\tab year of publication \cr
\bold{\emph{h.caf}}\tab Number of participants with headaches (caffeine group) \cr
\bold{\emph{n.caf}}\tab Number of participants (caffeine group) \cr
\bold{\emph{h.decaf}}\tab Number of participants with headaches (decaf group) \cr
\bold{\emph{n.decaf}}\tab Number of participants (decaf group) \cr
\bold{\emph{D1}}\tab Domain 1 of risk of bias 2 tool (RoB 2) \cr
\bold{\emph{D2}}\tab Domain 2 (RoB 2) \cr
\bold{\emph{D3}}\tab Domain 3 (RoB 2) \cr
\bold{\emph{D4}}\tab Domain 4 (RoB 2) \cr
\bold{\emph{D5}}\tab Domain 5 (RoB 2) \cr
\bold{\emph{rob}}\tab Overall RoB 2 assessment
}
}
\description{
Caffeine for daytime drowsiness (Cochrane Practice review)
}
\details{
Data come from the Cochrane Practice review on caffeine for daytime
drowsiness. Eight fictitous studies evaluate the risk of headaches
after drinking either caffeinated or decaffeinated coffee.
}
\examples{
oldset <- settings.meta("RevMan5")

data(caffeine)
head(caffeine)

m1 <- metabin(h.caf, n.caf, h.decaf, n.decaf, sm = "OR",
  data = caffeine, studlab = paste(study, year))

# Add risk of bias assessment to meta-analysis
m1 <- rob(D1, D2, D3, D4, D5, overall = rob, data = m1, tool = "rob2")

# Print risk of bias assessment
rob(m1)

# Forest plot with risk of bias assessment
forest(m1)

settings.meta(oldset)
}
\references{
Higgins JPT, Savović J, Page MJ, Sterne JA on behalf of the RoB2
Development Group (2019):
Revised Cochrane risk-of-bias tool for randomized trials.
\url{https://www.riskofbias.info/welcome/rob-2-0-tool}
}
\seealso{
\code{\link{metabin}}, \code{\link{rob}}
}
\keyword{datasets}
