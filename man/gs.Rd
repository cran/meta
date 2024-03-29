% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gs.R
\name{gs}
\alias{gs}
\title{Get default for a meta-analysis setting.}
\usage{
gs(x = NULL, unname = NULL)
}
\arguments{
\item{x}{A character string or vector with setting name(s).}

\item{unname}{A logical whether to remove names from attributes.}
}
\description{
Get default for a meta-analysis setting in R package \bold{meta}.
}
\details{
This function can be used to get the default for a meta-analysis
setting defined using \code{\link{settings.meta}}.

This function is primarily used to define default settings in
meta-analysis functions, e.g., \code{\link{metabin}} or
\code{\link{metacont}}. A list of all arguments with current
settings is printed using the command
\code{settings.meta("print")}.
}
\examples{
# Get default setting for confidence interval of random effects
# model
#
gs("method.random.ci")

# Get default setting for summary measure in metabin()
#
gs("smbin")

}
\seealso{
\code{\link{settings.meta}}
}
\author{
Guido Schwarzer \email{guido.schwarzer@uniklinik-freiburg.de}
}
