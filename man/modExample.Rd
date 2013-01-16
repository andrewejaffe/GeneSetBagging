\name{modExample}
\alias{modExample}
\docType{data}
\title{
Example model matrix
}
\description{
Example model matrix for GSE17913, with smoking status as the covariate 
of interest and adjusting for 10 surrogate variables
}
\usage{data(modExample)}
\format{
model.matrix, samples down the rows, predictors across the columns
}
\details{
}
\source{
http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17913
}
\references{
}
\examples{
data(modExample)
head(modExample)
}
\keyword{datasets}
