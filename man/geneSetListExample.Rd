\name{geneSetListExample}
\alias{geneSetListExample}
\docType{data}
\title{
First 100 GO gene sets
}
\description{
List containing the Entrez IDs corresponding to the first 100 GO gene sets
}
\usage{data(geneSetListExample)}
\format{
list, where elements contain Entrez Gene IDs, and names are gene set names
}
\details{
}
\source{
http://www.geneontology.org/
}
\references{
}
\examples{
data(geneSetListExample)
lapply(geneSetListExample[1:5],head)
}
\keyword{datasets}
