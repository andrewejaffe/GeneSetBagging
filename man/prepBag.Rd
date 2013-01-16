\name{prepBag}
\alias{prepBag}
\title{Creates bagging indices across all iterations}
\usage{
  prepBag(coi, B)
}
\arguments{
  \item{coi}{covarite of interest, typically binary, but
  can be categorical}

  \item{B}{number of bagging iterations}
}
\value{
  number of subjects by number of bagging iterations matrix
}
\description{
  Creates bagging indices across all iterations
}
\examples{
#prepBag(rep(1:2,each=10), 10)
}

