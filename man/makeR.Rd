\name{makeR}
\alias{makeR}
\title{Calculates replication proportion, R}
\usage{
  makeR(nullp, alpha = 0.05)
}
\arguments{
  \item{nullp}{indexes pertaining to set (rows) by
  bootstrap number (columns)}

  \item{alpha}{desired significance for the replication
  proportion}
}
\value{
  numeric vector giving the replication proportion for each
  gene set
}
\description{
  Calculates replication proportion, R
}

