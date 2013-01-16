\name{hyperBag}
\alias{hyperBag}
\title{Performs gene set bagging algorithm using the hypergeometric test for 2x2 tables}
\usage{
  hyperBag(exprsMat, geneSetList, mod, alpha = 0.05,
    geneids = rownames(exprsMat),
    geneSetUniverse = unique(geneids), coiIndex = 2,
    B = 100, seed = NULL)
}
\arguments{
  \item{exprsMat}{is the gene x sample expression matrix}

  \item{geneSetList}{is a list, where the names are the
  gene set names, and each element contains the gene ids
  for that set}

  \item{mod}{the model matrix for the experiment}

  \item{alpha}{is the cutoff used for the hypergeometric
  test and for generating the replication proportion}

  \item{geneids}{gene identifiers that correspond to the
  rows of the expression matrix}

  \item{geneSetUniverse}{vector of all gene ids in the
  universe}

  \item{coiIndex}{the index for the binary covariate of
  interest in the model matrix}

  \item{B}{the number of bagging iterations}

  \item{seed}{sets the seed for random sampling}
}
\value{
  R the replication proportion for each inputted gene set

  seed the seed used for the random sampling, for
  reproducibility
}
\description{
  Performs gene set bagging algorithm using the
  hypergeometric test for 2x2 tables
}
\examples{
# hyperBag(pExample, geneSetListExample, modExample, alpha=0.05,B=20)
}

