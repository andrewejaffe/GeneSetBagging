\name{gseaBag}
\alias{gseaBag}
\title{Performs gene set bagging algorithm using the wilcoxon mean rank gene set test}
\usage{
  gseaBag(exprsMat, geneSetList, mod, alpha = 0.05,
    wilcoxTestAlternative = "up",
    geneids = rownames(exprsMat), coiIndex = 2, B = 100,
    seed = NULL)
}
\arguments{
  \item{exprsMat}{the gene x sample expression matrix}

  \item{geneSetList}{a list, where the names are the gene
  set names, and each element contains the gene ids for
  that set}

  \item{mod}{the model matrix for the experiment}

  \item{alpha}{the cutoff used for the hypergeometric test
  and for generating the replication proportion}

  \item{wilcoxTestAlternative}{the alternative test for the
  geneSetTest in the limma package}

  \item{geneids}{gene identifiers that correspond to the
  rows of the expression matrix}

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
  Performs gene set bagging algorithm using the wilcoxon
  mean rank gene set test
}
\examples{
#gseaBag(pExample, geneSetListExample, modExample, alpha=0.05,B=20)
}
\seealso{
  limma::geneSetTest
}

