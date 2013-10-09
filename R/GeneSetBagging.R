###############
# gene set bagging
# andrew jaffe
# updated 10/9/2013
###########


#' Calculates replication proportion, R
#'
#' @param nullp indexes pertaining to set (rows) by bootstrap number (columns) 
#' @param alpha desired significance for the replication proportion
#' @return numeric vector giving the replication proportion for each gene set
#' @export
makeR = function(nullp, alpha=0.05) rowMeans(nullp < alpha)


#' Creates bagging indices across all iterations
#'
#' @param coi covarite of interest, typically binary, but can be categorical
#' @param B number of bagging iterations
#' @return number of subjects by number of bagging iterations matrix
#' @examples
#' prepBag(rep(1:2,each=10), 10) 
#' @export
prepBag = function(coi,B) {
	oIndexes = split(1:length(coi), coi)
	tmp = lapply(oIndexes, function(x) replicate(B, sample(x, replace=TRUE)))
	ind = do.call("rbind", tmp)
	return(ind)
}

#' Performs gene set bagging algorithm using the hypergeometric test for 2x2 tables
#' 
#' @param exprsMat is the gene x sample expression matrix
#' @param geneSetList is a list, where the names are the gene set names, and each element contains the gene ids for that set
#' @param mod the model matrix for the experiment
#' @param alpha is the cutoff used for the hypergeometric test and for generating the replication proportion
#' @param geneids gene identifiers that correspond to the rows of the expression matrix
#' @param geneSetUniverse vector of all gene ids in the universe
#' @param coiIndex the index for the binary covariate of interest in the model matrix
#' @param B the number of bagging iterations
#' @param seed sets the seed for random sampling
#' @return R the replication proportion for each inputted gene set
#' @return seed the seed used for the random sampling, for reproducibility
#' @examples
#' data(pExample, geneSetListExample, modExample)
#' hyperBag(pExample, geneSetListExample, modExample, alpha=0.05,B=20)
hyperBag = function(exprsMat, geneSetList, mod, alpha=0.05,
	geneids = rownames(exprsMat), geneSetUniverse=unique(geneids),
	coiIndex=2, B=100,seed=NULL) {
	
	if(ncol(exprsMat) != nrow(mod)) stop("The number of samples in the gene matrix does not match the model matrix.\n")
	
	if(!is.list(geneSetList)) stop("The gene set list must be a list.\n")
	
	if(length(unique(mod[,coiIndex]))!=2) stop("Only binary covariates of interest are supported.\n")
	
	if(!is.list(geneSetList)) stop("The gene set list must be a list.\n")
		
	if(is.null(seed)) seed= runif(1,min=0,max=1e8)
	
	set.seed(seed)
	
	ind = prepBag(mod[,coiIndex], B)

	nullt = apply(ind, 2, function(i) {
		fit = lmFit(exprsMat[,i], mod[i,])
		eb = ebayes(fit)
		eb$t[,coiIndex]})

	degfree = nrow(mod) - ncol(mod)
	pmat = 2*pt(-abs(nullt),df=degfree)

	pvals = apply(pmat, 2, function(x) {
		sig = geneids[x < alpha]
		hyperp = sapply(geneSetList, function(y) {
				base = sum(y %in% geneSetUniverse)
				enr = sum(y %in% sig)
				thep = phyper(enr-1, base, length(geneSetUniverse)-base,
					length(sig),lower.tail = F)
				return(thep)
		})
		return(hyperp)
	})
	
	R = makeR(pvals, alpha)
	
	outList = list(R = R, seed=seed)
	return(outList)
}

#' Performs gene set bagging algorithm using the wilcoxon mean rank gene set test
#' 
#' @param exprsMat the gene x sample expression matrix
#' @param geneSetList a list, where the names are the gene set names, and each element contains the gene ids for that set
#' @param mod the model matrix for the experiment
#' @param alpha the cutoff used for the hypergeometric test and for generating the replication proportion
#' @param wilcoxTestAlternative the alternative test for the geneSetTest in the limma package
#' @param geneids gene identifiers that correspond to the rows of the expression matrix
#' @param coiIndex the index for the binary covariate of interest in the model matrix
#' @param B the number of bagging iterations
#' @param seed sets the seed for random sampling
#' @seealso \code{\link[limma]{geneSetTest}}
#' @return R the replication proportion for each inputted gene set
#' @return seed the seed used for the random sampling, for reproducibility
#' @examples
#' data(pExample, geneSetListExample, modExample)
#' gseaBag(pExample, geneSetListExample, modExample, alpha=0.05,B=20)
gseaBag = function(exprsMat, geneSetList, mod, alpha=0.05,
	wilcoxTestAlternative="up",	geneids = rownames(exprsMat),
	coiIndex=2, B=100,  seed=NULL) {
	
	if(ncol(exprsMat) != nrow(mod)) stop("The number of samples in the gene matrix does not match the model matrix.\n")
	
	if(!is.list(geneSetList)) stop("The gene set list must be a list.\n")
	
	if(length(unique(mod[,coiIndex]))!=2) stop("Only binary covariates of interest are supported.\n")
	
	if(!is.list(geneSetList)) stop("The gene set list must be a list.\n")
		
	if(is.null(seed)) seed= runif(1,min=0,max=1e8)
	
	set.seed(seed)
	
	ind = prepBag(mod[,coiIndex], B)
		
	nullt = apply(ind, 2, function(i) {
		fit = lmFit(exprsMat[,i], mod[i,])
		eb = ebayes(fit)
		eb$t[,coiIndex]})

	pvals = apply(nullt, 2, function(x) {
		
			wilcoxp = sapply(geneSetList, function(y) {
					Index = which(geneids %in% y)
					return(geneSetTest(Index, x, wilcoxTestAlternative))
			})
			return(wilcoxp)
		
		})
	
	R = makeR(pvals, alpha)
	
	outList = list(R = R, seed=seed)
	return(outList)
}

#' @name pExample
#' Example expression dataset
#' @description Example gene expression dataset, first 5000 probes with Entrez Gene IDs from  GEO dataset: GSE17913
#' @docType data
#' @format matrix, genes down rows, samples across columns
#' @source http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17913
NULL

#' @name modExample
#' Example model matrix
#' @description Example model matrix for GSE17913, with smoking status as the covariate of interest
#' and adjusting for 10 surrogate variables
#' @docType data
#' @usage modExample
#' @format model.matrix, samples down the rows, predictors across the columns
#' @source http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17913
NULL

#' @name geneSetListExample
#' First 100 GO gene sets
#' @description List containing the Entrez IDs corresponding to the first 100 GO gene sets
#' @docType data
#' @usage geneSetListExample
#' @format list, where elements contain Entrez Gene IDs, and names are gene set names
#' @source http://www.geneontology.org/
NULL