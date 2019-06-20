#' genesetdnds
#'
#' Function to estimate global dN/dS values for a gene set when using whole-exome 
#' data. Global dN/dS values for a set of genes can also be obtained using dndscv 
#' (gene_list argument), but that option estimates the trinucleotide mutation rates 
#' exclusively from the list of genes of interest. This may be prefereable in large
#' datasets, but in small datasets, the genesetdnds option estimates global dN/dS
#' values for a set of genes while using all genes in the exome to fit the 
#' substitution model. Usage: genesetdnds(dndsout, gene_list).
#'
#' @author Inigo Martincorena (Wellcome Sanger Institute)
#' @details Martincorena I, et al. (2017) Universal patterns of selection in cancer and somatic tissues. Cell. 171(5):1029-1041.
#' 
#' @param dndsout Output object from dndscv. To generate a valid input object for this function, use outmats=T when running dndscv.
#' @param gene_list List of genes to restrict the analysis (gene set).
#' @param sm Substitution model (precomputed models are available in the data directory)
#'
#' @return 'genesetdnds' returns a list of objects:
#' @return globaldnds_geneset: Global dN/dS estimates in the gene set.
#' @return globaldnds_rest: Global dN/dS estimates in all other genes.
#'
#' @export

genesetdnds = function(dndsout, gene_list, sm = "192r_3w") {
    
    ## 1. Input
    
    if (is.null(dndsout$N)) { stop(sprintf("Invalid input: dndsout must be generated using outmats=T in dndscv.")) }
    
    allg = as.vector(dndsout$genemuts$gene) # All genes in the dndsout object
    nonex = gene_list[!(gene_list %in% allg)]
    if (length(nonex)>0) { 
        stop(sprintf("The following input gene names are not in the dndsout object: %s. To see the list of genes available use as.vector(dndsout$genemuts$gene).", paste(nonex,collapse=", ")))
    }
    if (length(gene_list)<2) {
        stop("The gene_list argument needs to contain at least two genes")
    }
    
    # Substitution model (The user can also input a custom substitution model as a matrix)
    if (length(sm)==1) {
        data(list=sprintf("submod_%s",sm), package="dndscv")
    } else {
        substmodel = sm
    }
    
    ## 2. Estimation of the global rate and selection parameters
    
    Lall = dndsout$L
    Nall = dndsout$N
    geneind = which(allg %in% gene_list) # Genes in the gene set
    L = rbind(apply(Lall[,,geneind], c(1,2), sum), apply(Lall[,,-geneind], c(1,2), sum))
    N = rbind(apply(Nall[,,geneind], c(1,2), sum), apply(Nall[,,-geneind], c(1,2), sum))
    
    # Subfunction: fitting substitution model
    
    fit_substmodel = function(N, L, substmodel) {
        
        l = c(L); n = c(N); r = c(substmodel)
        n = n[l!=0]; r = r[l!=0]; l = l[l!=0]
        
        params = unique(base::strsplit(x=paste(r,collapse="*"), split="\\*")[[1]])
        indmat = as.data.frame(array(0, dim=c(length(r),length(params))))
        colnames(indmat) = params
        for (j in 1:length(r)) {
            indmat[j, base::strsplit(r[j], split="\\*")[[1]]] = 1
        }
        
        model = glm(formula = n ~ offset(log(l)) + . -1, data=indmat, family=poisson(link=log))
        mle = exp(coefficients(model)) # Maximum-likelihood estimates for the rate params
        ci = exp(confint.default(model)) # Wald confidence intervals
        par = data.frame(name=gsub("\`","",rownames(ci)), mle=mle[rownames(ci)], cilow=ci[,1], cihigh=ci[,2])
        return(list(par=par, model=model))
    }
    
    syneqs = substmodel[,1] # Rate model for synonymous sites
    
    # Model 1: Fitting all mutation rates and the 3 global selection parameters
    
    rmatrix = array("",dim=dim(L))
    rmatrix[,1] = c(paste(syneqs,"*r_rel",sep=""), syneqs) # This adds an extra parameter (r_rel) to account for a different mutation rate (synonymous density) in the gene set
    rmatrix[,2] = c(paste(syneqs,"*r_rel*wmis_geneset",sep=""), paste(syneqs,"*wmis_rest",sep=""))
    rmatrix[,3] = c(paste(syneqs,"*r_rel*wnon_geneset",sep=""), paste(syneqs,"*wnon_rest",sep=""))
    rmatrix[,4] = c(paste(syneqs,"*r_rel*wspl_geneset",sep=""), paste(syneqs,"*wspl_rest",sep=""))
    poissout = fit_substmodel(N, L, rmatrix) # Original substitution model
    par1 = poissout$par

    # Model 2: Fitting all mutation rates and the 2 global selection parameters
    
    rmatrix = array("",dim=dim(L))
    rmatrix[,1] = c(paste(syneqs,"*r_rel",sep=""), syneqs) # This adds an extra parameter (r_rel) to account for a different mutation rate (synonymous density) in the gene set
    rmatrix[,2] = c(paste(syneqs,"*r_rel*wmis_geneset",sep=""), paste(syneqs,"*wmis_rest",sep=""))
    rmatrix[,3] = c(paste(syneqs,"*r_rel*wtru_geneset",sep=""), paste(syneqs,"*wtru_rest",sep=""))
    rmatrix[,4] = c(paste(syneqs,"*r_rel*wtru_geneset",sep=""), paste(syneqs,"*wtru_rest",sep=""))
    poissout = fit_substmodel(N, L, rmatrix) # Original substitution model
    par2 = poissout$par
    
    # Model 2: Fitting all mutation rates and the 2 global selection parameters
    
    rmatrix = array("",dim=dim(L))
    rmatrix[,1] = c(paste(syneqs,"*r_rel",sep=""), syneqs) # This adds an extra parameter (r_rel) to account for a different mutation rate (synonymous density) in the gene set
    rmatrix[,2] = c(paste(syneqs,"*r_rel*wall_geneset",sep=""), paste(syneqs,"*wall_rest",sep=""))
    rmatrix[,3] = c(paste(syneqs,"*r_rel*wall_geneset",sep=""), paste(syneqs,"*wall_rest",sep=""))
    rmatrix[,4] = c(paste(syneqs,"*r_rel*wall_geneset",sep=""), paste(syneqs,"*wall_rest",sep=""))
    poissout = fit_substmodel(N, L, rmatrix) # Original substitution model
    par3 = poissout$par
    
    globaldnds_geneset = rbind(par1, par2, par3)[c("wmis_geneset","wnon_geneset","wspl_geneset","wtru_geneset","wall_geneset"),-1]
    globaldnds_rest = rbind(par1, par2, par3)[c("wmis_rest","wnon_rest","wspl_rest","wtru_rest","wall_rest"),-1]
    return(list(globaldnds_geneset=globaldnds_geneset, globaldnds_rest=globaldnds_rest))
}