#' sitednds
#'
#' Function to estimate site-wise dN/dS values and p-values against neutrality. This function is in testing, please interpret the results with caution. Also note that recurrent artefacts or SNP contamination can violate the null model and dominate the list of sites under apparent selection. A considerable number of significant synonymous sites may reflect a problem with the data. Be very critical of the results and if suspicious sites appear recurrently mutated consider refining the variant calling (e.g. using a better unmatched normal panel). In the future, this function may be extended to perform inferences at a codon level instead of at a single-base level.
#'
#' @author Inigo Martincorena (Wellcome Sanger Institute)
#' 
#' @param dndsout Output object from dndscv. To generate a valid input object for this function, use outmats=T when running dndscv.
#' @param min_recurr Minimum number of mutations per site to estimate site-wise dN/dS ratios. [default=2]
#' @param gene_list List of genes to restrict the analysis (only needed if the user wants to restrict the analysis to a subset of the genes in dndsout) [default=NULL, sitednds will be run on all genes in dndsout]
#' @param theta_option 2 options: "mle" (uses the MLE of the negative binomial size parameter) or "conservative" (uses the lower bound of the CI95). Values other than "mle" will lead to the conservative option. [default="mle"]
#' @param syn_drivers Vector with a list of known synonymous driver mutations to exclude from the background model [default="TP53:T125T"]. See Martincorena et al., Cell, 2017 (PMID:29056346).
#'
#' @return 'sitednds' returns a table of recurrently mutated sites and the estimates of the size parameter:
#' @return - recursites: Table of recurrently mutated sites with site-wise dN/dS values and p-values
#' @return - theta: Maximum likelihood estimate and CI95% for the size parameter of the negative binomial distribution. The lower this value the higher the variation of the mutation rate across sites not captured by the trinucleotide change or by variation across genes.
#' 
#' @export

sitednds = function(dndsout, min_recurr = 2, gene_list = NULL, theta_option = "mle", syn_drivers = "TP53:T125T") {
    
    ## 1. Fitting a negative binomial distribution at the site level considering the background mutation rate of the gene and of each trinucleotide
    message("[1] Site-wise negative binomial model accounting for trinucleotides and relative gene mutability...")
    
    # N and L matrices for synonymous mutations
    N = dndsout$N[,1,]
    L = dndsout$L[,1,]
    if (length(N)==0) { stop(sprintf("Invalid input: dndsout must be generated using outmats=T in dndscv.")) }
    if (nrow(dndsout$mle_submodel)!=195) { stop(sprintf("Invalid input: dndsout must be generated using the default trinucleotide substitution model in dndscv."))}
    
    # Restricting the analysis to an input list of genes
    if (!is.null(gene_list)) {
        g = as.vector(dndsout$genemuts$gene_name)
        nonex = gene_list[!(gene_list %in% g)]
        if (length(nonex)>0) {
            warning(sprintf("The following input gene names are not in dndsout input object and will not be analysed: %s.", paste(nonex,collapse=", ")))
        }
        dndsout$annotmuts = dndsout$annotmuts[which(dndsout$annotmuts$gene %in% gene_list), ]
        dndsout$genemuts = dndsout$genemuts[which(g %in% gene_list), ]
        N = N[,which(g %in% gene_list)]
        L = L[,which(g %in% gene_list)]
    }
    
    # Counts of observed mutations
    annotsubs = dndsout$annotmuts[which(dndsout$annotmuts$impact %in% c("Synonymous","Missense","Nonsense","Essential_Splice")),]
    annotsubs$trisub = paste(annotsubs$chr,annotsubs$pos,annotsubs$ref,annotsubs$mut,annotsubs$gene,annotsubs$aachange,annotsubs$impact,annotsubs$ref3_cod,annotsubs$mut3_cod,sep=":")
    annotsubs = annotsubs[which(annotsubs$ref!=annotsubs$mut),]
    freqs = sort(table(annotsubs$trisub), decreasing=T)
    
    # Relative mutation rate per gene
    # Note that this assumes that the gene order in genemuts has not been altered with respect to the N and L matrices, as it is currently the case in dndscv
    relmr = dndsout$genemuts$exp_syn_cv/dndsout$genemuts$exp_syn
    names(relmr) = dndsout$genemuts$gene_name
    
    # Substitution rates (192 trinucleotide rates, strand-specific)
    sm = setNames(dndsout$mle_submodel$mle, dndsout$mle_submodel$name)
    sm["TTT>TGT"] = 1 # Adding the TTT>TGT rate (which is arbitrarily set to 1 relative to t)
    sm = sm*sm["t"] # Absolute rates
    sm = sm[setdiff(names(sm),c("wmis","wnon","wspl","t"))] # Removing selection parameters
    sm = sm[order(names(sm))] # Sorting
    mat_trisub = array(sm, dim=c(192,ncol(L))) # Relative mutation rates by trinucleotide
    mat_relmr = t(array(relmr, dim=c(ncol(L),192))) # Relative mutation rates by gene
    R = mat_trisub * mat_relmr # Expected rate for each mutation type in each gene
    
    # Expanded vectors: full vectors of observed and expected mutations per site across all sites considered in dndsout
    rvec = rep(R, times=L) # Expanded vector of expected mutation counts per site
    nvec = array(0, length(rvec)) # Initialising the vector with observed mutation counts per site
    
    mutsites = read.table(text=names(freqs), header=0, sep=":", stringsAsFactors=F) # Frequency table of mutations
    colnames(mutsites) = c("chr","pos","ref","mut","gene","aachange","impact","ref3_cod","mut3_cod")
    mutsites$freq = freqs
    trindex = setNames(1:192, names(sm))
    geneindex = setNames(1:length(names(relmr)), names(relmr))
    mutsites$trindex = trindex[paste(mutsites$ref3_cod, mutsites$mut3_cod, sep=">")]
    mutsites$geneindex = geneindex[mutsites$gene]
    
    synsites = mutsites[which(mutsites$impact=="Synonymous"),]
    synsites = synsites[!(paste(synsites$gene,synsites$aachange,sep=":") %in% syn_drivers),]
    Lcum = array(cumsum(L), dim=dim(L)) # Cumulative L indicating the position to place a given mutation in the nvec following rvec
    synsites$vecindex = apply(as.matrix(synsites[,c("trindex","geneindex")]), 1, function(x) Lcum[x[1], x[2]]) # Index for the mutation
    synsites = synsites[order(synsites$vecindex), ] # Sorting by index in the nvec
    
    # Correcting the index when there are multiple synonymous mutations in the same gene and trinucleotide class
    s = synsites$vecindex
    for (j in 2:nrow(synsites)) {
        if (s[j]<=s[j-1]) {
            s[j] = s[j-1] + 1
        }
    }
    synsites$vecindex2 = s
    
    nvec[synsites$vecindex2] = synsites$freq # Expanded nvec for the negative binomial regression
    rvec = rvec * sum(nvec) / sum(rvec) # Minor correction ensuring that global observed and expected rates are identical
    
    # Estimation of overdispersion modelling rates per site as negative binomially distributed (i.e. quantifying uncertainty above Poisson using a Gamma) 
    # Using optimize appears to yield reliable results. Problems experienced with fitdistr, glm.nb and theta.ml. Consider using grid search if problems appear with optimize.
    nbin = function(theta, n=nvec, r=rvec) { -sum(dnbinom(x=n, mu=r, log=T, size=theta)) } # nbin loglik function for optimisation
    ml = optimize(nbin, interval=c(0,1000))
    theta_ml = ml$minimum
    
    
    # CI95% for theta using profile likelihood and iterative grid search (this yields slightly conservative CI95)
    
    grid_proflik = function(bins=5, iter=5) {
        for (j in 1:iter) {
            if (j==1) {
                thetavec = sort(c(0, 10^seq(-3,3,length.out=bins), theta_ml, theta_ml*10, 1e4)) # Initial vals
            } else {
                thetavec = sort(c(seq(thetavec[ind[1]], thetavec[ind[1]+1], length.out=bins), seq(thetavec[ind[2]-1], thetavec[ind[2]], length.out=bins))) # Refining previous iteration
            }
            
            proflik = sapply(thetavec, function(theta) -sum(dnbinom(x=nvec, mu=rvec, size=theta, log=T))-ml$objective) < qchisq(.95,1)/2 # Values of theta within CI95%
            ind = c(which(proflik[1:(length(proflik)-1)]==F & proflik[2:length(proflik)]==T)[1],
                    which(proflik[1:(length(proflik)-1)]==T & proflik[2:length(proflik)]==F)[1]+1)
            if (is.na(ind[1])) { ind[1] = 1 }
            if (is.na(ind[2])) { ind[2] = length(thetavec) }
        }
        return(thetavec[ind])
    }
    
    theta_ci95 = grid_proflik(bins=5, iter=5)
    
    
    
    ## 2. Calculating site-wise dN/dS ratios and P-values for recurrently mutated sites (P-values are based on the Gamma assumption underlying the negative binomial modelling)
    message("[2] Calculating site-wise dN/dS ratios and p-values...")
    
    recursites = mutsites[mutsites$freq>=min_recurr, c("chr","pos","ref","mut","gene","aachange","impact","ref3_cod","mut3_cod","freq")]
    recursites$mu = relmr[recursites$gene] * sm[paste(recursites$ref3_cod,recursites$mut3_cod,sep=">")]
    
    if (theta_option=="mle") {
        theta = theta_ml
    } else { # Conservative
        theta = theta_ci95[1]
    }
    
    recursites$dnds = recursites$freq / recursites$mu # Site-wise dN/dS (point estimate)
    recursites$pval = pnbinom(q=recursites$freq-0.5, mu=recursites$mu, size=theta, lower.tail=F)
    recursites = recursites[order(recursites$pval, -recursites$freq), ] # Sorting by p-val and frequency
    recursites$qval = p.adjust(recursites$pval, method="BH", n=sum(L))
    rownames(recursites) = NULL
    thetaout = setNames(c(theta_ml, theta_ci95), c("MLE","CI95low","CI95_high"))
    
    return(list(recursites=recursites, theta=thetaout))

}