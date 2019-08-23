#' codondnds
#'
#' Function to estimate codon-wise dN/dS values and p-values against neutrality. To generate a valid RefCDS input object for this function, use the buildcodon function. Note that recurrent artefacts or SNP contamination can violate the null model and dominate the list of codons under apparent selection. Be very critical of the results and if suspicious codons appear recurrently mutated consider refining the variant calling (e.g. using a better unmatched normal panel).
#'
#' @author Inigo Martincorena (Wellcome Sanger Institute)
#' 
#' @param dndsout Output object from dndscv.
#' @param refcds RefCDS object annotated with codon-level information using the buildcodon function.
#' @param min_recurr Minimum number of mutations per codon to estimate codon-wise dN/dS ratios. [default=2]
#' @param gene_list List of genes to restrict the p-value and q-value calculations (Restricted Hypothesis Testing). Note that q-values are only valid if the list of genes is decided a priori. [default=NULL, codondnds will be run on all genes in dndsout]
#' @param codon_list List of hotspot codons to restrict the p-value and q-value calculations (Restricted Hypothesis Testing). Note that q-values are only valid if the list of codons is decided a priori. [default=NULL, codondnds will be run on all genes in dndsout]
#' @param theta_option 2 options: "mle" (uses the MLE of the negative binomial size parameter) or "conservative" (uses the lower bound of the CI95). Values other than "mle" will lead to the conservative option. [default="conservative"]
#' @param syn_drivers Vector with a list of known synonymous driver mutations to exclude from the background model [default="TP53:T125T"]. See Martincorena et al., Cell, 2017 (PMID:29056346).
#' @param method Overdispersion model: NB = Negative Binomial (Gamma-Poisson), LNP = Poisson-Lognormal (see Hess et al., BiorXiv, 2019). [default="NB"]
#' @param numbins Number of bins to discretise the rvec vector [default=1e4]. This enables fast execution of the LNP model in datasets of arbitrarily any size.
#'
#' @return 'codondnds' returns a table of recurrently mutated codons and the estimates of the size parameter:
#' @return - recurcodons: Table of recurrently mutated codons with codon-wise dN/dS values and p-values
#' @return - recurcodons_ext: The same table of recurrently mutated codons, but including additional information on the contribution of different changes within a codon.
#' @return - overdisp: Maximum likelihood estimate and CI95% for the overdispersion parameter (the size parameter of the negative binomial distribution or the sigma parameter of the lognormal distribution). The lower the size value or the higher the sigma value the higher the variation of the mutation rate across codons not captured by the trinucleotide change or by variation across genes.
#' @return - LL: Log-likelihood of the fit of the overdispersed model (see "method" argument) to all synonymous sites at a codon level.
#'
#' @export

codondnds = function(dndsout, refcds, min_recurr = 2, gene_list = NULL, codon_list = NULL, theta_option = "conservative", syn_drivers = "TP53:T125T", method = "NB", numbins = 1e4) {
    
    ## 1. Fitting an overdispersed distribution at the codon level considering the background mutation rate of the gene and of each trinucleotide
    message("[1] Codon-wise overdispersed model accounting for trinucleotides and relative gene mutability...")
    
    if (nrow(dndsout$mle_submodel)!=195) { stop("Invalid input: dndsout must be generated using the default trinucleotide substitution model in dndscv.") }
    if (is.null(refcds[[1]]$codon_impact)) { stop("Invalid input: the input RefCDS object must contain codon-level annotation. Use the buildcodon function to add this information.") }
    
    # Restricting refcds to genes in the dndsout object
    refcds = refcds[sapply(refcds, function(x) x$gene_name) %in% dndsout$genemuts$gene_name] # Only input genes
    
    # Restricting the analysis to an input list of genes
    if (!is.null(gene_list)) {
        g = as.vector(dndsout$genemuts$gene_name)
        # Correcting CDKN2A if required (hg19)
        if (any(g %in% c("CDKN2A.p14arf","CDKN2A.p16INK4a")) & any(gene_list=="CDKN2A")) {
            gene_list = unique(c(setdiff(gene_list,"CDKN2A"),"CDKN2A.p14arf","CDKN2A.p16INK4a"))
        }
        nonex = gene_list[!(gene_list %in% g)]
        if (length(nonex)>0) {
            warning(sprintf("The following input gene names are not in dndsout input object and will not be analysed: %s.", paste(nonex,collapse=", ")))
        }
        refaux = refcds[sapply(refcds, function(x) x$gene_name) %in% gene_list] # Only input genes
        numtests = sum(sapply(refaux, function(x) x$CDS_length))/3 # Number of codons in genes listed in gene_list
    } else {
        numtests = sum(sapply(refcds, function(x) x$CDS_length))/3 # Number of codons in all genes
    }
    
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
    
    # Annotated mutations per gene
    annotsubs = dndsout$annotmuts[which(dndsout$annotmuts$impact=="Synonymous"),]
    if (nrow(annotsubs)<2) {
        stop("Too few synonymous mutations found in the input. codondnds cannot run without synonymous mutations.")
    }
    annotsubs = annotsubs[!(paste(annotsubs$gene,annotsubs$aachange,sep=":") %in% syn_drivers),]
    annotsubs$codon = as.numeric(substr(annotsubs$aachange,2,nchar(annotsubs$aachange)-1)) # Numeric codon position
    annotsubs = split(annotsubs, f=annotsubs$gene)
    
    # Calculating observed and expected mutation rates per codon for every gene
    numcodons = sum(sapply(refcds, function(x) x$CDS_length))/3 # Number of codons in all genes
    nvec = rvec = array(NA, numcodons)
    pos = 1
    
    for (j in 1:length(refcds)) {
        
        nvec_syn = rvec_syn = rvec_ns = array(0,refcds[[j]]$CDS_length/3) # Initialising the obs and exp vectors
        gene = refcds[[j]]$gene_name
        sm_rel = sm * relmr[gene]
        
        # Expected rates
        ind = rep(1:(refcds[[j]]$CDS_length/3), each=9)
        syn = which(refcds[[j]]$codon_impact==1) # Synonymous changes
        ns = which(refcds[[j]]$codon_impact %in% c(2,3)) # Missense and nonsense changes
        
        aux = sapply(split(refcds[[j]]$codon_rates[syn], f=ind[syn]), function(x) sum(sm_rel[x]))
        rvec_syn[as.numeric(names(aux))] = aux
        
        aux = sapply(split(refcds[[j]]$codon_rates[ns], f=ind[ns]), function(x) sum(sm_rel[x]))
        rvec_ns[as.numeric(names(aux))] = aux
        
        # Observed mutations
        subs = annotsubs[[gene]]
        if (!is.null(subs)) {
            obs_syn = table(subs$codon)
            nvec_syn[as.numeric(names(obs_syn))] = obs_syn
        }
        
        rvec[pos:(pos+refcds[[j]]$CDS_length/3-1)] = rvec_syn
        nvec[pos:(pos+refcds[[j]]$CDS_length/3-1)] = nvec_syn
        pos = pos + refcds[[j]]$CDS_length/3
        
        refcds[[j]]$codon_rvec_ns = rvec_ns
        
        if (round(j/2000)==(j/2000)) { message(sprintf('    %0.3g%% ...', round(j/length(refcds),2)*100)) }
    }
    
    rvec = rvec * sum(nvec) / sum(rvec) # Small correction ensuring that global observed and expected rates are identical
    
    
    message("[2] Estimating overdispersion and calculating codon-wise dN/dS ratios...")
    
    # Estimation of overdispersion: Using optimize appears to yield reliable results. Problems experienced with fitdistr, glm.nb and theta.ml. Consider using grid search if problems appear with optimize.
    if (method=="LNP") { # Modelling rates per codon with a Poisson-Lognormal mixture
        
        lnp_est = fitlnpbin(nvec, rvec, theta_option = theta_option, numbins = numbins)
        theta_ml = lnp_est$ml$minimum
        theta_ci95 = lnp_est$sig_ci95
        LL = -lnp_est$ml$objective # LogLik
        thetaout = setNames(c(theta_ml, theta_ci95), c("MLE","CI95_high"))
        
    } else { # Modelling rates per codon as negative binomially distributed (i.e. quantifying uncertainty above Poisson using a Gamma)
        
        nbin = function(theta, n=nvec, r=rvec) { -sum(dnbinom(x=n, mu=r, log=T, size=theta)) } # nbin loglik function for optimisation
        ml = optimize(nbin, interval=c(0,1000))
        theta_ml = ml$minimum
        LL = -ml$objective # LogLik
        
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
        thetaout = setNames(c(theta_ml, theta_ci95), c("MLE","CI95low","CI95_high"))
    }
    
    
    ## 2. Calculating codon-wise dN/dS ratios and P-values for recurrently mutated codons
    
    # Theta option
    if (theta_option=="mle" | theta_option=="MLE") {
        theta = theta_ml
    } else { # Conservative
        message("    Using the conservative bound of the confidence interval of the overdispersion parameter.")
        theta = theta_ci95[1]
    }
    
    # Creating the recurcodons object
    annotsubs = dndsout$annotmuts[which(dndsout$annotmuts$impact %in% c("Missense","Nonsense")),]
    annotsubs$codon = substr(annotsubs$aachange,1,nchar(annotsubs$aachange)-1) # Codon position
    annotsubs$codonsub = paste(annotsubs$chr,annotsubs$gene,annotsubs$codon,sep=":")
    annotsubs = annotsubs[which(annotsubs$ref!=annotsubs$mut),]
    freqs = sort(table(annotsubs$codonsub), decreasing=T)

    recurcodons = read.table(text=names(freqs), header=0, sep=":", stringsAsFactors=F) # Frequency table of mutations
    colnames(recurcodons) = c("chr","gene","codon")
    recurcodons$freq = freqs
    
    # Gene RHT
    if (!is.null(gene_list)) {
        message("    Peforming Restricted Hypothesis Testing on the input list of a-priori genes")
        recurcodons = recurcodons[which(recurcodons$gene %in% gene_list), ] # Restricting the p-value and q-value calculations to gene_list
    }
    
    # Codon RHT
    if (!is.null(codon_list)) {
        message("    Peforming Restricted Hypothesis Testing on the input list of a-priori codons (numtests = length(codon_list))")
        mutstr = paste(recurcodons$gene,recurcodons$codon,sep=":")
        if (!any(mutstr %in% codon_list)) {
            stop("No mutation was observed in the restricted list of known hotspots. Codon-RHT cannot be run.")
        }
        recurcodons = recurcodons[which(mutstr %in% codon_list), ] # Restricting the p-value and q-value calculations to codon_list
        numtests = length(codon_list)
        
        # Calculating global dN/dS ratios at known hotcodons
        auxcodons = as.data.frame(do.call("rbind",strsplit(codon_list,split=":")), stringsAsFactors=F)
        auxcodons$V3 = as.numeric(substr(auxcodons$V2,2,nchar(auxcodons$V2)))
        auxcodons = auxcodons[auxcodons$V1 %in% names(relmr), ]
        colnames(auxcodons) = c("gene","codon","numcodon")
        auxcodons$mu = NA
        geneind = setNames(1:length(refcds), sapply(refcds, function(x) x$gene_name))
        for (j in 1:nrow(auxcodons)) {
            auxcodons$mu[j] = refcds[[geneind[auxcodons$gene[j]]]]$codon_rvec_ns[auxcodons$numcodon[j]] # Background non-synonymous rate for this codon
        }
        neutralexp = sum(auxcodons$mu) # Number of mutations expected at known hotspots expected under neutrality
        numobs = sum(recurcodons$freq) # Number observed
        poistest = poisson.test(numobs, T=neutralexp)
        globaldnds_knowncodons = setNames(c(numobs, neutralexp, poistest$estimate, poistest$conf.int), c("obs","exp","dnds","cilow","cihigh"))
        message(sprintf("    Mutations at known hotspots: %0.0f observed, %0.3g expected, obs/exp~%0.3g (CI95:%0.3g,%0.3g).", globaldnds_knowncodons[1], globaldnds_knowncodons[2], globaldnds_knowncodons[3], globaldnds_knowncodons[4], globaldnds_knowncodons[5]))
    }
    
    # Restricting the recurcodons output by min_recurr
    recurcodons = recurcodons[recurcodons$freq>=min_recurr, ] # Restricting the output to codons with min_recurr
    
    if (nrow(recurcodons)>1) {
    
        recurcodons$mu = NA
        codonnumeric = as.numeric(substr(recurcodons$codon,2,nchar(recurcodons$codon))) # Numeric codon position
        geneind = setNames(1:length(refcds), sapply(refcds, function(x) x$gene_name))
    
        for (j in 1:nrow(recurcodons)) {
            recurcodons$mu[j] = refcds[[geneind[recurcodons$gene[j]]]]$codon_rvec_ns[codonnumeric[j]] # Background non-synonymous rate for this codon
        }
        
        recurcodons$dnds = recurcodons$freq / recurcodons$mu # Codon-wise dN/dS (point estimate)
        
        if (method=="LNP") { # Modelling rates per codon with a Poisson-Lognormal mixture
            
            # Cumulative Lognormal-Poisson using poilog::dpoilog
            dpoilog = poilog::dpoilog
            ppoilog = function(n, mu, sig) {
                p = sum(dpoilog(n=floor(n+1):floor(n*10+1000), mu=log(mu)-sig^2/2, sig=sig))
                return(p)
            }
            
            message(sprintf("    Modelling substitution rates using a Lognormal-Poisson: sig = %0.3g (upperbound = %0.3g)", theta_ml, theta_ci95))
            recurcodons$pval = apply(recurcodons, 1, function(x) ppoilog(n=as.numeric(x["freq"])-0.5, mu=as.numeric(x["mu"]), sig=theta))
            
        } else { # Negative binomial model
            
            message(sprintf("    Modelling substitution rates using a Negative Binomial: theta = %0.3g (CI95:%0.3g,%0.3g)", theta_ml, theta_ci95[1], theta_ci95[2]))
            recurcodons$pval = pnbinom(q=recurcodons$freq-0.5, mu=recurcodons$mu, size=theta, lower.tail=F)
        }
        
        recurcodons = recurcodons[order(recurcodons$pval, -recurcodons$freq), ] # Sorting by p-val and frequency
        recurcodons$qval = p.adjust(recurcodons$pval, method="BH", n=numtests) # P-value adjustment for all possible changes
        rownames(recurcodons) = NULL
        
        # Additional annotation
        annotsubs$mutaa = substr(annotsubs$aachange,nchar(annotsubs$aachange),nchar(annotsubs$aachange))
        annotsubs$simplent = paste(annotsubs$ref,annotsubs$mut,sep=">")
        annotsubs$mutnt = paste(annotsubs$chr,annotsubs$pos,annotsubs$simplent,annotsubs$mutaa,sep="_")
        aux = split(annotsubs, f=annotsubs$codonsub)
        recurcodons_ext = recurcodons
        recurcodons_ext$codonsub = paste(recurcodons_ext$chr,recurcodons_ext$gene,recurcodons_ext$codon,sep=":")
        recurcodons_ext$mutnt = recurcodons_ext$mutaa = NA
        for (j in 1:nrow(recurcodons_ext)) {
            x = aux[[recurcodons_ext$codonsub[j]]]
            f = sort(table(x$mutaa),decreasing=T)
            recurcodons_ext$mutaa[j] = paste(names(f),f,sep=":",collapse="|")
            f = sort(table(x$mutnt),decreasing=T)
            recurcodons_ext$mutnt[j] = paste(names(f),f,sep=":",collapse="|")
        }
        
    } else {
        recurcodons = recurcodons_ext = NULL
        warning("No codon was found with the minimum recurrence requested [default min_recurr=2]")
    }
    
    if (is.null(codon_list)) {
        return(list(recurcodons=recurcodons, recurcodons_ext=recurcodons_ext, overdisp=thetaout, LL=LL))
    } else {
        return(list(recurcodons=recurcodons, recurcodons_ext=recurcodons_ext, overdisp=thetaout, LL=LL, globaldnds_knowncodons=globaldnds_knowncodons))
    }
}