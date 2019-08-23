#' sitednds
#'
#' Function to estimate site-wise dN/dS values and p-values against neutrality. To generate a valid input object for this function, use outmats=T when running dndscv. This function is in testing, please interpret the results with caution. Also note that recurrent artefacts or SNP contamination can violate the null model and dominate the list of sites under apparent selection. A considerable number of significant synonymous sites may reflect a problem with the data. Be very critical of the results and if suspicious sites appear recurrently mutated consider refining the variant calling (e.g. using a better unmatched normal panel). In the future, this function may be extended to perform inferences at a codon level instead of at a single-base level.
#'
#' @author Inigo Martincorena (Wellcome Sanger Institute)
#' 
#' @param dndsout Output object from dndscv. To generate a valid input object for this function, use outmats=T when running dndscv.
#' @param min_recurr Minimum number of mutations per site to estimate site-wise dN/dS ratios. [default=2]
#' @param gene_list List of genes to restrict the p-value and q-value calculations (Restricted Hypothesis Testing). Note that q-values are only valid if the list of genes is decided a priori. [default=NULL, sitednds will be run on all genes in dndsout]
#' @param site_list List of hotspot sites to restrict the p-value and q-value calculations (Restricted Hypothesis Testing). Note that q-values are only valid if the list of sites is decided a priori. [default=NULL, sitednds will be run on all genes in dndsout]
#' @param trinuc_list List of trinucleotide substitution to restrict the analysis of sitednds. This is used to estimate separate overdispersion parameters for different substitution contexts [default=NULL, sitednds will be run on all substitution contexts]
#' @param theta_option 2 options: "mle" (uses the MLE of the overdispersion parameter) or "conservative" (uses the conservative bound of the CI95). Values other than "mle" will lead to the conservative option [default="conservative"]
#' @param syn_drivers Vector with a list of known synonymous driver mutations to exclude from the background model [default="TP53:T125T"]. See Martincorena et al., Cell, 2017 (PMID:29056346).
#' @param method Overdispersion model: NB = Negative Binomial (Gamma-Poisson), LNP = Poisson-Lognormal (see Hess et al., BiorXiv, 2019). [default="NB"]
#' @param numbins Number of bins to discretise the rvec vector [default=1e4]. This enables fast execution of the LNP model in datasets of arbitrarily any size.
#' @param kc List of a-priori known cancer genes (to be excluded when fitting the background model)
#'
#' @return 'sitednds' returns a table of recurrently mutated sites and the estimates of the size parameter:
#' @return - recursites: Table of recurrently mutated sites with site-wise dN/dS values and p-values
#' @return - overdisp: Maximum likelihood estimate and CI95% for the overdispersion parameter (the size parameter of the negative binomial distribution or the sigma parameter of the lognormal distribution). The lower the size value or the higher the sigma value the higher the variation of the mutation rate across sites not captured by the trinucleotide change or by variation across genes.
#' @return - fpr_nonsyn_q05: Fraction of the significant non-synonymous sites (qval<0.05) that are estimated to be false positives. This assumes that all synonymous mutations (except those in TP53 and CDKN2A) are false positives, thus offering a conservative estimate of the false positive rate.
#' @return - LL: Log-likelihood of the fit of the overdispersed model (see "method" argument) to all synonymous sites.
#'
#' @export

sitednds = function(dndsout, min_recurr = 2, gene_list = NULL, site_list = NULL, trinuc_list = NULL, theta_option = "conservative", syn_drivers = "TP53:T125T", method = "NB", numbins = 1e4, kc = "cgc81") {
    
    ## 1. Fitting a negative binomial distribution at the site level considering the background mutation rate of the gene and of each trinucleotide
    message("[1] Site-wise overdispersed model accounting for trinucleotides and relative gene mutability...")
    
    # N and L matrices for synonymous mutations
    if (length(dndsout$N)==0) { stop(sprintf("Invalid input: dndsout must be generated using outmats=T in dndscv.")) }
    if (nrow(dndsout$mle_submodel)!=195) { stop("Invalid input: dndsout must be generated using the default trinucleotide substitution model in dndscv.") }
    
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
        numtests = sum(dndsout$L[,,which(g %in% gene_list)])
    } else {
        numtests = sum(dndsout$L)
    }
    
    # Input: known cancer genes to exclude from the background model fitting (the user can input a gene list as a character vector)
    if (is.null(kc)) {
        known_cancergenes = ""
    } else if (kc[1] %in% c("cgc81")) {
        data(list=sprintf("cancergenes_%s",kc), package="dndscv")
    } else {
        known_cancergenes = kc
    }
    
    # L matrix
    L = dndsout$L
    L[,2:4,which(as.vector(dndsout$genemuts$gene_name) %in% known_cancergenes)] = 0 # Removing non-synonymous sites in known_cancergenes
    L = apply(L, c(1,3), sum) # Total number of sites to be considered in the background model
    
    
    # Counts of observed mutations
    
    annotsubs = dndsout$annotmuts[which(dndsout$annotmuts$impact %in% c("Synonymous","Missense","Nonsense","Essential_Splice")),]
    num_syn_drivers_masked = sum(paste(annotsubs$gene,annotsubs$aachange,sep=":") %in% syn_drivers)
    
    if (!is.null(trinuc_list)) { # Restricting sitednds to certain trinucleotide changes
        annotsubs = annotsubs[which(paste(annotsubs$ref3_cod,annotsubs$mut3_cod,sep=">") %in% trinuc_list), ]
        if (nrow(annotsubs)==0) { stop("No mutations left after restricting by trinuc_list. Please review your input arguments and your mutation table (dndsout$annotmuts).") }
    }
    
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
    
    if (!is.null(trinuc_list)) { # Restricting sitednds to certain trinucleotide changes 
        sm[!(names(sm) %in% trinuc_list)] = 0 # Setting all other rates to zero
    }
    
    mat_trisub = array(sm, dim=c(192,nrow(dndsout$genemuts))) # Relative mutation rates by trinucleotide
    mat_relmr = t(array(relmr, dim=c(nrow(dndsout$genemuts),192))) # Relative mutation rates by gene
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
    
    # Mutations for the background model (excluding non-synonymous mutations in known_cancergenes)
    synsites = mutsites[!(mutsites$impact!="Synonymous" & mutsites$gene %in% known_cancergenes),]
    synsites = synsites[!(paste(synsites$gene,synsites$aachange,sep=":") %in% syn_drivers),]
    Lcum = array(cumsum(L), dim=dim(L)) # Cumulative L indicating the position to place a given mutation in the nvec following rvec
    synsites$vecindex = apply(as.matrix(synsites[,c("trindex","geneindex")]), 1, function(x) Lcum[x[1], x[2]]) # Index for the mutation
    synsites = synsites[order(synsites$vecindex), ] # Sorting by index in the nvec

    # Stop execution with an error if there are no synonymous mutations
    if (nrow(synsites)<2) {
        stop("Too few synonymous mutations found in the input. sitednds cannot run without synonymous mutations.")
    }
    
    # Correcting the index when there are multiple synonymous mutations in the same gene and trinucleotide class
    s = snew = synsites$vecindex
    sameind = 0
    for (j in 2:nrow(synsites)) {
        if (s[j]<=s[j-1]) {
            sameind = sameind + 1 # Annotating a run of elements
            snew[j] = s[j-1] - sameind # We assign it an earlier position in the vector
        } else {
            sameind = 0
        }
    }
    synsites$vecindex2 = snew
    
    nvec[synsites$vecindex2] = synsites$freq # Expanded nvec for the negative binomial regression
    rvec = rvec * (sum(dndsout$genemuts$n_syn)-num_syn_drivers_masked) / sum(dndsout$genemuts$exp_syn_cv) # Minor correction ensuring that global observed and expected rates are identical (this works after subsetting substitutions)
    
    # Estimation of overdispersion: Using optimize appears to yield reliable results. Problems experienced with fitdistr, glm.nb and theta.ml. Consider using grid search if problems appear with optimize.
    if (method=="LNP") { # Modelling rates per site with a Poisson-Lognormal mixture
        
        lnp_est = fitlnpbin(nvec, rvec, theta_option = theta_option, numbins = numbins)
        theta_ml = lnp_est$ml$minimum
        theta_ci95 = lnp_est$sig_ci95
        LL = -lnp_est$ml$objective # LogLik
        thetaout = setNames(c(theta_ml, theta_ci95), c("MLE","CI95_high"))
        
    } else { # Modelling rates per site as negative binomially distributed (i.e. quantifying uncertainty above Poisson using a Gamma)
        
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
    
    
    ## 2. Calculating site-wise dN/dS ratios and P-values for recurrently mutated sites
    message("[2] Calculating site-wise dN/dS ratios and p-values...")
    
    # Theta option
    if (theta_option=="mle" | theta_option=="MLE") {
        theta = theta_ml
    } else { # Conservative
        message("    Using the conservative bound of the confidence interval of the overdispersion parameter.")
        theta = theta_ci95[1]
    }
    
    # Creating the recursites object
    recursites = mutsites[, c("chr","pos","ref","mut","gene","aachange","impact","ref3_cod","mut3_cod","freq")]
    recursites$mu = relmr[recursites$gene] * sm[paste(recursites$ref3_cod,recursites$mut3_cod,sep=">")]
    
    # Gene RHT
    if (!is.null(gene_list)) {
        message("    Peforming Restricted Hypothesis Testing on the input list of a-priori genes")
        recursites = recursites[which(recursites$gene %in% gene_list), ] # Restricting the p-value and q-value calculations to gene_list
    }
    
    # Site RHT
    if (!is.null(site_list)) {
        message("    Peforming Restricted Hypothesis Testing on the input list of a-priori sites (numtests = length(site_list))")
        mutstr = paste(recursites$chr,recursites$pos,recursites$ref,recursites$mut,recursites$gene,recursites$aachange,recursites$ref3_cod,recursites$mut3_cod,sep=":")
        if (!any(mutstr %in% site_list)) {
            stop("No mutation was observed in the restricted list of known hotspots. Site-RHT cannot be run.")
        }
        recursites = recursites[which(mutstr %in% site_list), ] # Restricting the p-value and q-value calculations to site_list
        numtests = length(site_list)
        
        # Calculating global dN/dS ratios at known hotspots
        auxsites = as.data.frame(do.call("rbind",strsplit(site_list,split=":")), stringsAsFactors=F)
        auxsites = auxsites[auxsites$V5 %in% names(relmr), ]
        neutralexp = sum(relmr[auxsites$V5]*sm[paste(auxsites$V7,auxsites$V8,sep=">")]) # Number of mutations expected at known hotspots expected under neutrality
        numobs = sum(recursites$freq) # Number observed
        poistest = poisson.test(numobs, T=neutralexp)
        globaldnds_knownsites = setNames(c(numobs, neutralexp, poistest$estimate, poistest$conf.int), c("obs","exp","dnds","cilow","cihigh"))
        message(sprintf("    Mutations at known hotspots: %0.0f observed, %0.3g expected, obs/exp~%0.3g (CI95:%0.3g,%0.3g).", globaldnds_knownsites[1], globaldnds_knownsites[2], globaldnds_knownsites[3], globaldnds_knownsites[4], globaldnds_knownsites[5]))
    }
    
    # Restricting the recursites output by min_recurr
    recursites = recursites[recursites$freq>=min_recurr, ] # Restricting the output to sites with min_recurr
    
    if (nrow(recursites)>0) {
        
        recursites$dnds = recursites$freq / recursites$mu # Site-wise dN/dS (point estimate)
        
        if (method=="LNP") { # Modelling rates per site with a Poisson-Lognormal mixture
            
            # Cumulative Lognormal-Poisson using poilog::dpoilog
            dpoilog = poilog::dpoilog
            ppoilog = function(n, mu, sig) {
                p = sum(dpoilog(n=floor(n+1):floor(n*10+1000), mu=log(mu)-sig^2/2, sig=sig))
                return(p)
            }
            
            message(sprintf("    Modelling substitution rates using a Lognormal-Poisson: sig = %0.3g (upperbound = %0.3g)", theta_ml, theta_ci95))
            recursites$pval = apply(recursites, 1, function(x) ppoilog(n=as.numeric(x["freq"])-0.5, mu=as.numeric(x["mu"]), sig=theta))

        } else { # Negative binomial model
            
            message(sprintf("    Modelling substitution rates using a Negative Binomial: theta = %0.3g (CI95:%0.3g,%0.3g)", theta_ml, theta_ci95[1], theta_ci95[2]))
            recursites$pval = pnbinom(q=recursites$freq-0.5, mu=recursites$mu, size=theta, lower.tail=F)
        }
        
        recursites = recursites[order(recursites$pval, -recursites$freq), ] # Sorting by p-val and frequency
        recursites$qval = p.adjust(recursites$pval, method="BH", n=numtests) # P-value adjustment for all possible changes
        rownames(recursites) = NULL
        
        # Estimating False Positive Rates based on the observed number of significant synonymous hits
        
        qcutoff = 0.05 # q-value cutoff to estimate false positive rates
        if (any(recursites$qval<qcutoff)) {
            
            signifsites = recursites[recursites$qval<qcutoff, ]
            obs_hits = c(sum(signifsites$impact=="Synonymous" & !(signifsites$gene %in% c("TP53","CDKN2A","CDKN2A.p14arf"))), sum(signifsites$impact!="Synonymous"))
            exp_frac = c(sum(dndsout$genemuts$exp_syn), sum(dndsout$genemuts$exp_mis+dndsout$genemuts$exp_non+dndsout$genemuts$exp_spl))
            obsexp = (obs_hits[2]/obs_hits[1])/(exp_frac[2]/exp_frac[1])
            fpr_nonsyn = list()
            fpr_nonsyn$estimate = 1 / pmax(1,obsexp) # i.e. 1 - driver_fraction
            fpr_nonsyn$conf.int = rev(1 / pmax(1,as.vector(poisson.test(x=obs_hits[2:1], T=exp_frac[2:1])$conf.int)))
            
            if (fpr_nonsyn$estimate>qcutoff) {
                warning(sprintf("The estimated false positive rate for nonsynonymous hits (qval<0.05) is %0.3f (CI95:%0.3f,%0.3f). High false positive rates (>>0.05) evidence problems with the data or the model and mean that the results are not reliable.", fpr_nonsyn$estimate, fpr_nonsyn$conf.int[1], fpr_nonsyn$conf.int[2]))
            }
            
            if (!is.null(site_list)) { # We do not compute FPRs when restricting the analysis to a list of a priori hotspots
                fpr_nonsyn = NULL
            }
            
        } else {
            fpr_nonsyn = NULL
        }
        
    } else {
        recursites = fpr_nonsyn = lnp_est = NULL
        warning("No site was found with the minimum recurrence requested [default min_recurr=2]")
    }
    
    if (is.null(site_list)) {
        return(list(recursites=recursites, overdisp=thetaout, fpr_nonsyn_q05=fpr_nonsyn, LL=LL))
    } else {
        return(list(recursites=recursites, overdisp=thetaout, fpr_nonsyn_q05=fpr_nonsyn, LL=LL, globaldnds_knownsites=globaldnds_knownsites))
    }
}
