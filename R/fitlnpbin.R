#' fitlnpbin
#'
#' Function to fit a Lognormal-Poisson model to estimate overdispersion on synonymous changes for sitednds and codondnds.
#'
#' @author Inigo Martincorena (Wellcome Sanger Institute)
#' 
#' @param nvec Vector of observed counts of mutations per site.
#' @param rvec Vector of expected counts of mutations per site.
#' @param level Confidence level desired for the confidence interval of the overdispersion parameter [defaul=0.95]
#' @param theta_option 2 options: "mle" (uses the MLE of the overdispersion parameter) or "conservative" (uses the conservative bound of the CI95). Values other than "mle" will lead to the conservative option [default="conservative"]
#' @param numbins Number of bins to discretise the rvec vector [default=1e4]. This enables fast execution of the LNP model in datasets of arbitrarily any size.
#'
#' @return 'fitlnpbin' returns the maximum likelihood estimate and confidence intervals of the "sig" overdispersion parameter of the LNP model:
#' 
#' @export

fitlnpbin = function(nvec, rvec, level = 0.95, theta_option = "conservative", numbins = 1e4) {
    
    # 1. Binning the r vector
    minrate = 1e-10 # Values <<1/exome_length are due to 0 observed counts for a given trinucleotide
    rvec = pmax(minrate, rvec) # Setting values below minrate to minrate
    br = cut(log(rvec),breaks=numbins) # Binning rvec in log space
    binmeans = tapply(rvec, br, mean) # Mean value per bin
    rvecbinned = as.numeric(binmeans[br]) # Binned values for rvec
    message(sprintf("    Binning the rate vector: maximum deviation of %0.3g", max(abs(rvecbinned-rvec)/rvec)))
    rvec = rvecbinned # Using the binned values
    freqs = as.matrix(plyr::count(cbind(rvec,nvec))) # Frequency table
    
    # Excluding sites with rates < minrate from the calculation (they should yield LL=0)
    freqs = freqs[which(freqs[,1] > minrate), ]
    
    # 2. Vectorising dpoilog
    dpoilog = poilog::dpoilog
    lnp = function(freqs, sig) {
        -sum(apply(freqs, 1, function(x) x[3]*log(dpoilog(n=x[2], mu=log(x[1])-sig^2/2, sig=sig)))) # vectorised dpoilog with fixed expected rates and log-transformed
    }
    
    # 3. Estimating the MLE: grid search followed by optim within reasonable bounds
    lnp_mle = function(minsig=1e-2, maxsig=5, bins=5, iter=8) {
        
        lls = sigs = NULL # Saving the log-likelihoods calculated
        
        # 1. Grid search to identify reasonable bounds
        for (j in 1:iter) {
            if (j==1) {
                sigvec = sort(exp(seq(log(minsig), log(maxsig), length.out=bins))) # Initial vals
            } else {
                sigvec = seq(sigvec[pmax(1,ind-1)], sigvec[pmin(bins,ind+1)], length.out=bins) # Refining previous iteration
            }
            proflik = sapply(sigvec, function(sig) lnp(freqs=freqs, sig=sig))
            ind = which.min(proflik)
            
            sigs = c(sigs, sigvec) # Saving the result
            lls = c(lls, proflik) # Saving the result
        }
        
        # 2. Optim for precise estimation of the MLE (optim without narrow bounds tends to fail)
        f = function(sig, n=nvec, r=rvec) { lnp(freqs=freqs, sig=sig) }
        ml = optimize(f, interval=c(sigvec[pmax(1,ind-1)], sigvec[pmin(bins,ind+1)]))
        
        sigs = c(sigs, ml$minimum) # Saving the result
        lls = c(lls, ml$objective) # Saving the result
        ll = cbind(sigs,lls)
        
        return(list(ml=ml, ll=unique(ll[order(ll[,1]),])))
    }
    
    ml = lnp_mle(minsig=1e-2, maxsig=5, bins=5, iter=8) # Maximum likelihood estimate of the overdispersion
    
    
    # 4. Estimating the lower bound of the CI95% using profile likelihood
    #    This is done exploiting the points already evaluated for the MLE
    
    if (theta_option == "mle") {
        ml$sig_ci95 = NA # We only estimate the lower bound of sig if requested by the user
    } else { 
        grid_proflik = function(minsig=1e-2, maxsig=5, bins=5, iter=8) {
            for (j in 1:iter) {
                if (j==1) {
                    sigvec = ml$ll[,1]
                    ind = min(which(ml$ll[,1]>ml$ml$minimum & (ml$ll[,2]-ml$ml$objective)>qchisq(.95,1)/2)) # First value outside of bounds
                }
                sigvec = seq(sigvec[pmax(1,ind-1)], sigvec[pmin(length(sigvec),ind)], length.out=bins) # New grid based on the previous iteration
                proflik = sapply(sigvec, function(sig) lnp(freqs=freqs, sig=sig)) # Calculating log-likelihoods
                ind = min(which(sigvec>ml$ml$minimum & (proflik-ml$ml$objective)>qchisq(.95,1)/2)) # First value outside of bounds
            }
            return(sigvec[ind]) # Conservative estimate for the lower bound of the CI95% for sig
        }
        ml$sig_ci95 = grid_proflik(minsig=1e-2, maxsig=5, bins=5, iter=8)
    }
    
    return(ml)
}
