#' geneci
#'
#' Function to calculate confidence intervals for dN/dS values per gene under the dNdScv model using profile likelihood. To generate a valid dndsout input object for this function, use outmats=T when running dndscv.
#'
#' @author Inigo Martincorena (Wellcome Sanger Institute)
#' @details Martincorena I, et al. (2017) Universal patterns of selection in cancer and somatic tissues. Cell. 171(5):1029-1041.
#' 
#' @param dndsout Output object from dndscv.
#' @param gene_list List of genes to restrict the analysis (by default, all genes in dndsout will be analysed)
#' @param level Confidence level desired [default = 0.95]
#'
#' @return 'dndscv' returns a list of objects:
#' @return - ci: Confidence intervals for the dN/dS ratios per gene under the dNdScv model.
#' 
#' @export

geneci = function(dndsout, gene_list = NULL, level = 0.95) {

    # Ensuring valid level value
    if (level > 1) {
        warning("Confidence level must be lower than 1, using 0.95 as default")
        level = 0.95
    }
    
    # N and L matrices
    N = dndsout$N
    L = dndsout$L
    if (length(N)==0) { stop(sprintf("Invalid input: the dndsout input object must be generated using outmats=T as an argument to dndscv.")) }
    if (nrow(dndsout$mle_submodel)!=195) { stop(sprintf("Invalid input: dndsout must be generated using the default trinucleotide substitution model in dndscv."))}
    
    # Restricting the analysis to an input list of genes
    if (!is.null(gene_list)) {
        g = as.vector(dndsout$genemuts$gene_name) # Genes in the input object
        nonex = gene_list[!(gene_list %in% g)] # Excluding genes from the input gene_list if they are not present in the input dndsout object
        if (length(nonex)>0) {
            warning(sprintf("The following input gene names are not in dndsout input object and will not be analysed: %s.", paste(nonex,collapse=", ")))
        }
        dndsout$annotmuts = dndsout$annotmuts[which(dndsout$annotmuts$gene %in% gene_list), ] # Restricting to genes of interest
        dndsout$genemuts = dndsout$genemuts[which(g %in% gene_list), ] # Restricting to genes of interest
        N = N[,,which(g %in% gene_list)] # Restricting to genes of interest
        L = L[,,which(g %in% gene_list)] # Restricting to genes of interest
    }
    gene_list = as.vector(dndsout$genemuts$gene_name)
    
    wnoneqspl = all(dndsout$sel_cv$wnon_cv==dndsout$sel_cv$wspl_cv) # Deciding on wnon==wspl based on the input object
    
    ## Subfunction: Analytical opt_t (aka tML) given fixed w values
    mle_tcvgivenw = function(n, theta, exp_neutral_cv, E) {
        shape = theta; scale = exp_neutral_cv/theta
        tml = (n+shape-1)/(1+E+(1/scale))
        if (shape<=1) { # i.e. when theta<=1
            tml = max(shape*scale,tml) # i.e. tml is bounded to the mean of the gamma (i.e. y[9]) when theta<=1
        }
        return(pmax(tml,1e-6))
    }
    
    ## Subfunction: Log-Likelihood of the model given fixed w values (requires finding MLEs for t and the free w values given the fixed w values)
    loglik_givenw = function(w,x,y,mutrates,theta,wtype,wnoneqspl) {
    
        # 1. tML given w
        exp_neutral_cv = y[9]
        exp_rel = y[6:8]/y[5]
        n = y[1] + sum(y[wtype+1])
        E = sum(exp_rel[wtype])*w
        tML = mle_tcvgivenw(n, theta, exp_neutral_cv, E)
        mrfold = max(1e-10, tML/y[5]) # Correction factor of "t" under the model
    
        # 2. Calculating the MLEs of the unconstrained w values
        if (!wnoneqspl) {
          wfree = y[2:4]/y[6:8]/mrfold; wfree[y[2:4]==0] = 0 # MLEs for w given tval
        } else {
          wmisfree = y[2]/y[6]/mrfold; wmisfree[y[2]==0] = 0
          wtruncfree = sum(y[3:4])/sum(y[7:8])/mrfold; wtruncfree[sum(y[3:4])==0] = 0
          wfree = c(wmisfree,wtruncfree,wtruncfree) # MLEs for w given tval
        }
        wfree[wtype] = w # Replacing free w values by fixed input values
    
        # 2. loglik of the model under tML and w
        llpois = sum(dpois(x=x$n, lambda=x$l*mutrates*mrfold*t(array(c(1,wfree),dim=c(4,length(mutrates)))), log=T))
        llgamm = dgamma(x=tML, shape=theta, scale=exp_neutral_cv/theta, log=T)
        return(-(llpois+llgamm))
    }
    
    ## Subfunction: Working with vector inputs
    loglik_vec = function(wfixed,x,y,mutrates,theta,wtype,wnoneqspl) {
        sapply(wfixed, function(w) loglik_givenw(w,x,y,mutrates,theta,wtype,wnoneqspl))
    }
    
    
    ## Subfunction: iterative search for the CI95% boundaries for wvec
    iterative_search_ci95 = function(wtype,x,y,mutrates,theta,wmle,ml,grid_size=10,iter=10,wnoneqspl=T,wmax = 10000) {
    
      if (wmle[wtype][1]<wmax) {
    
        # Iteratively searching for the lower bound of the CI95% for "t"
        if (wmle[wtype][1]>0) {
          search_range = c(1e-9, wmle[wtype][1])
          for (it in 1:iter) {
            wvec = seq(search_range[1], search_range[2],length.out=grid_size)
            ll = -loglik_vec(wvec,x,y,mutrates,theta,wtype,wnoneqspl)
            lr = 2*(ml-ll) > qchisq(p=level,df=1)
            ind = max(which(wvec<=wmle[wtype][1] & lr))
            search_range = c(wvec[ind], wvec[ind+1])
          }
          w_low = wvec[ind]
        } else {
          w_low = 0
        }
        
        # Iteratively searching for the higher bound of the CI95% for "t"
        search_range = c(wmle[wtype][1], wmax)
        llhighbound = -loglik_vec(wmax,x,y,mutrates,theta,wtype,wnoneqspl)
        outofboundaries = !(2*(ml-llhighbound) > qchisq(p=level,df=1))
        if (!outofboundaries) {
          for (it in 1:iter) {
            wvec = seq(search_range[1], search_range[2],length.out=grid_size)
            ll = -loglik_vec(wvec,x,y,mutrates,theta,wtype,wnoneqspl)
            lr = 2*(ml-ll) > qchisq(p=level,df=1)
            ind = min(which(wvec>=wmle[wtype][1] & lr))
            search_range = c(wvec[ind-1], wvec[ind])
          }
          w_high = wvec[ind]
        } else {
          w_high = wmax
        }
    
      } else {
        wmle[wtype] = w_low = w_high = wmax
      }
    
      return(c(wmle[wtype][1],w_low,w_high))
    }
    
    
    ## Subfunction: calculate the MLEs and CI95% of each independent w value (unconstraining the other values)
    ci95cv_intt = function(x,y,mutrates,theta,grid_size=10,iter=10,wnoneqspl=T) {
    
      # MLE
      exp_neutral_cv = y[9]
      n = y[1]; E = 0 # Only synonymous mutations are considered
      tML = mle_tcvgivenw(n, theta, exp_neutral_cv, E)
      mrfold = max(1e-10, tML/y[5])
      if (!wnoneqspl) {
        wmle = y[2:4]/y[6:8]/mrfold; wmle[y[2:4]==0] = 0 # MLEs for w given tval
      } else {
        wmisfree = y[2]/y[6]/mrfold; wmisfree[y[2]==0] = 0
        wtruncfree = sum(y[3:4])/sum(y[7:8])/mrfold; wtruncfree[sum(y[3:4])==0] = 0
        wmle = c(wmisfree,wtruncfree,wtruncfree) # MLEs for w given tval
      }
      llpois = sum(dpois(x=x$n, lambda=x$l*mutrates*mrfold*t(array(c(1,wmle),dim=c(4,length(mutrates)))), log=T))
      llgamm = dgamma(x=tML, shape=theta, scale=y[9]/theta, log=T)
      ml = llpois+llgamm
    
      # Iteratively searching for the lower bound of the CI95% for "t"
      w_ci95 = array(NA,c(3,3))
      colnames(w_ci95) = c("mle","low","high")
      rownames(w_ci95) = c("mis","non","spl")
      if (!wnoneqspl) {
        for (h in 1:3) {
          w_ci95[h,] = iterative_search_ci95(wtype=h,x,y,mutrates,theta,wmle,ml,grid_size,iter,wnoneqspl)
        }
      } else {
        w_ci95[1,] = iterative_search_ci95(wtype=1,x,y,mutrates,theta,wmle,ml,grid_size,iter,wnoneqspl)
        w_ci95[2,] = iterative_search_ci95(wtype=c(2,3),x,y,mutrates,theta,wmle,ml,grid_size,iter,wnoneqspl)
        w_ci95[3,] = w_ci95[2,]
      }
      return(w_ci95)
    }
    
    
    ## Calculating CI95% across all genes
    
    message("Calculating CI95 across all genes...")
    
    ci95 = array(NA, dim=c(length(gene_list),9))
    colnames(ci95) = c("mis_mle","non_mle","spl_mle","mis_low","non_low","spl_low","mis_high","non_high","spl_high")
    theta = dndsout$nbreg$theta
    
    data("submod_192r_3w", package="dndscv")
    parmle =  setNames(dndsout$mle_submodel[,2], dndsout$mle_submodel[,1])
    mutrates = sapply(substmodel[,1], function(x) prod(parmle[base::strsplit(x,split="\\*")[[1]]])) # Expected rate per available site
    
    for (j in 1:length(gene_list)) {
        geneind = which(dndsout$genemuts$gene_name==gene_list[j])
        y = as.numeric(dndsout$genemuts[geneind,-1])
        if (length(gene_list)==1) {
            x = list(n=N, l=L)
        } else {
            x = list(n=N[,,geneind], l=L[,,geneind])
        }
        ci95[j,] = c(ci95cv_intt(x,y,mutrates,theta,grid_size=10,iter=10,wnoneqspl=wnoneqspl))
        if (round(j/1000)==(j/1000)) { print(j/length(gene_list), digits=2) } # Progress
    }
    
    ci95df = cbind(gene=gene_list, as.data.frame(ci95))
    
    # Restricting columns if we forced wnon==wspl
    if (wnoneqspl==T) {
        ci95df = ci95df[,-c(4,7,10)]
        colnames(ci95df) = c("gene","mis_mle","tru_mle","mis_low","tru_low","mis_high","tru_high")
    }

    return(ci95df)
    
} # EOF
