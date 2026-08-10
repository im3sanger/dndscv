#' geneciloc
#'
#' Function to calculate confidence intervals for dN/dS values per gene under the dNdSloc model using a Poisson ratio.
#'
#' @author Inigo Martincorena (Wellcome Sanger Institute)
#' @details Martincorena I, et al. (2017) Universal patterns of selection in cancer and somatic tissues. Cell. 171(5):1029-1041.
#' 
#' @param dndsout Output object from dndscv.
#' @param gene_list List of genes to restrict the analysis (by default, all genes in dndsout will be analysed)
#' @param level Confidence level desired [default = 0.95]
#'
#' @return ci: Dataframe with the confidence intervals for dN/dS ratios per gene under the dNdSloc model.
#' 
#' @export

geneciloc = function(dndsout, gene_list = NULL, level = 0.95) {
  
  # Ensuring valid level value
  if (level > 1) {
    warning("Confidence level must be lower than 1, using 0.95 as default")
    level = 0.95
  }
  
  # Subfunction to run the poisson.test function on each gene (r = expected neutral ratio given the substitution model and the gene sequence)
  poissonci = function(ns, s, r) {
    return(t(rbind((ns/s)/r, sapply(1:length(ns), function(j) poisson.test(x=c(ns[j],s[j]), conf.level = level)$conf.int / r[j]))))
  }
  
  ## Calculating CI95% across all genes
  
  if (!is.null(gene_list)) {
    genemuts = dndsout$genemuts[dndsout$genemuts$gene_name %in% gene_list,]
  } else {
    genemuts = dndsout$genemuts
  }
  ci95df = data.frame(gene_name = genemuts$gene_name, mis_mle=NA, tru_mle=NA, non_mle=NA, spl_mle = NA, mis_low = NA, tru_low=NA, non_low = NA, spl_low = NA, mis_high = NA, tru_high=NA, non_high = NA, spl_high = NA)
  
  # Missense
  ci95df[,c("mis_mle","mis_low","mis_high")] = poissonci(genemuts$n_mis, genemuts$n_syn, genemuts$exp_mis/genemuts$exp_syn)
  # Truncating (nonsense + essential splice sites)
  ci95df[,c("tru_mle","tru_low","tru_high")] = poissonci(genemuts$n_non+genemuts$n_spl, genemuts$n_syn, (genemuts$exp_non+genemuts$exp_spl)/genemuts$exp_syn)
  # Nonsense
  ci95df[,c("non_mle","non_low","non_high")] = poissonci(genemuts$n_non, genemuts$n_syn, genemuts$exp_non/genemuts$exp_syn)
  # Essential splice
  ci95df[,c("spl_mle","spl_low","spl_high")] = poissonci(genemuts$n_spl, genemuts$n_syn, genemuts$exp_spl/genemuts$exp_syn)
  
  return(ci95df)
}

  