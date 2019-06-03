#' buildcodon
#' 
#' Function to build a refcodon object from a RefCDS object. This function takes a RefCDS object as input and adds to it two fields required to run the codondnds function. Recommended usage: RefCDS = buildcodon(RefCDS)
#' 
#' @author Inigo Martincorena (Wellcome Sanger Institute)
#' @details Martincorena I, et al. (2017) Universal patterns of selection in cancer and somatic tissues. Cell. 171(5):1029-1041.
#' 
#' @param refcds Input RefCDS object
#' @param numcode NCBI genetic code number (default = 1; standard genetic code). To see the list of genetic codes supported use: ? seqinr::translate
#'
#' @export

buildcodon = function(refcds, numcode = 1) {
    
    ## 1. Valid chromosomes and reference CDS per gene
    message("Adding codon-level information to RefCDS to run codondnds...")

    nt = c("A","C","G","T")
    trinuc_list = paste(rep(nt,each=16,times=1), rep(nt,each=4,times=4), rep(nt,each=1,times=16), sep="")
    trinuc_ind = structure(1:64, names=trinuc_list)
    trinuc_subs = NULL; for (j in 1:length(trinuc_list)) { trinuc_subs = c(trinuc_subs, paste(trinuc_list[j], paste(substr(trinuc_list[j],1,1), setdiff(nt,substr(trinuc_list[j],2,2)), substr(trinuc_list[j],3,3), sep=""), sep=">")) }
    trinuc_subsind = structure(1:192, names=trinuc_subs)
    
    # Precalculating a 64x64 matrix with the functional impact of each codon transition (1=Synonymous, 2=Missense, 3=Nonsense)
    impact_matrix = array(NA, dim=c(64,64))
    colnames(impact_matrix) = rownames(impact_matrix) = trinuc_list
    for (j in 1:64) {
        for (h in 1:64) {
            from_aa = seqinr::translate(strsplit(trinuc_list[j],"")[[1]], numcode = numcode)
            to_aa = seqinr::translate(strsplit(trinuc_list[h],"")[[1]], numcode = numcode)
            # Annotating the impact of the mutation
            if (to_aa == from_aa){ 
                impact_matrix[j,h] = 1
            } else if (to_aa == "*"){
                impact_matrix[j,h] = 3
            } else if ((to_aa != "*") & (from_aa != "*") & (to_aa != from_aa)){
                impact_matrix[j,h] = 2
            } else if (from_aa=="*") {
                impact_matrix[j,h] = NA
            }
        }
    }
        
    # Initialising and populating the Refcodon object
    refcodon = array(list(NULL), length(refcds)) # Initialising empty object

    for (j in 1:length(refcds)) {
        
        cdsseq = as.character(as.vector(refcds[[j]]$seq_cds))
        cdsseq1up = as.character(as.vector(refcds[[j]]$seq_cds1up))
        cdsseq1down = as.character(as.vector(refcds[[j]]$seq_cds1down))
        
        # Exonic mutations
        
        ind = rep(1:length(cdsseq), each=3)
        old_trinuc = paste(cdsseq1up[ind], cdsseq[ind], cdsseq1down[ind], sep="")
        new_base = c(sapply(cdsseq, function(x) nt[nt!=x]))
        new_trinuc = paste(cdsseq1up[ind], new_base, cdsseq1down[ind], sep="")
        codon_start = rep(seq(1,length(cdsseq),by=3),each=9)
        old_codon = paste(cdsseq[codon_start], cdsseq[codon_start+1], cdsseq[codon_start+2], sep="")
        pos_in_codon = rep(rep(1:3, each=3), length.out=length(old_codon))
        aux = strsplit(old_codon,"")
        new_codon = sapply(1:length(old_codon), function(x) { new_codonx = aux[[x]]; new_codonx[pos_in_codon[x]] = new_base[x]; return(new_codonx) } )
        new_codon = paste(new_codon[1,], new_codon[2,], new_codon[3,], sep="")
        
        imp = impact_matrix[(trinuc_ind[new_codon]-1)*64 + trinuc_ind[old_codon]]
        matrind = trinuc_subsind[paste(old_trinuc, new_trinuc, sep=">")]
        
        refcds[[j]]$codon_impact = imp
        refcds[[j]]$codon_rates = matrind
        
        if (round(j/1000)==(j/1000)) { message(sprintf('    %0.3g%% ...', round(j/length(refcds),2)*100)) }
    }
    
    return(refcds)
    
} # EOF
