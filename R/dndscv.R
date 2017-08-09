#' dNdScv
#'
#' Analyses of selection using the dNdScv and dNdSloc models. Default parameters typically increase the performance of the method on cancer genomic studies. Reference files are currently only available for the GRCh37/hg19 version of the human genome.
#'
#' @author Inigo Martincorena (Wellcome Trust Sanger Institute)
#' @details Martincorena I, et al. (2017) Universal patterns of selection in cancer and somatic tissues. Under revision. Preprint available in BioRxiv: https://doi.org/10.1101/132324
#' 
#' @param mutations Table of mutations (5 columns: sampleID, chr, pos, ref, alt). Only list independent events as mutations.
#' @param gene_list List of genes to restrict the analysis (use for targeted sequencing studies)
#' @param refdb Reference database (path to .rda file)
#' @param sm Substitution model (precomputed models are available in the data directory)
#' @param kc List of a-priori known cancer genes (to be excluded from the indel background model)
#' @param cv Covariates (a matrix of covariates -columns- for each gene -rows-) [default: reference covariates] [cv=NULL runs dndscv without covariates]
#' @param max_muts_per_gene_per_sample If n<Inf, arbitrarily the first n mutations by chr position will be kept
#' @param max_coding_muts_per_sample Hypermutator samples often reduce power to detect selection
#' @param use_indel_sites Use unique indel sites instead of the total number of indels (it tends to be more robust)
#' @param min_indels Minimum number of indels required to run the indel recurrence module
#' @param maxcovs Maximum number of covariates that will be considered (additional columns in the matrix of covariates will be excluded)
#' @param constrain_wnon_wspl This constrains wnon==wspl (this typically leads to higher power to detect selection)
#' @param outp Output: 1 = Global dN/dS values; 2 = Global dN/dS and dNdSloc; 3 = Global dN/dS, dNdSloc and dNdScv
#'
#' @return 'dndscv' returns a list of objects:
#' @return - globaldnds: Global dN/dS estimates across all genes.
#' @return - sel_cv: Gene-wise selection results using dNdScv.
#' @return - sel_loc: Gene-wise selection results using dNdSloc.
#' @return - annotmuts: Annotated coding mutations.
#' @return - genemuts: Observed and expected numbers of mutations per gene.
#' @return - mle_submodel: MLEs of the substitution model.
#' @return - exclsamples: Samples excluded from the analysis.
#' @return - exclmuts: Coding mutations excluded from the analysis.
#' @return - nbreg: Negative binomial regression model for substitutions.
#' @return - nbregind: Negative binomial regression model for indels.
#' @return - poissmodel: Poisson regression model used to fit the substitution model and the global dNdS values.
#' 
#' @export

dndscv = function(mutations, gene_list = NULL, refdb = "hg19", sm = "192r_3w", kc = "cgc81", cv = "hg19", max_muts_per_gene_per_sample = 3, max_coding_muts_per_sample = 3000, use_indel_sites = T, min_indels = 5, maxcovs = 20, constrain_wnon_wspl = T, outp = 3) {

    ## 1. Environment
    message("[1] Loading the environment...")

    mutations[,c(1,2,4,5)] = lapply(mutations[,c(1,2,4,5)], as.character) # Factors to character
    
    # [Input] Reference database
    if (refdb == "hg19") {
        data("refcds_hg19", package="dndscv")
    } else {
        load(refdb)
    }
    
    # [Input] Gene list (The user can input a gene list as a character vector)
    if (is.null(gene_list)) {
        gene_list = sapply(RefCDS, function(x) x$gene_name) # All genes [default]
    } else { # Using only genes in the input gene list
        allg = sapply(RefCDS,function(x) x$gene_name)
        nonex = gene_list[!(gene_list %in% allg)]
        if (length(nonex)>0) { stop(sprintf("The following input gene names are not in the RefCDS database: %s", paste(nonex,collapse=", "))) }
        RefCDS = RefCDS[allg %in% gene_list] # Only input genes
        gr_genes = gr_genes[gr_genes$names %in% gene_list] # Only input genes
    }

    # [Input] Covariates (The user can input a custom set of covariates as a matrix)
    if (is.character(cv)) {
        data(list=sprintf("covariates_%s",cv), package="dndscv")
    } else {
        covs = cv
    }
    
    # [Input] Known cancer genes (The user can input a gene list as a character vector)
    if (kc[1] %in% c("cgc81")) {
        data(list=sprintf("cancergenes_%s",kc), package="dndscv")
    } else {
        known_cancergenes = kc
    }
    
    # [Input] Substitution model (The user can also input a custom substitution model as a matrix)
    if (length(sm)==1) {
        data(list=sprintf("submod_%s",sm), package="dndscv")
    } else {
        substmodel = sm
    }
    
    # Expanding the reference sequences [for faster access]
    for (j in 1:length(RefCDS)) {
        RefCDS[[j]]$seq_cds = base::strsplit(as.character(RefCDS[[j]]$seq_cds), split="")[[1]]
        RefCDS[[j]]$seq_cds1up = base::strsplit(as.character(RefCDS[[j]]$seq_cds1up), split="")[[1]]
        RefCDS[[j]]$seq_cds1down = base::strsplit(as.character(RefCDS[[j]]$seq_cds1down), split="")[[1]]
        if (!is.null(RefCDS[[j]]$seq_splice)) {
            RefCDS[[j]]$seq_splice = base::strsplit(as.character(RefCDS[[j]]$seq_splice), split="")[[1]]
            RefCDS[[j]]$seq_splice1up = base::strsplit(as.character(RefCDS[[j]]$seq_splice1up), split="")[[1]]
            RefCDS[[j]]$seq_splice1down = base::strsplit(as.character(RefCDS[[j]]$seq_splice1down), split="")[[1]]
        }
    }
    
    
    ## 2. Mutation annotation
    message("[2] Annotating the mutations...")
    
    colnames(mutations) = c("sampleID","chr","pos","ref","mut")
    nt = c("A","C","G","T")
    trinucs = paste(rep(nt,each=16,times=1),rep(nt,each=4,times=4),rep(nt,each=1,times=16), sep="")
    trinucinds = setNames(1:64, trinucs)
    
    trinucsubs = NULL
    for (j in 1:length(trinucs)) {
        trinucsubs = c(trinucsubs, paste(trinucs[j], paste(substr(trinucs[j],1,1), setdiff(nt,substr(trinucs[j],2,2)), substr(trinucs[j],3,3), sep=""), sep=">"))
    }
    trinucsubsind = setNames(1:192, trinucsubs)
    
    ind = setNames(1:length(RefCDS), sapply(RefCDS,function(x) x$gene_name))
    gr_genes_ind = ind[gr_genes$names]
    
    # Warning about possible unannotated dinucleotide substitutions
    if (any(diff(mutations$pos)==1)) {
        warning("Mutations observed in contiguous sites within a sample. Please annotate or remove dinucleotide or complex substitutions for best results.")
    }
    
    # Warning about multiple instances of the same mutation in different sampleIDs
    if (nrow(unique(mutations[,2:5])) < nrow(mutations)) {
        warning("Same mutations observed in different sampleIDs. Please verify that these are independent events and remove duplicates otherwise.")
    }
    
    # Mapping mutations to genes
    gr_muts = GRanges(mutations$chr, IRanges(mutations$pos,mutations$pos))
    ol = as.matrix(findOverlaps(gr_muts, gr_genes, type="any", select="all"))
    mutations = mutations[ol[,1],] # Duplicating subs if they hit more than one gene
    mutations$geneind = gr_genes_ind[ol[,2]]
    mutations$gene = sapply(RefCDS,function(x) x$gene_name)[mutations$geneind]
    
    # Optional: Excluding samples exceeding the limit of mutations/sample [see Default parameters]
    nsampl = sort(table(mutations$sampleID))
    exclsamples = NULL
    if (any(nsampl>max_coding_muts_per_sample)) {
        message(sprintf('    Note: %0.0f samples excluded for exceeding the limit of mutations per sample',sum(nsampl>max_coding_muts_per_sample)))
        exclsamples = names(nsampl[nsampl>max_coding_muts_per_sample])
        mutations = mutations[!(mutations$sampleID %in% names(nsampl[nsampl>max_coding_muts_per_sample])),]
    }
    
    # Optional: Limiting the number of mutations per gene per sample (to minimise the impact of unannotated kataegis and other mutation clusters) [see Default parameters]
    mutrank = ave(mutations$pos, paste(mutations$sampleID,mutations$gene), FUN = function(x) rank(x))
    exclmuts = NULL
    if (any(mutrank>max_muts_per_gene_per_sample)) {
        message(sprintf('    Note: %0.0f mutations removed for exceeding the limit of mutations per gene per sample',sum(mutrank>max_muts_per_gene_per_sample)))
        exclmuts = mutations[mutrank>max_muts_per_gene_per_sample,]
        mutations = mutations[mutrank<=max_muts_per_gene_per_sample,]
    }
    
    # Additional annotation of substitutions
    
    snv = (mutations$ref %in% nt & mutations$mut %in% nt)
    indels = mutations[!snv,]
    mutations = mutations[snv,]
    mutations$ref_cod = mutations$ref
    mutations$mut_cod = mutations$mut
    compnt = setNames(rev(nt), nt)
    
    mutations$strand = sapply(RefCDS,function(x) x$strand)[mutations$geneind]
    isminus = (mutations$strand==-1)
    mutations$ref_cod[isminus] = compnt[mutations$ref[isminus]]
    mutations$mut_cod[isminus] = compnt[mutations$mut[isminus]]
    
    for (j in 1:length(RefCDS)) {
        RefCDS[[j]]$N = array(0, dim=c(192,4)) # Initialising the N matrices
    }
    
    # Subfunction: obtaining the codon positions of a coding mutation given the exon intervals
    
    chr2cds = function(pos,cds_int,strand) {
        if (strand==1) {
            return(which(pos==unlist(apply(cds_int, 1, function(x) x[1]:x[2]))))
        } else if (strand==-1) {
            return(which(pos==rev(unlist(apply(cds_int, 1, function(x) x[1]:x[2])))))
        }
    }
    
    # Annotating the functional impact of each substitution and populating the N matrices
    
    ref3_cod = mut3_cod = wrong_ref = aachange = ntchange = impact = array(NA, nrow(mutations))
    
    for (j in 1:nrow(mutations)) {
    
        geneind = mutations$geneind[j]
        pos = mutations$pos[j]
        if (any(pos == RefCDS[[geneind]]$intervals_splice)) { # Essential splice-site substitution
        
            impact[j] = "Essential_Splice"; impind = 4
            pos_ind = (pos==RefCDS[[geneind]]$intervals_splice)
            cdsnt = RefCDS[[geneind]]$seq_splice[pos_ind]
            ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], RefCDS[[geneind]]$seq_splice[pos_ind], RefCDS[[geneind]]$seq_splice1down[pos_ind])
            mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_splice1up[pos_ind], mutations$mut_cod[j], RefCDS[[geneind]]$seq_splice1down[pos_ind])
            aachange[j] = ntchange[j] = "-"
        
        } else { # Coding substitution
        
            pos_ind = chr2cds(pos, RefCDS[[geneind]]$intervals_cds, RefCDS[[geneind]]$strand)
            cdsnt = RefCDS[[geneind]]$seq_cds[pos_ind]
            ref3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], RefCDS[[geneind]]$seq_cds[pos_ind], RefCDS[[geneind]]$seq_cds1down[pos_ind])
            mut3_cod[j] = sprintf("%s%s%s", RefCDS[[geneind]]$seq_cds1up[pos_ind], mutations$mut_cod[j], RefCDS[[geneind]]$seq_cds1down[pos_ind])
            codon_pos = c(ceiling(pos_ind/3)*3-2, ceiling(pos_ind/3)*3-1, ceiling(pos_ind/3)*3)
            old_codon = as.character(as.vector(RefCDS[[geneind]]$seq_cds[codon_pos]))
            pos_in_codon = pos_ind-(ceiling(pos_ind/3)-1)*3
            new_codon = old_codon; new_codon[pos_in_codon] = mutations$mut_cod[j]
            old_aa = seqinr::translate(old_codon)
            new_aa = seqinr::translate(new_codon)
            aachange[j] = sprintf('%s%0.0f%s',old_aa,ceiling(pos_ind/3),new_aa)
            ntchange[j] = sprintf('%s%0.0f%s',mutations$ref_cod[j],pos_ind,mutations$mut_cod[j])
        
            # Annotating the impact of the mutation
            if (new_aa == old_aa){ 
                impact[j] = "Synonymous"; impind = 1
            } else if (new_aa == "*"){
                impact[j] = "Nonsense"; impind = 3
            } else if (old_aa != "*"){
                impact[j] = "Missense"; impind = 2
            } else if (old_aa=="*") {
                impact[j] = "Stop_loss"; impind = NA
            }
        }
        
        if (mutations$ref_cod[j] != as.character(cdsnt)) { # Incorrect base annotation in the input mutation file (the mutation will be excluded with a warning)
            wrong_ref[j] = 1
        } else if (!is.na(impind)) { # Correct base annotation in the input mutation file
            trisub = trinucsubsind[ paste(ref3_cod[j], mut3_cod[j], sep=">") ]
            RefCDS[[geneind]]$N[trisub,impind] = RefCDS[[geneind]]$N[trisub,impind] + 1 # Adding the mutation to the N matrices
        }
      
        if (round(j/1e4)==(j/1e4)) { message(sprintf('    %0.3g %%...', round(j/nrow(mutations),2)*100)) }
    }
    
    mutations$ref3_cod = ref3_cod
    mutations$mut3_cod = mut3_cod
    mutations$aachange = aachange
    mutations$ntchange = ntchange
    mutations$impact = impact
    mutations$pid = sapply(RefCDS,function(x) x$protein_id)[mutations$geneind]
    
    if (any(!is.na(wrong_ref))) {
        stop(sprintf('%0.0f mutations have a wrong reference base, please correct and rerun.',sum(!is.na(wrong_ref)))) # This can be made into a mere warning and the rest of the code will work
        wrong_refbase = mutations[!is.na(wrong_ref),]
        mutations = mutations[is.na(wrong_ref),]
    }
    
    if (any(nrow(indels))) { # If there are indels we concatenate the tables of subs and indels
        indels = cbind(indels, data.frame(ref_cod=".", mut_cod=".", strand=".", ref3_cod=".", mut3_cod=".", aachange=".", ntchange=".", impact="no-SNV", pid=sapply(RefCDS,function(x) x$protein_id)[indels$geneind]))
        annot = rbind(mutations, indels)
    } else {
        annot = mutations
    }
    annot = annot[order(annot$sampleID, annot$chr, annot$pos),]
    
    
    ## 3. Estimation of the global rate and selection parameters
    message("[3] Estimating global rates...")
    
    Lall = array(sapply(RefCDS, function(x) x$L), dim=c(192,4,length(RefCDS)))
    Nall = array(sapply(RefCDS, function(x) x$N), dim=c(192,4,length(RefCDS)))
    L = apply(Lall, c(1,2), sum)
    N = apply(Nall, c(1,2), sum)
    
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
    
    # Fitting all mutation rates and the 3 global selection parameters
     
    poissout = fit_substmodel(N, L, substmodel) # Original substitution model
    par = poissout$par
    poissmodel = poissout$model
    parmle =  setNames(par[,2], par[,1])
    mle_submodel = par
    rownames(mle_submodel) = NULL
    
    # Fitting models with 1 and 2 global selection parameters
    
    s1 = gsub("wmis","wall",gsub("wnon","wall",gsub("wspl","wall",substmodel)))
    par1 = fit_substmodel(N, L, s1)$par # Substitution model with 1 selection parameter
    s2 = gsub("wnon","wtru",gsub("wspl","wtru",substmodel))
    par2 = fit_substmodel(N, L, s2)$par # Substitution model with 1 selection parameter
    globaldnds = rbind(par, par1, par2)[c("wmis","wnon","wspl","wtru","wall"),]
    sel_loc = sel_cv = NULL
    
    ## 4. dNdSloc: variable rate dN/dS model (gene mutation rate inferred from synonymous subs in the gene only)
    
    genemuts = data.frame(gene_name = sapply(RefCDS, function(x) x$gene_name), n_syn=NA, n_mis=NA, n_non=NA, n_spl=NA, exp_syn=NA, exp_mis=NA, exp_non=NA, exp_spl=NA)
    genemuts[,2:5] = t(sapply(RefCDS, function(x) colSums(x$N)))
    mutrates = sapply(substmodel[,1], function(x) prod(parmle[base::strsplit(x,split="\\*")[[1]]])) # Expected rate per available site
    genemuts[,6:9] = t(sapply(RefCDS, function(x) colSums(x$L*mutrates)))
    numrates = length(mutrates)
    
    if (outp > 1) {
        message("[4] Running dNdSloc...")
        
        selfun_loc = function(j) {
            y = as.numeric(genemuts[j,-1])
            x = RefCDS[[j]]
            
            # a. Neutral model: wmis==1, wnon==1, wspl==1
            mrfold = sum(y[1:4])/sum(y[5:8]) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
            ll0 = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,1,1,1),dim=c(4,numrates))), log=T)) # loglik null model
            
            # b. Missense model: wmis==1, free wnon, free wspl
            mrfold = max(1e-10, sum(y[c(1,2)])/sum(y[c(5,6)])) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
            wfree = y[3:4]/y[7:8]/mrfold; wfree[y[3:4]==0] = 0
            llmis = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,1,wfree),dim=c(4,numrates))), log=T)) # loglik free wmis
            
            # c. free wmis, wnon and wspl
            mrfold = max(1e-10, y[1]/y[5]) # Correction factor of "t"
            w = y[2:4]/y[6:8]/mrfold; w[y[2:4]==0] = 0 # MLE of dN/dS based on the local rate (using syn muts as neutral)
            llall = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,w),dim=c(4,numrates))), log=T)) # loglik free wmis, wnon, wspl
            w[w>1e4] = 1e4
            
            p = 1-pchisq(2*(llall-c(llmis,ll0)),df=c(1,3))
            return(c(w,p))
        }
        
        sel_loc = as.data.frame(t(sapply(1:nrow(genemuts), selfun_loc)))
        colnames(sel_loc) = c("wmis_loc","wnon_loc","wspl_loc","pmis_loc","pall_loc")
        sel_loc$qmis_loc = p.adjust(sel_loc$pmis_loc, method="BH")
        sel_loc$qall_loc = p.adjust(sel_loc$pall_loc, method="BH")
        sel_loc = cbind(genemuts[,1:5],sel_loc)
        sel_loc = sel_loc[order(sel_loc$pall_loc,sel_loc$pmis_loc,-sel_loc$wmis_loc),]
    }
    
    
    ## 5. dNdScv: Negative binomial regression (with or without covariates) + local synonymous mutations
    
    nbreg = nbregind = NULL
    if (outp > 2) {
        
        message("[5] Running dNdScv...")
    
         # Covariates
        if (is.null(cv)) {
            nbrdf = genemuts[,c("n_syn","exp_syn")]
            model = glm.nb(n_syn ~ offset(log(exp_syn)) - 1 , data = nbrdf)
            message(sprintf("    Regression model for substitutions: no covariates were used (theta = %0.3g).", model$theta))
        } else {
            covs = covs[genemuts$gene_name,]
            if (ncol(covs) > maxcovs) {
                warning(sprintf("More than %s input covariates. Only the first %s will be considered.", maxcovs, maxcovs))
                covs = covs[,1:maxcovs]
            }
            nbrdf = cbind(genemuts[,c("n_syn","exp_syn")], covs)
            
            # Negative binomial regression
            model = suppressWarnings(glm.nb(n_syn ~ offset(log(exp_syn)) + . , data = nbrdf))
            if (!is.null(model$th.warn) | nrow(genemuts)<500) { # If there are warnings or if <500 genes, we run the regression without covariates
                model = glm.nb(n_syn ~ offset(log(exp_syn)) - 1 , data = nbrdf)
                message(sprintf("    Regression model for substitutions: no covariates were used (theta = %0.3g).", model$theta))
            } else {
                message(sprintf("    Regression model for substitutions: all covariates were used (theta = %0.3g).", model$theta))
            }
        }
        if (all(model$y==genemuts$n_syn)) {
            genemuts$exp_syn_cv = model$fitted.values
        }
        theta = model$theta
        nbreg = model
        
        # Subfunction: Analytical opt_t using only neutral subs
        mle_tcv = function(n_neutral, exp_rel_neutral, shape, scale) {
            tml = (n_neutral+shape-1)/(exp_rel_neutral+(1/scale))
            if (shape<=1) { # i.e. when theta<=1
                tml = max(shape*scale,tml) # i.e. tml is bounded to the mean of the gamma (i.e. y[9]) when theta<=1, since otherwise it takes meaningless values
            }
            return(tml)
        }
        
        # Subfunction: dNdScv per gene
        selfun_cv = function(j) {
            y = as.numeric(genemuts[j,-1])
            x = RefCDS[[j]]
            exp_rel = y[5:8]/y[5]
            # Gamma
            shape = theta
            scale = y[9]/theta
            
            # a. Neutral model
            indneut = 1:4 # vector of neutral mutation types under this model (1=synonymous, 2=missense, 3=nonsense, 4=essential_splice)
            opt_t = mle_tcv(n_neutral=sum(y[indneut]), exp_rel_neutral=sum(exp_rel[indneut]), shape=shape, scale=scale)
            mrfold = max(1e-10, opt_t/y[5]) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
            ll0 = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,1,1,1),dim=c(4,numrates))), log=T)) + dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik null model
            
            # b. Missense model: wmis==1, free wnon, free wspl
            indneut = 1:2
            opt_t = mle_tcv(n_neutral=sum(y[indneut]), exp_rel_neutral=sum(exp_rel[indneut]), shape=shape, scale=scale)
            mrfold = max(1e-10, opt_t/sum(y[5])) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
            wfree = y[3:4]/y[7:8]/mrfold; wfree[y[3:4]==0] = 0
            llmis = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,1,wfree),dim=c(4,numrates))), log=T)) + dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik free wmis
            
            # c. Truncating muts model: free wmis, wnon==wspl==1
            indneut = c(1,3,4)
            opt_t = mle_tcv(n_neutral=sum(y[indneut]), exp_rel_neutral=sum(exp_rel[indneut]), shape=shape, scale=scale)
            mrfold = max(1e-10, opt_t/sum(y[5])) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
            wfree = y[2]/y[6]/mrfold; wfree[y[2]==0] = 0
            lltrunc = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,wfree,1,1),dim=c(4,numrates))), log=T)) + dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik free wmis
            
            # d. Free selection model: free wmis, free wnon, free wspl
            indneut = 1
            opt_t = mle_tcv(n_neutral=sum(y[indneut]), exp_rel_neutral=sum(exp_rel[indneut]), shape=shape, scale=scale)
            mrfold = max(1e-10, opt_t/sum(y[5])) # Correction factor of "t" based on the obs/exp ratio of "neutral" mutations under the model
            wfree = y[2:4]/y[6:8]/mrfold; wfree[y[2:4]==0] = 0
            llall_unc = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,wfree),dim=c(4,numrates))), log=T)) + dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik free wmis
          
            if (constrain_wnon_wspl == 0) {
            
                p = 1-pchisq(2*(llall_unc-c(llmis,lltrunc,ll0)),df=c(1,2,3))
                return(c(wfree,p))
                
            } else { # d2. Free selection model: free wmis, free wnon==wspl
            
                wmisfree = y[2]/y[6]/mrfold; wmisfree[y[2]==0] = 0
                wtruncfree = sum(y[3:4])/sum(y[7:8])/mrfold; wtruncfree[sum(y[3:4])==0] = 0
                llall = sum(dpois(x=x$N, lambda=x$L*mutrates*mrfold*t(array(c(1,wmisfree,wtruncfree,wtruncfree),dim=c(4,numrates))), log=T)) + dgamma(opt_t, shape=shape, scale=scale, log=T) # loglik free wmis, free wnon==wspl
                p = 1-pchisq(2*c(llall_unc-llmis,llall-c(lltrunc,ll0)),df=c(1,1,2))
                return(c(wmisfree,wtruncfree,wtruncfree,p))
            }
        }
        
        sel_cv = as.data.frame(t(sapply(1:nrow(genemuts), selfun_cv)))
        colnames(sel_cv) = c("wmis_cv","wnon_cv","wspl_cv","pmis_cv","ptrunc_cv","pallsubs_cv")
        sel_cv$qmis_cv = p.adjust(sel_cv$pmis_cv, method="BH")
        sel_cv$qtrunc_cv = p.adjust(sel_cv$ptrunc_cv, method="BH")
        sel_cv$qallsubs_cv = p.adjust(sel_cv$pallsubs_cv, method="BH")
        sel_cv = cbind(genemuts[,1:5],sel_cv)
        sel_cv = sel_cv[order(sel_cv$pallsubs_cv, sel_cv$pmis_cv, sel_cv$ptrunc_cv, -sel_cv$wmis_cv),] # Sorting genes in the output file
    
        ## Indel recurrence: based on a negative binomial regression (ideally fitted excluding major known driver genes)
        
        if (nrow(indels) >= min_indels) {
        
            geneindels = as.data.frame(array(0,dim=c(length(RefCDS),8)))
            colnames(geneindels) = c("gene_name","n_ind","n_induniq","n_indused","cds_length","excl","exp_unif","exp_indcv")
            geneindels$gene_name = sapply(RefCDS, function(x) x$gene_name)
            geneindels$n_ind = as.numeric(table(indels$gene)[geneindels[,1]]); geneindels[is.na(geneindels[,2]),2] = 0
            geneindels$n_induniq = as.numeric(table(unique(indels[,-1])$gene)[geneindels[,1]]); geneindels[is.na(geneindels[,3]),3] = 0
            
            if (use_indel_sites) {
                geneindels$n_indused = geneindels[,3]
            } else {
                geneindels$n_indused = geneindels[,2]
            }
            geneindels$cds_length = sapply(RefCDS, function(x) x$CDS_length)
            geneindels$excl = (geneindels[,1] %in% known_cancergenes)
            if (sum(geneindels[!geneindels$excl,"n_indused"]) == 0) { # If there are no indels for the background model we do not exclude any gene
                geneindels$excl = F 
            }
            geneindels$exp_unif = sum(geneindels[!geneindels$excl,"n_indused"]) / sum(geneindels[!geneindels$excl,"cds_length"]) * geneindels$cds_length
          
            # Negative binomial regression for indels
            
            if (is.null(cv)) {
                nbrdf = geneindels[,c("n_indused","exp_unif")][!geneindels[,6],] # We exclude known drivers from the fit
                model = suppressWarnings(glm.nb(n_indused ~ offset(log(exp_unif)) + . , data = nbrdf))
                model = glm.nb(n_indused ~ offset(log(exp_unif)) - 1 , data = nbrdf)
                message(sprintf("    Regression model for indels: no covariates were used (theta = %0.3g)", model$theta))
                nbrdf_all = geneindels[,c("n_indused","exp_unif")]
            } else {
                nbrdf = cbind(geneindels[,c("n_indused","exp_unif")], covs)[!geneindels[,6],] # We exclude known drivers from the fit
                model = suppressWarnings(glm.nb(n_indused ~ offset(log(exp_unif)) + . , data = nbrdf))
                if (!is.null(model$th.warn) | nrow(genemuts)<500) { # If there are warnings or if <500 genes, we run the regression without covariates
                    model = glm.nb(n_indused ~ offset(log(exp_unif)) - 1 , data = nbrdf)
                    message(sprintf("    Regression model for indels: no covariates were used (theta = %0.3g)", model$theta))
                } else {
                    message(sprintf("    Regression model for indels: all covariates were used (theta = %0.3g)", model$theta))
                }
                nbrdf_all = cbind(geneindels[,c("n_indused","exp_unif")], covs)
            }
             
            theta_indels = model$theta
            nbregind = model
            geneindels$exp_indcv = exp(predict(model,nbrdf_all))
            geneindels$wind = geneindels$n_indused / geneindels$exp_indcv
        
            # Statistical testing for indel recurrence per gene
        
            geneindels$pind = pnbinom(q=geneindels$n_indused-1, mu=geneindels$exp_indcv, size=theta_indels, lower.tail=F)
            geneindels$qind = p.adjust(geneindels$pind, method="BH")
        
            # Fisher combined p-values (substitutions and indels)
            
            sel_cv = merge(sel_cv, geneindels, by="gene_name")[,c("gene_name","n_syn","n_mis","n_non","n_spl","n_indused","wmis_cv","wnon_cv","wspl_cv","wind","pmis_cv","ptrunc_cv","pallsubs_cv","qallsubs_cv","pind")]
            colnames(sel_cv) = c("gene_name","n_syn","n_mis","n_non","n_spl","n_ind","wmis_cv","wnon_cv","wspl_cv","wind_cv","pmis_cv","ptrunc_cv","pallsubs_cv","qallsubs_cv","pind_cv")
            sel_cv$pglobal_cv = 1 - pchisq(-2 * (log(sel_cv$pallsubs_cv) + log(sel_cv$pind_cv)), df = 4)
            sel_cv$qglobal_cv = p.adjust(sel_cv$pglobal, method="BH")
        
            sel_cv = sel_cv[order(sel_cv$pglobal_cv, sel_cv$pallsubs_cv, sel_cv$pmis_cv, sel_cv$ptrunc_cv, -sel_cv$wmis_cv),] # Sorting genes in the output file
        }
    }
    
    dndscvout = list(globaldnds = globaldnds, sel_cv = sel_cv, sel_loc = sel_loc, annotmuts = annot, genemuts = genemuts, mle_submodel = mle_submodel, exclsamples = exclsamples, exclmuts = exclmuts, nbreg = nbreg, nbregind = nbregind, poissmodel = poissmodel)
    
} # EOF