#' Generate a RefCDS object
#' 
#' Function to build a RefCDS object from a reference genome and a table of transcripts. The RefCDS object has to be precomputed for any new species or assembly prior to running dndscv. This function generates an .rda file that needs to be input into dndscv using the refdb argument. Note that when multiple CDS share the same gene name (second column of cdsfile), the longest coding CDS will be chosen for the gene. CDS with ambiguous bases (N) will not be considered.
#' 
#' @author Inigo Martincorena (Wellcome Sanger Institute)
#' @details Martincorena I, et al. (2017) Universal patterns of selection in cancer and somatic tissues. Cell. 171(5):1029-1041.
#' 
#' @param cdsfile Path to the reference transcript table.
#' @param genomefile Path to the indexed reference genome file.
#' @param outfile Output file name (default = "RefCDS.rda").
#' @param numcode NCBI genetic code number (default = 1; standard genetic code). To see the list of genetic codes supported use: ? seqinr::translate
#' @param excludechrs Vector or string with chromosome names to be excluded from the RefCDS object (default: no chromosome will be excluded). The mitochondrial chromosome should be excluded as it has different genetic code and mutation rates, either using the excludechrs argument or not including mitochondrial transcripts in cdsfile.
#' @param onlychrs Vector of valid chromosome names (default: all chromosomes will be included)
#' 
#' @export

buildref = function(cdsfile, genomefile, outfile = "RefCDS.rda", numcode = 1, excludechrs = NULL, onlychrs = NULL) {
    
    ## 1. Valid chromosomes and reference CDS per gene
    message("[1/3] Preparing the environment...")
    
    reftable = read.table(cdsfile, header=1, sep="\t", stringsAsFactors=F)
    colnames(reftable) = c("gene.id","gene.name","cds.id","chr","chr.coding.start","chr.coding.end","cds.start","cds.end","length","strand")
    
    validchrs = sapply(strsplit(system(sprintf("grep '^>' %s", genomefile), intern=T), split=" "), function(x) substr(x[1],2,nchar(x[1])))
    validchrs = setdiff(validchrs, excludechrs)
    if (length(onlychrs)>0) {
        validchrs = validchrs[validchrs %in% onlychrs]
    }
    
    # Restricting to chromosomes present in both the genome file and the CDS table
    if (any(validchrs %in% unique(reftable$chr))) {
        validchrs = validchrs[validchrs %in% unique(reftable$chr)]
    } else { # Try adding a chr prefix
        reftable$chr = paste("chr", reftable$chr, sep="")
        validchrs = validchrs[validchrs %in% unique(reftable$chr)]
        if (length(validchrs)==0) { # No matching chromosome names
            stop("No chromosome names in common between the genome file and the CDS table")
        }
    }

    # Selecting the longest complete CDS for every gene (required when there are multiple alternative transcripts per unique gene name)
    
    reftable = reftable[reftable[,1]!="" & reftable[,2]!="" & reftable[,3]!="" & !is.na(reftable[,5]) & !is.na(reftable[,6]),] # Removing invalid entries
    reftable = reftable[which(reftable$chr %in% validchrs),] # Only valid chromosomes
    
    cds_table = unique(reftable[,c(1:3,9)])
    cds_table = cds_table[order(cds_table$gene.name, -cds_table$length), ] # Sorting CDS from longest to shortest
    cds_table = cds_table[(cds_table$length %% 3)==0, ] # Removing CDS of length not multiple of 3
    fullcds = intersect(reftable$cds.id[reftable$cds.start==1], reftable$cds.id[reftable$cds.end==reftable$length]) # List of complete CDS
    
    cds_table = cds_table[cds_table$cds.id %in% fullcds, ] # Complete CDS
    reftable = reftable[reftable$cds.id %in% fullcds, ] # Complete CDS
    gene_list = unique(cds_table$gene.name)
    
    reftable = reftable[order(reftable$chr, reftable$chr.coding.start), ]
    cds_split = split(reftable, f=reftable$cds.id)
    gene_split = split(cds_table, f=cds_table$gene.name)
    
    
    ## 2. Building the RefCDS object
    message("[2/3] Building the RefCDS object...")
    
    # Subfunction to extract the coding sequence
    get_CDSseq = function(gr, strand) {
        cdsseq = strsplit(paste(as.vector(Rsamtools::scanFa(genomefile, gr)),collapse=""),"")[[1]]
        if (strand==-1) {
            cdsseq = rev(seqinr::comp(cdsseq,forceToLower=F))
        }
        return(cdsseq)
    }
    
    # Subfunction to extract essential splice site sequences
    # Definition of essential splice sites: (5' splice site: +1,+2,+5; 3' splice site: -1,-2)
    get_splicesites = function(cds) {
        splpos = numeric(0)
        if (nrow(cds)>1) { # If the CDS has more than one exon
            if (cds[1,10]==1) { # + strand
                spl5prime = cds[-nrow(cds),6] # Exon end before splice site
                spl3prime = cds[-1,5] # Exon start after splice site
                splpos = unique(sort(c(spl5prime+1, spl5prime+2, spl5prime+5, spl3prime-1, spl3prime-2)))
            } else if (cds[1,10]==-1) { # - strand
                spl5prime = cds[-1,5] # Exon end before splice site
                spl3prime = cds[-nrow(cds),6] # Exon start after splice site
                splpos = unique(sort(c(spl5prime-1, spl5prime-2, spl5prime-5, spl3prime+1, spl3prime+2)))
            }
        }
        return(splpos)
    }
    
    # Subfunction to extract the essential splice site sequence
    get_spliceseq = function(gr, strand) {
        spliceseq = unname(as.vector(Rsamtools::scanFa(genomefile, gr)))
        if (strand==-1) {
            spliceseq = seqinr::comp(spliceseq,forceToLower=F)
        }
        return(spliceseq)
    }
    
    # Initialising and populating the RefCDS object
    
    RefCDS = array(list(NULL), length(gene_split)) # Initialising empty object
    invalid_genes = rep(0, length(gene_split)) # Initialising empty object
    
    for (j in 1:length(gene_split)) {
        
        gene_cdss = gene_split[[j]]
        h = keeptrying = 1
        
        while (h<=nrow(gene_cdss) & keeptrying) {
            
            pid = gene_cdss[h,3]
            cds = cds_split[[pid]]
            strand = cds[1,10]
            chr = cds[1,4]
            gr = GenomicRanges::GRanges(chr, IRanges::IRanges(cds[,5], cds[,6]))
            cdsseq = get_CDSseq(gr,strand)
            pseq = seqinr::translate(cdsseq, numcode = numcode)
            
            if (all(pseq[-length(pseq)]!="*") & all(cdsseq!="N")) { # A valid CDS has been found (no stop codons inside the CDS excluding the last codon) and no "N" nucleotides
                
                # Essential splice sites
                splpos = get_splicesites(cds) # Essential splice sites
                if (length(splpos)>0) { # CDSs with a single exon do not have splice sites
                    gr_spl = GenomicRanges::GRanges(chr, IRanges::IRanges(splpos, splpos))
                    splseq = get_spliceseq(gr_spl, strand)
                }
                
                # Obtaining the splicing sequences and the coding and splicing sequence contexts
                if (strand==1) {
                    
                    cdsseq1up = get_CDSseq(GenomicRanges::GRanges(chr, IRanges::IRanges(cds[,5]-1, cds[,6]-1)), strand)
                    cdsseq1down = get_CDSseq(GenomicRanges::GRanges(chr, IRanges::IRanges(cds[,5]+1, cds[,6]+1)), strand)
                    if (length(splpos)>0) {
                        splseq1up = get_spliceseq(GenomicRanges::GRanges(chr, IRanges::IRanges(splpos-1, splpos-1)), strand)
                        splseq1down = get_spliceseq(GenomicRanges::GRanges(chr, IRanges::IRanges(splpos+1, splpos+1)), strand)
                    }
                    
                } else if (strand==-1) {
                    
                    cdsseq1up = get_CDSseq(GenomicRanges::GRanges(chr, IRanges::IRanges(cds[,5]+1, cds[,6]+1)), strand)
                    cdsseq1down = get_CDSseq(GenomicRanges::GRanges(chr, IRanges::IRanges(cds[,5]-1, cds[,6]-1)), strand)
                    if (length(splpos)>0) {
                        splseq1up = get_spliceseq(GenomicRanges::GRanges(chr, IRanges::IRanges(splpos+1, splpos+1)), strand)
                        splseq1down = get_spliceseq(GenomicRanges::GRanges(chr, IRanges::IRanges(splpos-1, splpos-1)), strand)
                    }
                    
                }
                
                # Annotating the CDS in the RefCDS database
                
                RefCDS[[j]]$gene_name = gene_cdss[h,2]
                RefCDS[[j]]$gene_id = gene_cdss[h,1]
                RefCDS[[j]]$protein_id = gene_cdss[h,3]
                RefCDS[[j]]$CDS_length = gene_cdss[h,4]
                RefCDS[[j]]$chr = cds[1,4]
                RefCDS[[j]]$strand = strand
                RefCDS[[j]]$intervals_cds = unname(as.matrix(cds[,5:6]))
                RefCDS[[j]]$intervals_splice = splpos
                
                RefCDS[[j]]$seq_cds = Biostrings::DNAString(paste(cdsseq, collapse=""))
                RefCDS[[j]]$seq_cds1up = Biostrings::DNAString(paste(cdsseq1up, collapse=""))
                RefCDS[[j]]$seq_cds1down = Biostrings::DNAString(paste(cdsseq1down, collapse=""))
                
                if (length(splpos)>0) { # If there are splice sites in the gene
                    RefCDS[[j]]$seq_splice = Biostrings::DNAString(paste(splseq, collapse=""))
                    RefCDS[[j]]$seq_splice1up = Biostrings::DNAString(paste(splseq1up, collapse=""))
                    RefCDS[[j]]$seq_splice1down = Biostrings::DNAString(paste(splseq1down, collapse=""))
                }
                
                keeptrying = 0 # Stopping the while loop
            }
            h = h+1
        }
        if (keeptrying) {
            invalid_genes[j] = 1 # No valid CDS was found for this gene and the gene will be removed from the RefCDS object
        }
        if (round(j/1000)==(j/1000)) { message(sprintf('    %0.3g%% ...', round(j/length(gene_split),2)*100)) }
    }
    
    RefCDS = RefCDS[!invalid_genes] # Removing genes without a valid CDS
    
    
    ## 3. L matrices: number of synonymous, missense, nonsense and splice sites in each CDS at each trinucleotide context
    message("[3/3] Calculating the impact of all possible coding changes...")
    
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
    
    for (j in 1:length(RefCDS)) {
        
        L = array(0, dim=c(192,4))
        cdsseq = as.character(as.vector(RefCDS[[j]]$seq_cds))
        cdsseq1up = as.character(as.vector(RefCDS[[j]]$seq_cds1up))
        cdsseq1down = as.character(as.vector(RefCDS[[j]]$seq_cds1down))
        
        # 1. Exonic mutations
        
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
        
        # Synonymous
        matrix_ind = table(matrind[which(imp==1)])
        L[as.numeric(names(matrix_ind)), 1] = matrix_ind
        
        # Missense
        matrix_ind = table(matrind[which(imp==2)])
        L[as.numeric(names(matrix_ind)), 2] = matrix_ind
        
        # Nonsense
        matrix_ind = table(matrind[which(imp==3)])
        L[as.numeric(names(matrix_ind)), 3] = matrix_ind
        
        # 2. Splice site mutations
        if (length(RefCDS[[j]]$intervals_splice)>0) {
            splseq = as.character(as.vector(RefCDS[[j]]$seq_splice))
            splseq1up = as.character(as.vector(RefCDS[[j]]$seq_splice1up))
            splseq1down = as.character(as.vector(RefCDS[[j]]$seq_splice1down))
            old_trinuc = rep(paste(splseq1up, splseq, splseq1down, sep=""), each=3)
            new_trinuc = paste(rep(splseq1up, each=3), c(sapply(splseq, function(x) nt[nt!=x])), rep(splseq1down,each=3), sep="")
            matrind = trinuc_subsind[paste(old_trinuc, new_trinuc, sep=">")]
            matrix_ind = table(matrind)
            L[as.numeric(names(matrix_ind)), 4] = matrix_ind
        }
        
        RefCDS[[j]]$L = L # Saving the L matrix
        if (round(j/1000)==(j/1000)) { message(sprintf('    %0.3g%% ...', round(j/length(gene_split),2)*100)) }
    }
    
    ## Saving the reference GenomicRanges object

    aux = unlist(sapply(1:length(RefCDS), function(x) t(cbind(x,rbind(RefCDS[[x]]$intervals_cds,cbind(RefCDS[[x]]$intervals_splice,RefCDS[[x]]$intervals_splice))))))
    df_genes = as.data.frame(t(array(aux,dim=c(3,length(aux)/3))))
    colnames(df_genes) = c("ind","start","end")
    df_genes$chr = unlist(sapply(1:length(RefCDS), function(x) rep(RefCDS[[x]]$chr,nrow(RefCDS[[x]]$intervals_cds)+length(RefCDS[[x]]$intervals_splice))))
    df_genes$gene = sapply(RefCDS, function(x) x$gene_name)[df_genes$ind]
    
    gr_genes = GenomicRanges::GRanges(df_genes$chr, IRanges::IRanges(df_genes$start, df_genes$end))
    GenomicRanges::mcols(gr_genes)$names = df_genes$gene
    
    save(RefCDS, gr_genes, file=outfile)
    
} # EOF
