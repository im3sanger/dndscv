# Package definitions for dndscv
#
# Author: Inigo Martincorena
###############################################################################

#' Detection of selection in cancer and somatic evolution
#'
#' The dNdScv R package is a suite of maximum-likelihood dN/dS methods designed
#' to quantify selection in cancer and somatic evolution (Martincorena et al., 2017).
#' The package contains functions to quantify dN/dS ratios for missense, nonsense
#' and essential splice mutations, at the level of individual genes, groups of
#' genes or at whole-genome level. The dNdScv method was designed to detect cancer
#' driver genes (i.e. genes under positive selection in cancer) on datasets ranging
#' from a few samples to thousands of samples, in whole-exome/genome or targeted
#' sequencing studies.
#' @name dndscv-package
#' @docType package
#' @title Detection of selection in cancer and somatic evolution
#' @author Inigo Martincorena, Wellcome Trust Sanger Institute, \email{im3@@sanger.ac.uk}
#' @references Martincorena I, et al. (2017) Universal Patterns of Selection in Cancer and Somatic Tissues. Cell.
#' @keywords package
#' @seealso \code{\link{dndscv}}
#' @seealso \code{\link{buildref}}
#' @import seqinr
#' @import MASS
#' @import GenomicRanges
#' @import IRanges
#' @import Rsamtools
#' @import Biostrings
NA
