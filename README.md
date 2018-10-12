dndscv
=====

Description
---
The **dNdScv** R package is a suite of maximum-likelihood dN/dS methods designed to 
	quantify selection in cancer and somatic evolution (Martincorena *et al.*, 2017). The 
	package contains functions to quantify dN/dS ratios for missense, nonsense and 
	essential splice mutations, at the level of individual genes, groups of genes or at 
	whole-genome level. The *dNdScv* method was designed to detect cancer driver genes 
	(*i.e.* genes under positive selection in cancer) on datasets ranging from a few 
	samples to thousands of samples, in whole-exome/genome or targeted sequencing studies. 
	
The background mutation rate of each gene is estimated by combining local information 
	(synonymous mutations in the gene) and global information (variation of the mutation 
	rate across genes, exploiting epigenomic covariates), and controlling for the sequence 
	composition of the gene and mutational signatures. Unlike traditional implementations 
	of dN/dS, *dNdScv* uses trinucleotide context-dependent substitution matrices to 
	avoid common mutation biases affecting dN/dS (Greenman *et al.*, 2006).

Note
----
The latest version of this package (from 12 October 2018) includes support for other 
human genome assemblies and other species.

Installation
--------
You can use devtools::install_github() to install *dndscv* from this repository:

	> library(devtools); install_github("im3sanger/dndscv")

Tutorial
--------
For a tutorial on dNdScv see the vignette included with the package. This includes 
examples for whole-exome/genome data and for targeted data.

[Tutorial: getting started with dNdScv](http://htmlpreview.github.io/?http://github.com/im3sanger/dndscv/blob/master/vignettes/dNdScv.html)

By default, dNdScv assumes that mutation data is mapped to the GRCh37/hg19 assembly of the
human genome. Users interested in trying dNdScv on a different set of transcripts, a
different assembly or a different species can follow this [tutorial](http://htmlpreview.github.io/?http://github.com/im3sanger/dndscv/blob/master/vignettes/buildref.html).

Reference
----
Martincorena I, *et al*. (2017) Universal Patterns of Selection in Cancer and Somatic Tissues. *Cell*.
http://www.cell.com/cell/fulltext/S0092-8674(17)31136-4

Acknowledgements
--------

Moritz Gerstung and Peter Campbell.