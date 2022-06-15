dndscv
=====

Description
---
The **dNdScv** R package is a group of maximum-likelihood dN/dS methods designed to 
	quantify selection in cancer and somatic evolution (Martincorena *et al.*, 2017). The 
	package contains functions to quantify dN/dS ratios for missense, nonsense and 
	essential splice mutations, at the level of individual genes, groups of genes or at 
	whole-exome level. The *dndscv* function within the package was designed to detect cancer driver genes 
	(*i.e.* genes under positive selection in cancer) on datasets ranging from a few 
	samples to thousands of samples, in whole-exome/genome or targeted sequencing studies. 
	
Although initially designed for cancer genomic studies, this package can also be used to quantify
	selection in other resequencing studies, such as SNP analyses, mutation accumulation 
	studies in bacteria or for the discovery of mutations causing developmental disorders 
	using data from human trios.
	
The background mutation rate of each gene is estimated by combining local information 
	(synonymous mutations in the gene) and global information (variation of the mutation 
	rate across genes, exploiting epigenomic covariates), and controlling for the sequence 
	composition of the gene and mutational signatures. Unlike traditional implementations 
	of dN/dS, *dNdScv* uses trinucleotide context-dependent substitution matrices to 
	avoid common mutation biases affecting dN/dS (Greenman *et al.*, 2006).

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

Precomputed reference files (RefCDS objects) to run *dNdScv* on other popular assemblies 
(e.g. GRCh38/hg38) or species (e.g. mouse, rat, cow, dog, yeast or SARS-CoV-2) are
available for download from [this link](https://github.com/im3sanger/dndscv_data/tree/master/data).

Reference
----
Martincorena I, *et al*. (2017) Universal Patterns of Selection in Cancer and Somatic Tissues. *Cell*.
http://www.cell.com/cell/fulltext/S0092-8674(17)31136-4

Acknowledgements
--------

Moritz Gerstung and Peter Campbell.

Federico Abascal for extensive testing, feedback and ideas.
