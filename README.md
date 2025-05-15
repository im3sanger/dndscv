dndscv
=====

Description
---
The **dNdScv** R package is a group of maximum-likelihood dN/dS methods designed to 
	quantify selection in cancer and somatic evolution (Martincorena *et al.*, 2017). The
	package contains functions to quantify dN/dS ratios for missense, nonsense and 
	essential splice mutations, at the level of individual genes, groups of genes or at 
	whole-exome level. The *dndscv* function within the package was originally designed 
	to detect cancer driver genes (*i.e.* genes under positive selection in cancer genomes)
	on datasets ranging from a few samples to thousands of samples, in whole-exome/genome 
	or targeted sequencing studies. 
	
Although initially designed for cancer genomic studies, this package can also be used 
    with appropriate caution to study selection in other resequencing studies, such 
    as SNP analyses, mutation accumulation studies in bacteria or for the discovery 
    of mutations causing developmental disorders using data from human trios. Please 
    study the optional arguments carefully if you are using the dndscv function for 
    other applications.
	
When using the dndscv function in the package (sel_cv output object), the background 
    mutation rate of each gene is estimated by combining local information 
	(synonymous mutations in the gene) and global information (variation of the mutation 
	rate across genes, exploiting epigenomic covariates), and controlling for the sequence 
	composition of the gene and mutational signatures. This allows to perform more
	sensitive inferences of positive selection in sparse datasets. Unlike traditional 
	implementations of dN/dS using Markov-chain models, mutations in dNdScv are modelled
	as Poisson events (or negative binomial events), which allows the use of more
	complex substitution models and the estimation of dN/dS ratios for truncating
	mutations. By default, *dNdScv* uses a trinucleotide context-dependent substitution 
	model, which is important to avoid common biases affecting simpler substitution 
	models in dN/dS (Greenman *et al.*, 2006, and Martincorena *et al*, 2017).

Installation
--------
You can use devtools::install_github() to install *dndscv* from this repository:

	> library(devtools); install_github("im3sanger/dndscv")

Tutorial
--------
For a tutorial on dNdScv see the vignette included with the package. This includes 
examples for whole-exome/genome data and for targeted data.

[Tutorial: getting started with dNdScv](http://htmlpreview.github.io/?http://github.com/im3sanger/dndscv/blob/master/vignettes/dNdScv.html)

Genome assemblies and species
--------
By default, *dNdScv* assumes that mutation data is mapped to the GRCh37/hg19 assembly of the
human genome. If you are using human data mapped to the GRCh38/hg38 assembly, you can use 
refdb="hg38" as an argument in dndscv to use the default GRCh38/hg38 precomputed database
and epigenomic covariates (please ensure that you have downloaded the latest version of
dNdScv).

Users interested in trying *dNdScv* on a different set of transcripts, a different human assembly
or a different species can use the buildref function to create a custom RefCDS, as explained
in this [tutorial](http://htmlpreview.github.io/?http://github.com/im3sanger/dndscv/blob/master/vignettes/buildref.html).

Pre-computed RefCDS files (RefCDS objects) to run *dNdScv* on some popular species
(e.g. mouse, rat, cow, dog, yeast or SARS-CoV-2) are available from 
[this link](https://github.com/im3sanger/dndscv_data/tree/master/data).

Reference
----
Martincorena I, *et al*. (2017) Universal Patterns of Selection in Cancer and Somatic Tissues. *Cell*.
http://www.cell.com/cell/fulltext/S0092-8674(17)31136-4

Author
--------
Inigo Martincorena.

Acknowledgements
--------
Moritz Gerstung and Peter Campbell for advice and inspiration. Federico Abascal and Andrew Lawson for testing, feedback and ideas.
