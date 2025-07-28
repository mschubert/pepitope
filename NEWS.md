# git head

* Add a workaround for `EnsDb`/`BSGenome` mismatch on `UCSC` names (#3)
* Fix a bug where multi-exon dropping could subset incorrectly (#5)
* Re-introduce filter that CDS width must be equal variant width (#5)

# pepitope 0.3.1

* Multi-exon variants are now dropped with a warning (#5)
* FASTQ demultiplexing is now done in a clean directory (#7)
* `fqtk` and `guide-counter` are now supplied via R packages
* `fqtk` and `guide-counter` now handle spaces in paths correctly (#8)
* `filter_variants` will now error if provided more than one sample implicitly

# pepitope 0.3.0

* Added QC and screen incl. plotting functionality to the package
* Added command-line wrappers for `fqtk` demultiplexing and `guide-counter`
* Provided usage vignettes for Variant calling, QC and Co-culture screen

# pepitope 0.2.0

* Added support for RNA-based gene fusion using the `rnafusion` NF-core pipeline

# pepitope 0.1.0
 
* Provided simple variant annotation from the `sarek` NF-core pipeline
