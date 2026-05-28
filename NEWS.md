# pepitope 0.4.0

* Added `plot_read_structure()` to identify and visualize read features
* Read feature positions and strand are now identified automatically by default
* Deprecate `plot_reads()`, rename to `plot_read_count()`
* Deprecate `plot_distr()`, rename to `plot_read_distr()`
* Replace `fqtk` and `guide-counter` tools with internal one-pass `count_fastq()`

# pepitope 0.3.3

* Speed up runtime by caching `txdb`-derived objects and reusing them

# pepitope 0.3.2

* Add a workaround for `EnsDb`/`BSGenome` mismatch on `UCSC` names (#3)
* Fix a bug where multi-exon dropping could subset incorrectly (#5)
* Re-introduce filter that CDS width must be equal variant width (#5)
* Fix error message when sample sheet does not contain required fields
* Add clearer error messages about `seqlevels` mismatches in `annotate_coding()` (#11)
* Barcode counting now uses `--exact-match` by default

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
