# pepitope 0.4.3

* Use `minigenes` for sample construct libraries and barcode library membership
* Refactor `remove_cutsite()` into `Rcpp` utility function
* Refactor and 1.6x (3.5x cached) speed up `annodate_coding()` (1000 SNVs tested)
* Refactor and 2.2x speed up `count_fastq()` (400M .gz reads tested)

# pepitope 0.4.2

* Add unknown sample(s) to `plot_read_count()`
* Plotting functions now take `label_format` argument for sample labels
* Reduced the number of direct package dependencies

# pepitope 0.4.1

* Partial read structures can now be inferred (but not counted)
* `plot_read_structure()` now indicates if a barcode is valid
* `plot_read_structure()` now displays the barcode type counts on each position

# pepitope 0.4.0

* Added `plot_read_structure()` to identify and visualize read features
* Read feature positions and strand are now identified automatically by default
* `mut_id` is now shared among variants on the same position (#6)
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
