# Changelog

## pepitope 0.4.0

- Added [`plot_read_structure()`](../reference/plot_read_structure.md)
  to identify and visualize read features
- Read feature positions and strand are now identified automatically by
  default
- Deprecate [`plot_reads()`](../reference/plot_read_count.md), rename to
  [`plot_read_count()`](../reference/plot_read_count.md)
- Deprecate [`plot_distr()`](../reference/plot_read_distr.md), rename to
  [`plot_read_distr()`](../reference/plot_read_distr.md)
- Replace `fqtk` and `guide-counter` tools with internal one-pass
  [`count_fastq()`](../reference/count_fastq.md)

## pepitope 0.3.3

- Speed up runtime by caching `txdb`-derived objects and reusing them

## pepitope 0.3.2

- Add a workaround for `EnsDb`/`BSGenome` mismatch on `UCSC` names
  ([\#3](https://github.com/mschubert/pepitope/issues/3))
- Fix a bug where multi-exon dropping could subset incorrectly
  ([\#5](https://github.com/mschubert/pepitope/issues/5))
- Re-introduce filter that CDS width must be equal variant width
  ([\#5](https://github.com/mschubert/pepitope/issues/5))
- Fix error message when sample sheet does not contain required fields
- Add clearer error messages about `seqlevels` mismatches in
  [`annotate_coding()`](../reference/annotate_coding.md)
  ([\#11](https://github.com/mschubert/pepitope/issues/11))
- Barcode counting now uses `--exact-match` by default

## pepitope 0.3.1

- Multi-exon variants are now dropped with a warning
  ([\#5](https://github.com/mschubert/pepitope/issues/5))
- FASTQ demultiplexing is now done in a clean directory
  ([\#7](https://github.com/mschubert/pepitope/issues/7))
- `fqtk` and `guide-counter` are now supplied via R packages
- `fqtk` and `guide-counter` now handle spaces in paths correctly
  ([\#8](https://github.com/mschubert/pepitope/issues/8))
- `filter_variants` will now error if provided more than one sample
  implicitly

## pepitope 0.3.0

- Added QC and screen incl. plotting functionality to the package
- Added command-line wrappers for `fqtk` demultiplexing and
  `guide-counter`
- Provided usage vignettes for Variant calling, QC and Co-culture screen

## pepitope 0.2.0

- Added support for RNA-based gene fusion using the `rnafusion` NF-core
  pipeline

## pepitope 0.1.0

- Provided simple variant annotation from the `sarek` NF-core pipeline
