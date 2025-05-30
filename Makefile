.PHONY: all
all: doc vignettes test

R = R --no-save --no-restore -e
PKGVER = $(shell grep Version: < DESCRIPTION | sed "s/Version: //")

.PHONY: test
test:
	$(R) "devtools::test()"

.PHONY: check
check:
	$(R) "devtools::check()"

rmd_files=$(wildcard vignettes/*.rmd)
knit_results=$(patsubst vignettes/%.rmd,inst/doc/%.md,$(rmd_files))

.PHONY: vignettes
vignettes: inst/doc ${knit_results}
	$(R) "library(knitr); library(devtools); build_vignettes()"

inst/doc:
	mkdir -p $@

inst/doc/%.md: vignettes/%.rmd
	$(R) "knitr::knit('$<', '$@')"

.PHONY: doc
doc:
	$(R) "devtools::document()"

.PHONY: package
package: doc vignettes
	R CMD build .
	R CMD check pepitope_$(PKGVER).tar.gz
	$(R) "BiocCheck::BiocCheck()"

.PHONY: deploy
deploy:
	$(R) "pkgdown::deploy_to_branch()"

.PHONY: clean
clean:
	${RM} -r inst/doc
	${RM} -r man
