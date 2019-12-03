# GSEA

The original enrichment analysis method from the Broad Institute.

### Intro

* t-tests all the livelong day have no biological context. Also miss set/pathway effects (several small expression changes may have a larger impact than a single large expression change)
* GSEA generates an enrichment score of a gene set, rather than individual genes
* Gene sets as knowledge-based

### Proc

Given an input table with (nrows, ncols) = (n_features, n_samples) (NB "feature" == "gene"):
* the samples would be experimental groups - assume only 2 groups. Generate a "score" of some kind between groups, e.g. t-test.
* rank rows by this score. def. **reference signature** is the list of rownames ranked in descending order of score.
* compare ref sig to a query gene set, which is just a list of genes in a gene set. For each row in ref sig, check if row `i` is in query gene set.
  * If hit: ES += score[i]
  * If miss: ES -= some meaningful way to decrease ES, usually related to size of query set.
* **Filtering** - the first few rows in the ref sig are the genes that contribute the most to the ES.
* **nominal p-value** - Generate a null distribution by label permutation.
* **Normalized ES** - some formula to normalize between gene sets of different sizes.

~End~
