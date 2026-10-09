# Methods

All statistics follow the public code of the NextGen paper (Neiswender, Maffa, Brenan *et al.*, Nature 2026;
<https://github.com/broadinstitute/DepMap-NextGen-Public>, `src/generate.py` and `src/data_utils.py`). Where this project goes beyond the paper
(the co-dependency map itself, the 2D comparison, lineage decomposition), the same conventions are reused.

## Data

| Item | Definition |
|---|---|
| Gene effect | Screen-level Chronos gene effect, `screen_gene_effect.csv` (1,534 screens × 18,435 genes). Missing values are kept; every correlation uses the screens where both genes have data. |
| NextGen | `ScreenType ∈ {2DO, 3DO, 3DN, 2DN}`, `IsNextGen`, `OncotreeLineage ≠ Other` → 147 screens (organoids and neurospheres; all Cas12a Humagne-CD). |
| DepMap 2D | `IsTraditional2D` screens of the NextGen lineages (the paper's `adherent_next_gen_lineage_matched_screens`) → 597 screens (586 Cas9 Avana / KY, 11 Humagne-CD). |
| Cas12a 2D | All traditional 2D screens made with the Humagne-CD library (36), used only as a library check. |
| Expression | Model-level log2(TPM + 1), expanded to screens; screens without expression are dropped where expression is used. |

**Validation.** Step 1 recomputes the paper's published expression-addiction table from the raw matrices (17,889 genes; agreement to rounding
precision, identical *n*). Step 2 re-implements the paper's dependency classification and reproduces the published NextGen classes for 100 %
of 18,435 genes before applying the same rules to the DepMap 2D screens.

## Single-gene analyses (step 3)

* **Dependency classes** – Pan-dependency, Strongly Selective, High Variance, Weakly Selective, Non-dependency (paper rules; `min_n_dep = 3`).
* **Organoids vs adherent 2D** – the paper's table (`organoids_vs_2d_gene_dependency_tests_full.csv`; Mann–Whitney U, BH over all genes; genes
  essential in both groups are flagged `DropEssential` and excluded from calls). Call: FDR < 0.05 and |mean difference| > 0.3.
* **Lineage vs rest** – two-sided Mann–Whitney U per lineage (≥ 5 screens), one BH per dataset over the RWE × lineage table; call: FDR < 0.05
  and mean difference < −0.3.
* **Biomarkers** – the paper's expression-addiction and paralog tables with the paper's family (Strongly Selective / High Variance /
  Pan-dependency genes; paralog expression variance > 1), BH, and calls (FDR < 0.005 and r < −0.3, resp. r > 0.3).
* **Library check** – 24 models screened with both Cas9 and Cas12a: paired difference per gene (one-sample t-test, BH).
* **Set enrichment** – Fisher's exact test of dependency classes for RWE families / roles / caveat classes against the genome, BH per dataset.

## Co-dependency (step 4)

* **Family** – genes that are a dependency (Strongly Selective, High Variance or Pan-dependency) in NextGen *or* DepMap 2D: 5,693 genes,
  16,202,278 unique pairs. The same pairs are tested in both datasets.
* **Test** – Pearson r on pairwise-complete screens; two-sided p from a *t* distribution with *n* − 2 degrees of freedom (`scipy.stats.pearsonr`).
* **Multiplicity** – one Benjamini–Hochberg correction over all pairs of a dataset (`scipy.stats.false_discovery_control`).
* **Call** – FDR < 0.005 and |r| > 0.3 (the paper's thresholds for correlation-type associations). At these settings the NextGen boundary is
  p ≈ 2 × 10⁻⁵ and the DepMap 2D boundary p ≈ 3 × 10⁻⁴.
* **Comparison class** – called in both / NextGen only / DepMap 2D only; Fisher-z difference P between the two r values as an annotation.
* **2D library flag** – within-lineage KY-vs-Avana R² ≥ 0.1 marks genes whose 2D gene effect is largely a library difference (e.g. SNU13);
  their 2D edges are flagged.

## Lineage-restricted associations (step 5, NextGen only)

For pairs with an RWE gene, Pearson r inside each NextGen lineage with ≥ 14 screens (CNS/Brain 39, Esophagus/Stomach 27, Pancreas 23,
Bowel 22, Breast 14); one BH over all pair × lineage tests; call: FDR < 0.005, |r| > 0.3, and a Fisher-z difference P < 0.01 between the lineage
and the sample-size-weighted mean of the other large lineages. Classes: *restricted* (other lineages r ≤ 0.3 in the same direction),
*enriched* (0.3–0.5), *shared, stronger here* (> 0.5). Not done in 2D because the 2D lineages mix libraries in different proportions.

## Functional modules (step 6)

Two levels, defined identically for every gene (no hand-made or hypothesis-specific modules):

* **categories** – 16 broad functional categories assigned to each gene by regular expressions on its GO BP / GO CC / Reactome / KEGG / CORUM
  terms (first match in a fixed order); used for node colours and the overview heat map;
* **pathways** – every Reactome 2022 pathway and CORUM complex with 5–150 members among the dependency genes (identical member sets kept once;
  1,439 modules).

Test (the paper's `run_rank_hypergeometric` convention): the 100 most positively correlated partners of each RWE gene, universe = dependency
genes, one-sided hypergeometric test per module; BH per dataset and level; call FDR < 0.05. Run in NextGen, DepMap 2D and Cas12a 2D (library
check). As in the paper, the over-representation test does not model correlation among module members.

## Known partners (step 8)

Known partners are the gene symbols named in each RWE gene's curated complex annotation that are dependency genes (138 pairs, all modification
families). For each pair: r and call in both datasets; DepMap 2D down-sampled to 147 screens (300 draws) to separate sample size from data
differences; and NextGen r after regressing out the screen-level pan-essential depletion axis. Summaries per modification family.

## Explorer (steps 9–10)

The page embeds gene attributes, all called pairs with an RWE gene plus the 20 strongest called partners of every family gene in each dataset,
and quantised gene-effect vectors (int8, 0.025 units) of all bundle genes in NextGen, DepMap 2D and Cas12a 2D. Per-edge statistics in the
browser (lineage decomposition, scatter plots, Cochran Q) are computed from these vectors; displayed r values for edges come from the
pipeline (full precision).

## Caveats

* NextGen uses one library (Cas12a); a NextGen-only signal must be checked in Cas12a 2D lines before being attributed to 3D culture.
* No covariates (as in the paper): part of a pooled correlation can come from lineage composition or general essential machinery; the lineage
  decomposition shows this per edge.
* Co-dependency is not interaction. Many RWE dependencies reflect catalysis-independent structural roles (`caveat_class`).
