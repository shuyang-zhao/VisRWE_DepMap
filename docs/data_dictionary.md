# Data dictionary (`results/tables/`)

Gene keys are DepMap keys `SYMBOL (ENTREZ)`. Datasets: `NextGen` (147 organoid / neurosphere screens), `DepMap2D` (597 traditional 2D screens of
the NextGen lineages). r = Pearson correlation of gene effects across screens; FDR = Benjamini–Hochberg.

| File | Rows | Key columns |
|---|---|---|
| `rwe_gene_list.tsv` | RWE gene | `family` (modification), `role` (writer / writer-cofactor / reader / reader-cofactor / eraser / other), `evidence_tier`, `caveat_class`, `NextGenClass_paper` |
| `validation_expression_addiction.tsv` | gene | recomputed vs published r and n (validation of data handling) |
| `dependency_classes.tsv` | gene | `NextGenClass` (paper), `NextGenClass_reimplemented`, `DepMap2DClass` (same rules), `in_codependency_family` |
| `single_gene_rwe.tsv` | RWE gene | `<group>_mean / _sd / _frac_dep` for NextGen, Organoid3D, CNS, DepMap2D; `org_vs_2D_*` (paper table; `org_vs_2D_call`); `library_*` (dual-library models); `<dataset>_top_lineage*`; `biomarkers` |
| `lineage_tests_rwe.tsv` | RWE gene × lineage | `<dataset>_mean_difference`, `<dataset>_FDR`, `<dataset>_hit`, `class` |
| `biomarkers_rwe.tsv` | test | paper's expression-addiction / paralog rows for RWE genes with recomputed FDR and `Significant` |
| `set_enrichment.tsv` | set × outcome × dataset | Fisher odds ratio, p, FDR |
| `codep_overview.tsv` | dataset | pairs tested, calls, p at the FDR boundary |
| `codep_rwe_calls.tsv` | gene pair | `r_/p_/FDR_/n_/call_<dataset>`, `class` (both / NextGen only / DepMap 2D only), `rwe_rwe`, `DepMap2D_library_flag`, `p_difference` |
| `codep_summary_per_gene.tsv` | RWE gene | `n_calls_<dataset>`, `n_calls_both`, `top_partners_<dataset>` |
| `library_sensitivity_2D.tsv` | family gene | `DepMap2D_library_R2` (within-lineage KY vs Avana), `DepMap2D_library_sensitive` |
| `lineage_restricted_pairs.tsv` | pair × lineage | `r`, `FDR` within the lineage, `r_other_lineages`, `p_difference`, `restriction`, r in every large lineage |
| `gene_categories.tsv` | gene | functional `category` (16 + unannotated) |
| `module_enrichment_categories.tsv` | RWE gene × category | `overlap_ / odds_ / p_ / FDR_ / call_ / mean_r_<dataset>` for NextGen, DepMap2D, Cas12a2D; `class` |
| `module_enrichment_pathways.tsv` | RWE gene × Reactome / CORUM module | same columns; rows called in NextGen or DepMap 2D (all rows in `results/intermediate/`) |
| `known_partners.tsv` | pair | `family`, r and call in both datasets, DepMap 2D down-sampled to 147 screens (median r, fraction called), axis-adjusted NextGen r |
| `known_partners_by_family.tsv` | modification family | pairs, median r and fraction called per dataset, 2D at 147 screens |
| `literature_rwe.tsv` | RWE gene | `claim_stance`, `concordance_<dataset>`, audit verdict counts |
| `literature_claims.tsv` | claim | citation, PMID, claimed lineage / phenotype, verdict, current classes |

`results/intermediate/` (not versioned) holds the full pair tables and correlation matrices needed by steps 8–11.
