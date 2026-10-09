# VisRWE · DepMap

**English** · [中文](README.zh-CN.md)

Dependencies and co-dependencies of RNA-modification **writers, readers and erasers (RWE)** in the DepMap NextGen organoid / neurosphere
CRISPR screens, side by side with the original DepMap 2D screens.

**Interactive explorer → <https://shuyang-zhao.github.io/VisRWE_DepMap/>** (single static page; no installation)

![Explorer](docs/assets/explorer.png)

## What is here

* **217 RWE genes** (m6A, Ψ, m5C, m7G, 2′-O-methylation, tRNA modifications, …) with role, modification, substrates and caveats.
* **Single-gene dependencies** in NextGen (147 screens) and DepMap 2D (597 screens of the same lineages): dependency class, effect size,
  organoid vs 2D, lineage selectivity, biomarkers.
* **A co-dependency map of all 5,693 dependency genes**, computed the same way in both datasets (16.2 million pairs each), with every pair that
  involves an RWE gene labelled *called in both / NextGen only / DepMap 2D only*.
* **An explorer** that starts from the gene you enter, shown alone at the centre: click it to add its called partners (closer = larger |r|),
  then click partners to grow the network level by level (hover previews what a click adds), open any edge to see
  its lineage decomposition (forest plot and scatter plots, NextGen next to 2D, plus a Cas12a-library check), switch between NextGen, DepMap 2D
  and a comparison view, or look at functional-module enrichment as a heat map. Chinese / English interface.

The aim is practical: anyone working on an RNA-modification enzyme can look up its dependency profile and partners in both datasets without
processing DepMap themselves.

## Statistics (aligned with the NextGen paper)

The statistical conventions follow the paper's public code ([broadinstitute/DepMap-NextGen-Public](https://github.com/broadinstitute/DepMap-NextGen-Public)):
screen-level Chronos gene effect; NextGen and 2D screen sets as defined there; the paper's dependency classes (re-implemented and reproduced
for 100 % of genes, then applied to 2D); Pearson correlation on pairwise-complete screens; one Benjamini–Hochberg correction per results table;
a co-dependency is **called at FDR < 0.005 and |r| > 0.3**. The paper's published expression-addiction table is recomputed from the raw data as
a check of the data handling. Details: [docs/methods.md](docs/methods.md).

## Main observations

* About a third of RWE genes are pan-dependencies in both datasets (writer cofactors 57 %, erasers 0 %); classes agree closely between NextGen and 2D.
* Known partners of RWE genes (138 pairs named in the curated complex annotations, all modification families) are co-dependent far less often
  in NextGen (14 %) than in 2D (40 %); down-sampling 2D to 147 screens hardly changes the 2D rate (35 %), so the loss reflects the NextGen data,
  not sample size or the statistics. It is seen across modification families, not in one machinery ([report §4](docs/report.zh-CN.md)).
* Module over-representation among each RWE gene's top partners (16 broad categories and 1,439 Reactome / CORUM modules, tested identically for
  every gene) points mostly to the machineries the RWE genes themselves belong to: mitochondrial translation, tRNA processing / modification and
  rRNA modification.
* 35 associations are strong inside a single NextGen lineage only (33 in CNS/brain).

Full report (Chinese): [docs/report.zh-CN.md](docs/report.zh-CN.md) · figures: [results/figures](results/figures) · tables:
[results/tables](results/tables) ([data dictionary](docs/data_dictionary.md)) · workbook: [results/VisRWE_DepMap_tables.xlsx](results/VisRWE_DepMap_tables.xlsx)

## Repository layout

```
├── data/
│   ├── curated/          RWE gene list, control genes, literature claim audit (shipped)
│   ├── raw/              DepMap / figshare / gene-set downloads (not versioned; see data/README.md)
│   └── raw.md5           checksums of the raw inputs
├── src/rwe_nextgen/      data access (paper conventions), statistics, dependency classification
├── pipeline/             numbered analysis steps 01–12
├── web/                  explorer source (HTML template, CSS, JavaScript, i18n, vendored vis-network)
├── docs/                 GitHub Pages site (index.html) + methods, report, data dictionary
└── results/              tables, figures, workbook (results/intermediate/ is re-created by the pipeline)
```

## Reproduce

```bash
conda env create -f environment.yml && conda activate visrwe
# download the inputs listed in data/README.md into data/raw/, then
make check      # verify checksums
make all        # analysis -> explorer -> figures -> report (about 10 minutes, ~4 GB RAM)
```

| Step | Script | Output |
|---|---|---|
| 1 | `01_inputs_and_gene_list.py` | input check, validation against the paper, RWE gene list |
| 2 | `02_dependency_classes.py` | dependency classes (NextGen from the paper; 2D with the same rules) |
| 3 | `03_single_gene.py` | single-gene profiles, organoid vs 2D, lineage tests, biomarkers, library check |
| 4 | `04_codependency.py` | genome-wide co-dependency in both datasets |
| 5 | `05_lineage_restricted.py` | associations restricted to one NextGen lineage |
| 6 | `06_modules.py` | functional categories and Reactome / CORUM module over-representation (NextGen, 2D, Cas12a 2D) |
| 7 | `07_literature.py` | published claims vs data |
| 8 | `08_known_partners.py` | known partners of every RWE gene: NextGen vs 2D, down-sampling |
| 9–10 | `09_explorer_data.py`, `10_build_site.py` | data bundle and `docs/index.html` |
| 11 | `11_figures.py` | figures |
| 12 | `12_report.py` | report and workbook |

## Citation

If you use this resource, please cite it (see [CITATION.cff](CITATION.cff)):

> Zhao S. *VisRWE: dependencies and co-dependencies of RNA-modification writers, readers and erasers in DepMap NextGen organoid and DepMap 2D
> CRISPR screens.* GitHub, 2026. <https://github.com/shuyang-zhao/VisRWE_DepMap>

and the data it is built on:

* Neiswender, Maffa, Brenan *et al.* A dependency map enhanced with next-generation 3D cancer models. *Nature* (2026). doi:10.1038/s41586-026-10843-7
* Tsherniak A. *et al.* Defining a cancer dependency map. *Cell* 170, 564–576 (2017).
* Dempster J. M. *et al.* Chronos: a cell population dynamics model of CRISPR experiments that improves inference of gene fitness effects. *Genome Biology* 22, 343 (2021).
* DepMap portal, Broad Institute — <https://depmap.org>, file set "NextGen Model Manuscript 2026".

## Author

**Shuyang Zhao** — Department of Chemistry, University of California, Riverside.

## Licence

Code: MIT (see [LICENSE](LICENSE)). DepMap data are released by the Broad Institute under CC BY 4.0; the tables and the quantised gene-effect
vectors embedded in the explorer are derived from them. Gene-set libraries keep the licences of their sources. vis-network (bundled in the page)
is dual-licensed Apache-2.0 / MIT.
