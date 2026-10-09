# Input data

The raw inputs are public but large, so they are not part of the repository. Put them under `data/raw/` with the layout below and
check them with `md5sum -c raw.md5` (run inside `data/`). Step 1 of the pipeline also lists anything missing.

```
data/
├── curated/                              shipped with the repository
│   ├── rwe_genes.tsv                     217 RNA-modification writers / readers / erasers with DepMap gene keys and annotations
│   ├── rwe_control_genes.tsv             32 control genes (non-RNA methyltransferases, IGFBPs, ALKBH family)
│   └── literature_claim_audit.tsv        86 PMID-verified published cancer-dependency claims with the audit verdict
└── raw/                                  downloaded, not versioned
    ├── depmap_portal_nextgen/            DepMap portal, file set "NextGen Model Manuscript 2026"
    │   ├── screen_gene_effect.csv
    │   ├── screen_gene_dependency.csv
    │   ├── next_gen_expression.csv
    │   └── traditional_expression.csv
    ├── figshare_nextgen_29472362/        figshare doi:10.6084/m9.figshare.29472362 (processed tables of the paper)
    │   ├── screen_metadata.csv
    │   ├── tda_dependency_classification.csv
    │   ├── expression_addiction_tests_full.csv
    │   ├── paralog_dependency_tests_full.csv
    │   └── organoids_vs_2d_gene_dependency_tests_full.csv
    └── genesets/                         Enrichr gene-set libraries (https://maayanlab.cloud/Enrichr/#libraries)
        ├── GO_Biological_Process_2023.gmt
        ├── GO_Cellular_Component_2023.gmt
        ├── Reactome_2022.gmt
        ├── KEGG_2021_Human.gmt
        └── CORUM.gmt
```

## Where to get the files

| Folder | Source |
|---|---|
| `depmap_portal_nextgen/` | <https://depmap.org/portal/data_page/?tab=allData> → choose the file set **"NextGen Model Manuscript 2026"**. The portal blocks scripted downloads, so download the four files in a browser. |
| `figshare_nextgen_29472362/` | <https://doi.org/10.6084/m9.figshare.29472362> (Neiswender, Maffa, Brenan *et al.*, Nature 2026). |
| `genesets/` | Enrichr library downloads (text files renamed to `.gmt`). |

## Curated files

* **`rwe_genes.tsv`** – union of a curated RWE catalogue (modification, role, substrates, evidence tier, caveat class, cancer-dependency claims
  with references) and an independent lab list of RNA-modification regulators; `key` is the DepMap gene key `SYMBOL (ENTREZ)`.
* **`rwe_control_genes.tsv`** – genes the lab list excludes on purpose (non-RNA methyltransferases, IGFBPs, ALKBH family), used as a control set.
* **`literature_claim_audit.tsv`** – published claims that an RWE gene is a cancer dependency, each with citation / PMID, claimed lineage and
  phenotype, and the verdict of a claim-by-claim check against the NextGen single-gene data (supported / partially supported / not supported /
  contradicted / not testable).

## Licences

DepMap data are distributed by the Broad Institute under CC BY 4.0; please cite the DepMap and NextGen publications when using
derived tables (see the main README). Enrichr libraries carry the licences of their original sources (GO, Reactome, KEGG, CORUM).
