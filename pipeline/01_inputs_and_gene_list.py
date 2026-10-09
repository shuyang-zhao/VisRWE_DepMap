"""Step 1 - check the inputs, the screen sets and the RWE gene list.

* the NextGen and DepMap-2D screen sets are built exactly as in the paper's code (see src/rwe_nextgen/data.py)
* validation: the paper's own expression-addiction table (figshare, expression_addiction_tests_full.csv) is recomputed from the raw
  matrices; agreement shows that screen sets, gene keys and matrix handling are identical to the paper
* the curated RWE list (data/curated/rwe_genes.tsv) is mapped to DepMap gene keys
Outputs: results/tables/rwe_gene_list.tsv, results/tables/validation_expression_addiction.tsv
"""
import os

import pandas as pd
from _common import banner, data, paths, save, stats

banner('1. inputs, screen sets, RWE gene list')
need = [os.path.join(paths.PORTAL, f) for f in ('screen_gene_effect.csv', 'screen_gene_dependency.csv', 'next_gen_expression.csv', 'traditional_expression.csv')] + \
       [os.path.join(paths.FIGSHARE, f) for f in ('screen_metadata.csv', 'tda_dependency_classification.csv', 'expression_addiction_tests_full.csv', 'paralog_dependency_tests_full.csv',
                                                 'organoids_vs_2d_gene_dependency_tests_full.csv')] + \
       [os.path.join(paths.GENESETS, f + '.gmt') for f in ('GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'Reactome_2022', 'KEGG_2021_Human', 'CORUM')] + \
       [os.path.join(paths.CURATED, f) for f in ('rwe_genes.tsv', 'rwe_control_genes.tsv', 'literature_claim_audit.tsv')]
missing = [f for f in need if not os.path.exists(f)]
for f in need:
    print(('  ok      ' if f not in missing else '  MISSING ') + os.path.relpath(f, paths.ROOT))
if missing:
    raise SystemExit('missing inputs - see data/README.md')

ge = data.screen_effect()
sm = data.screen_meta()
ng, tr = list(data.nextgen_screens()), list(data.traditional_lineage_matched_screens())
print(f'\nscreen_gene_effect: {ge.shape[0]:,} screens x {ge.shape[1]:,} genes')
print(f'NextGen screens: {len(ng)}  ({sm.loc[ng, "ScreenType"].value_counts().to_dict()})')
print(f'DepMap 2D comparator (traditional 2D screens of the NextGen lineages): {len(tr)}  ({sm.loc[tr, "Library2"].value_counts().to_dict()})')
print('lineages (NextGen):', sm.loc[ng, 'OncotreeLineage'].value_counts().to_dict())
print('lineages (DepMap 2D):', sm.loc[tr, 'OncotreeLineage'].value_counts().to_dict())

# ---- validation against the paper's published expression-addiction correlations
paper = pd.read_csv(os.path.join(paths.FIGSHARE, 'expression_addiction_tests_full.csv'))
ex = data.to_screens(data.expression(), ng)
G = ge.loc[ng]
rows = []
for g in paper.Dependency:
    if g in G.columns and g in ex.columns:
        r, p, n = stats.pearson_vs_all(G[g].values, ex[[g]].values)
        rows.append((g, float(r[0]), int(n[0])))
chk = pd.DataFrame(rows, columns=['Dependency', 'r_recomputed', 'n_recomputed']).merge(paper[['Dependency', 'NextGenPearsonR', 'NextGenN']], on='Dependency').dropna()
dr = (chk.r_recomputed - chk.NextGenPearsonR).abs()
print(f'\nvalidation - paper expression-addiction table recomputed for {len(chk):,} genes: max |delta r| = {dr.max():.1e}, median {dr.median():.1e}; n identical for {(chk.n_recomputed == chk.NextGenN).mean():.1%}')
save(chk, 'tables', 'validation_expression_addiction.tsv')

# ---- RWE gene list
u = data.rwe_genes().copy()
tda = data.tda_classes()
u['in_gene_effect'] = u.key.isin(ge.columns)
u['NextGenClass_paper'] = u.key.map(tda.NextGenClass)
u['OrganoidsClass_paper'] = u.key.map(tda.OrganoidsClass)
cols = ['symbol', 'key', 'hgnc_symbol', 'entrez', 'family', 'role', 'evidence_tier', 'caveat_class', 'modifications', 'substrates', 'complex', 'source', 'in_lab_list',
        'in_gene_effect', 'NextGenClass_paper', 'OrganoidsClass_paper', 'cancer_dependency_claims', 'key_refs', 'notes']
save(u[cols], 'tables', 'rwe_gene_list.tsv')
print(f'RWE genes: {len(u)} | with a DepMap key: {u.key.notna().sum()} | in the gene-effect matrix: {int(u.in_gene_effect.sum())}')
print('roles:', u.role.value_counts().to_dict())
print('NextGen classes (paper):', u.NextGenClass_paper.value_counts().to_dict())
c = data.control_genes()
print(f'control genes (non-RNA methyltransferases / IGFBPs / ALKBH): {len(c)}, in matrix {int(c.key.isin(ge.columns).sum())}')
