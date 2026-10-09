"""Step 8 - known partners of every RWE gene: are they co-dependent in NextGen and in DepMap 2D?

Known partners are taken from the curated complex annotation of each RWE gene (data/curated/rwe_genes.tsv, column `complex`): every gene
symbol named there that belongs to the dependency-gene family (step 4).  The same rule is applied to every modification family.
For each pair:
  r and call in NextGen and DepMap 2D (step 4 statistics)
  DepMap 2D down-sampled to 147 screens (300 draws): median r and the fraction of draws that would be called at the NextGen call
  threshold (p <= NextGen BH boundary and |r| > 0.3)  -> separates sample size from data differences
  NextGen r after regressing out the screen-level pan-essential depletion axis (mean effect of the paper's Pan-dependency genes)
Outputs: results/tables/known_partners.tsv, known_partners_by_family.tsv
"""
import os
import re

import numpy as np
import pandas as pd
from _common import banner, data, paths, save, stats

DRAWS = 300
banner('8. known partners: NextGen vs DepMap 2D')
rng = np.random.default_rng(147)
ge = data.screen_effect()
tda = data.tda_classes()
cls = pd.read_csv(os.path.join(paths.TABLES, 'dependency_classes.tsv'), sep='\t').set_index('key')
fam_keys = [k for k in cls.index[cls.in_codependency_family] if k in ge.columns]
by_sym = {data.symbol_of(k): k for k in fam_keys}
genes = data.rwe_genes().dropna(subset=['key']).drop_duplicates('key')

pairs = set()
for g in genes.itertuples():
    if g.key not in set(fam_keys) or not isinstance(g.complex, str):
        continue
    for tok in set(re.findall(r'\b[A-Z][A-Z0-9]+\b', g.complex)):          # hyphens split names ("METTL3-METTL14")
        if tok in by_sym and by_sym[tok] != g.key:
            pairs.add(tuple(sorted((g.key, by_sym[tok]))))
pairs = sorted(pairs)
family_of = genes.set_index('key').family
print(f'known partner pairs from the curated complex annotations (both genes dependencies): {len(pairs)}')

ov = pd.read_csv(os.path.join(paths.TABLES, 'codep_overview.tsv'), sep='\t').set_index('dataset')
p_star = float(ov.loc['NextGen', 'p_boundary'])
ng, tr = list(data.nextgen_screens()), list(data.traditional_lineage_matched_screens())
pan = [k for k in tda.index[tda.NextGenClass == 'Pan-dependency'] if k in ge.columns]
axis = ge.loc[ng, pan].mean(1).values
G_ng, G_2d = ge.loc[ng], ge.loc[tr]
idx = np.array([rng.choice(len(tr), 147, replace=False) for _ in range(DRAWS)])
rows = []
for a, b in pairs:
    r = {}
    for d, G in (('NextGen', G_ng), ('DepMap2D', G_2d)):
        rr, p, n = stats.pearson_vs_all(G[a].values, G[[b]].values)
        r[d] = (rr[0], p[0], n[0])
    xa, xb = G_2d[a].values, G_2d[b].values
    sub = np.array([stats.pearson_vs_all(xa[i], xb[i][:, None])[0][0] for i in idx])
    called_sub = (stats.r_to_p(sub, 147) <= p_star) & (np.abs(sub) > stats.R_THRESHOLD)
    va, vb = G_ng[a].values, G_ng[b].values
    m = ~np.isnan(va) & ~np.isnan(vb)
    X = np.column_stack([np.ones(m.sum()), axis[m]])
    res = lambda v: v - X @ np.linalg.lstsq(X, v, rcond=None)[0]  # noqa: E731
    rwe_gene = a if a in family_of.index else b
    rows.append(dict(gene_a=data.symbol_of(a), gene_b=data.symbol_of(b), family=family_of.get(rwe_gene), r_NextGen=r['NextGen'][0], r_DepMap2D=r['DepMap2D'][0],
                     DepMap2D_n147_median_r=float(np.median(sub)), DepMap2D_n147_frac_called=float(called_sub.mean()),
                     r_NextGen_axis_adjusted=float(np.corrcoef(res(va[m]), res(vb[m]))[0, 1]), key_a=a, key_b=b))
K = pd.DataFrame(rows)
calls = pd.read_csv(os.path.join(paths.TABLES, 'codep_rwe_calls.tsv'), sep='\t')
called = {d: set(map(tuple, np.sort(calls.loc[calls[f'call_{d}'], ['gene_a', 'gene_b']].values, axis=1))) for d in ('NextGen', 'DepMap2D')}
for d in ('NextGen', 'DepMap2D'):
    K[f'call_{d}'] = [tuple(sorted((a, b))) in called[d] for a, b in zip(K.key_a, K.key_b)]
save(K.sort_values(['family', 'r_DepMap2D'], ascending=[True, False]), 'tables', 'known_partners.tsv', fmt='%.4g')
S = K.groupby('family').agg(pairs=('r_NextGen', 'size'), median_r_NextGen=('r_NextGen', 'median'), median_r_DepMap2D=('r_DepMap2D', 'median'),
                            called_NextGen=('call_NextGen', 'mean'), called_DepMap2D=('call_DepMap2D', 'mean'), called_DepMap2D_at_n147=('DepMap2D_n147_frac_called', 'mean'))
tot = pd.DataFrame([dict(pairs=len(K), median_r_NextGen=K.r_NextGen.median(), median_r_DepMap2D=K.r_DepMap2D.median(), called_NextGen=K.call_NextGen.mean(),
                         called_DepMap2D=K.call_DepMap2D.mean(), called_DepMap2D_at_n147=K.DepMap2D_n147_frac_called.mean())], index=['all families'])
S = pd.concat([S, tot])
S.index.name = 'family'
save(S.reset_index(), 'tables', 'known_partners_by_family.tsv', fmt='%.3g')
print(S.to_string(float_format=lambda v: f'{v:.2f}'))
co = K[K.r_DepMap2D >= 0.3]
print(f'\npairs co-dependent in 2D (r >= 0.3): {len(co)} | median r 2D {co.r_DepMap2D.median():.2f} -> NextGen {co.r_NextGen.median():.2f} | '
      f'2D down-sampled to 147 screens still called in {co.DepMap2D_n147_frac_called.mean():.0%} of draws (mean over pairs)')
