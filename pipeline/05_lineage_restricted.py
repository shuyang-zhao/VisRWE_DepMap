"""Step 5 - associations that are strong inside ONE NextGen lineage only ("lineage-restricted").

The pooled maps (step 4) remain the discovery layer; this step looks inside each lineage with >= 14 screens for pairs with an RWE gene:
  test   : Pearson r within the lineage, two-sided p; one BH over all (pair x lineage) tests of a dataset
  call   : FDR < 0.005 and |r| > 0.3 within the lineage, AND the lineage differs from the other large lineages
           (Fisher-z difference vs their sample-size-weighted mean, P < 0.01)
  class  : restricted (other lineages r <= 0.3 in the same direction) / enriched (0.3-0.5) / shared but stronger here (> 0.5)
NextGen only: every NextGen screen uses the same library, whereas the DepMap 2D lineages mix Avana, KY and Humagne screens in different
proportions, so within-lineage correlations in 2D are confounded by library composition (e.g. SNU13: 66% of its within-lineage 2D
variance is the KY-vs-Avana difference).  Per-lineage 2D correlations are still shown next to NextGen in the explorer.
Output: results/tables/lineage_restricted_pairs.tsv
"""
import os

import numpy as np
import pandas as pd
from _common import banner, data, paths, save, stats
from scipy import stats as ss

MIN_N = 14
banner('5. lineage-restricted associations')
cls = pd.read_csv(os.path.join(paths.TABLES, 'dependency_classes.tsv'), sep='\t').set_index('key')
ge, sm = data.screen_effect(), data.screen_meta()
fam = [k for k in cls.index[cls.in_codependency_family] if k in ge.columns]
pos = {k: i for i, k in enumerate(fam)}
rwe = [k for k in fam if k in set(data.rwe_genes().key.dropna())]
ri = np.array([pos[k] for k in rwe])
out = []
for dname, screens in (('NextGen', data.nextgen_screens),):
    scr = list(screens())
    lin = sm.loc[scr, 'OncotreeLineage']
    big = [L for L, n in lin.value_counts().items() if n >= MIN_N]
    Rs, Ns = {}, {}
    for L in big:
        X = ge.loc[[s for s in scr if lin[s] == L], fam].values
        Rs[L], Ns[L] = stats.pairwise_r(X, rows=ri)
        for j, i in enumerate(ri):
            Rs[L][j, i] = np.nan
    # p and BH over all (pair x lineage) tests; RWE-RWE pairs appear twice in the RWE-row matrices -> keep each unordered pair once
    keep = np.ones((len(ri), len(fam)), bool)
    for j, i in enumerate(ri):
        for j2, i2 in enumerate(ri):
            if i2 < i:
                keep[j, i2] = False
    P = {L: stats.r_to_p(Rs[L].astype(float), Ns[L].astype(float)) for L in big}
    flat = np.concatenate([np.where(keep, P[L], np.nan).ravel() for L in big])
    Q = stats.bh(flat).reshape(len(big), len(ri), len(fam))
    for li, L in enumerate(big):
        r, n, q = Rs[L].astype(float), Ns[L].astype(float), Q[li]
        others = [o for o in big if o != L]
        w = np.stack([np.maximum(Ns[o].astype(float) - 3, 0) for o in others])
        z = np.stack([np.arctanh(np.clip(np.nan_to_num(Rs[o].astype(float)), -0.999, 0.999)) for o in others])
        zr = (w * z).sum(0) / w.sum(0)
        r_rest = np.tanh(zr)
        se = np.sqrt(1 / np.maximum(n - 3, 1) + 1 / np.maximum(w.sum(0), 1))
        p_diff = 2 * ss.norm.sf(np.abs(np.arctanh(np.clip(r, -0.999, 0.999)) - zr) / se)
        hit = keep & (q < stats.FDR_THRESHOLD) & (np.abs(r) > stats.R_THRESHOLD) & (p_diff < 0.01)
        for j, i in zip(*np.where(hit)):
            s_ = np.sign(r[j, i])
            rr = r_rest[j, i] * s_
            out.append(dict(dataset=dname, lineage=L, n=int(n[j, i]), gene_a=rwe[j], gene_b=fam[i], r=r[j, i], p=P[L][j, i], FDR=q[j, i], r_other_lineages=r_rest[j, i],
                            p_difference=p_diff[j, i], restriction='restricted' if rr <= 0.3 else ('enriched' if rr <= 0.5 else 'shared, stronger here'),
                            **{f'r_{o}': Rs[o][j, i] for o in big}))
    print(f'{dname}: lineages tested {big} | (pair x lineage) tests {int(np.isfinite(flat).sum()):,} | hits {sum(1 for o in out if o["dataset"] == dname)}')
H = pd.DataFrame(out)
H['symbol_a'] = H.gene_a.map(data.symbol_of)
H['symbol_b'] = H.gene_b.map(data.symbol_of)
H = H[['dataset', 'lineage', 'n', 'symbol_a', 'symbol_b', 'r', 'p', 'FDR', 'r_other_lineages', 'p_difference', 'restriction', 'gene_a', 'gene_b'] +
      [c for c in H.columns if c.startswith('r_') and c != 'r_other_lineages']].sort_values(['dataset', 'p'])
save(H, 'tables', 'lineage_restricted_pairs.tsv', fmt='%.4g')
print(H.groupby(['dataset', 'lineage']).size().to_string())
print(H.groupby('dataset').head(12)[['dataset', 'lineage', 'n', 'symbol_a', 'symbol_b', 'r', 'FDR', 'r_other_lineages', 'restriction']].to_string(index=False, float_format=lambda v: f'{v:.3g}'))
