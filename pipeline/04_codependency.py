"""Step 4 - genome-wide co-dependency among dependency genes, NextGen and DepMap 2D, with the paper's statistics.

Family   : genes that are a dependency (Strongly Selective / High Variance / Pan-dependency) in NextGen OR in DepMap 2D (step 2);
           the same pairs are tested in both datasets so that the two maps can be compared pair by pair.
Test     : Pearson r on pairwise-complete screens (no covariates), two-sided p.
FDR      : one Benjamini-Hochberg correction over all unique gene pairs of the family, separately per dataset.
Call     : FDR < 0.005 and |r| > 0.3 (the paper's thresholds for correlation-type associations).
Library  : the DepMap 2D comparator mixes Avana, KY and Humagne screens (as in the paper); genes whose within-lineage 2D effect is
           largely a library difference (R^2 of KY vs Avana >= 0.1) are flagged, and so are their 2D calls.

Outputs
  results/intermediate/codep_<dataset>.npz        r (float16), packed call matrix (upper triangle), gene order  (re-creatable, not versioned)
  results/intermediate/codep_rwe_pairs.tsv.gz     every pair with at least one RWE gene: r, p, FDR, n and call in both datasets
  results/tables/codep_overview.tsv               genome-wide numbers per dataset (pairs, calls, p at the FDR boundary)
  results/tables/codep_rwe_calls.tsv              pairs with an RWE gene called in either dataset, with comparison class
  results/tables/codep_summary_per_gene.tsv       per RWE gene: number of calls in each dataset and the strongest partners
  results/intermediate/codep_neighbours.npz       top partners of every family gene (for the explorer)
"""
import os

import numpy as np
import pandas as pd
from _common import banner, data, paths, save, stats

TOP_K = 30          # partners kept per gene for the explorer

banner('4. co-dependency (paper statistics), NextGen and DepMap 2D')
cls = pd.read_csv(os.path.join(paths.TABLES, 'dependency_classes.tsv'), sep='\t').set_index('key')
ge = data.screen_effect()
fam = [k for k in cls.index[cls.in_codependency_family] if k in ge.columns]
rwe_keys = set(data.rwe_genes().key.dropna())
F = len(fam)
is_rwe = np.array([k in rwe_keys for k in fam])
print(f'family: {F:,} dependency genes (RWE among them: {int(is_rwe.sum())}) -> {F * (F - 1) // 2:,} unique pairs per dataset')


iu = np.triu_indices(F, 1)
res = {}
for dname, screens in data.DATASETS.items():
    X = ge.loc[list(screens()), fam].values.astype(np.float64)
    R, N = stats.pairwise_r(X)
    np.fill_diagonal(R, np.nan)
    r, n = R[iu].astype(np.float64), N[iu].astype(np.float64)
    p = stats.r_to_p(r, n)
    q = stats.bh(p)
    call = (q < stats.FDR_THRESHOLD) & (np.abs(r) > stats.R_THRESHOLD)
    pmax = np.nanmax(np.where(q < stats.FDR_THRESHOLD, p, np.nan))
    print(f'{dname}: {X.shape[0]} screens | calls {int(call.sum()):,} of {len(r):,} pairs ({call.mean():.2%}) | p at the FDR boundary {pmax:.1e} | median |r| of calls {np.median(np.abs(r[call])):.2f}')
    res[dname] = dict(R=R, N=N, r=r, p=p, q=q, call=call, overview=dict(dataset=dname, screens=X.shape[0], family_genes=F, pairs=len(r), calls=int(call.sum()),
                                                                     p_boundary=float(pmax), median_abs_r_calls=float(np.median(np.abs(r[call])))))

save(pd.DataFrame([d['overview'] for d in res.values()]), 'tables', 'codep_overview.tsv')

# ------------------------------------------------------------------------------- 2D library sensitivity per gene (within-lineage KY vs Avana)
sm = data.screen_meta()
tr = list(data.traditional_lineage_matched_screens())
lib, lin = sm.loc[tr, 'Library'], sm.loc[tr, 'OncotreeLineage']
m = lib.isin(['Avana', 'KY']).values
Xd = ge.loc[tr, fam]
Xd = (Xd - Xd.groupby(lin.values).transform('mean')).values[m]
ky = (lib.values[m] == 'KY').astype(float)
ky = (ky - ky.mean()) / ky.std()
Xc = np.where(np.isnan(Xd), 0, Xd - np.nanmean(Xd, 0))
lib_r2 = ((Xc * ky[:, None]).sum(0) / np.sqrt((Xc ** 2).sum(0) * (ky ** 2).sum())) ** 2
lib_flag = lib_r2 >= 0.1
print(f'DepMap 2D library-sensitive genes (R^2 KY vs Avana >= 0.1): {int(lib_flag.sum())}; RWE among them:',
      [data.symbol_of(k) for k, f in zip(fam, lib_flag) if f and k in rwe_keys])
save(pd.DataFrame({'key': fam, 'DepMap2D_library_R2': lib_r2, 'DepMap2D_library_sensitive': lib_flag}), 'tables', 'library_sensitivity_2D.tsv')

# ------------------------------------------------------------------------------- pairs with an RWE gene (both datasets)
a_idx, b_idx = iu
sel = is_rwe[a_idx] | is_rwe[b_idx]
T = pd.DataFrame({'gene_a': np.array(fam)[a_idx[sel]], 'gene_b': np.array(fam)[b_idx[sel]]})
for dname, d in res.items():
    T[f'r_{dname}'] = d['r'][sel]
    T[f'p_{dname}'] = d['p'][sel]
    T[f'FDR_{dname}'] = d['q'][sel]
    T[f'n_{dname}'] = d['N'][iu][sel]
    T[f'call_{dname}'] = d['call'][sel]
T['DepMap2D_library_flag'] = lib_flag[a_idx[sel]] | lib_flag[b_idx[sel]]
T['p_difference'] = stats.fisher_diff_p(T.r_NextGen, T.n_NextGen, T.r_DepMap2D, T.n_DepMap2D)
T['class'] = np.select([T.call_NextGen & T.call_DepMap2D, T.call_NextGen, T.call_DepMap2D], ['both', 'NextGen only', 'DepMap 2D only'], '')
for c in ('a', 'b'):
    T[f'symbol_{c}'] = T[f'gene_{c}'].map(data.symbol_of)
    T[f'rwe_{c}'] = T[f'gene_{c}'].isin(rwe_keys)
T['rwe_rwe'] = T.rwe_a & T.rwe_b
cols = ['symbol_a', 'symbol_b', 'gene_a', 'gene_b', 'rwe_a', 'rwe_b', 'rwe_rwe', 'class'] + [c for c in T.columns if c.endswith(('_NextGen', '_DepMap2D'))] + ['DepMap2D_library_flag', 'p_difference']
T = T[cols]
save(T, 'intermediate', 'codep_rwe_pairs.tsv.gz', fmt='%.4g')
C = T[T['class'] != ''].sort_values('FDR_NextGen')
save(C, 'tables', 'codep_rwe_calls.tsv', fmt='%.4g')
print('RWE pairs called:', C['class'].value_counts().to_dict(), '| RWE-RWE:', C[C.rwe_rwe]['class'].value_counts().to_dict())

# ------------------------------------------------------------------------------- per RWE gene summary
rows = []
gi = {k: i for i, k in enumerate(fam)}
for k in sorted(rwe_keys):
    row = dict(key=k, symbol=data.symbol_of(k), in_family=k in gi)
    if k in gi:
        i = gi[k]
        for dname, d in res.items():
            rr = d['R'][i].astype(float)
            sub = C[(C.gene_a == k) | (C.gene_b == k)]
            called = sub[sub[f'call_{dname}']]
            row[f'n_calls_{dname}'] = len(called)
            order = np.argsort(-np.nan_to_num(np.abs(rr), nan=-1))[:8]
            row[f'top_partners_{dname}'] = ', '.join(f'{data.symbol_of(fam[j])} {rr[j]:+.2f}' for j in order)
        row['n_calls_both'] = int((((C.gene_a == k) | (C.gene_b == k)) & (C['class'] == 'both')).sum())
    rows.append(row)
S = pd.DataFrame(rows)
for c in ('n_calls_NextGen', 'n_calls_DepMap2D', 'n_calls_both'):
    S[c] = S[c].astype('Int64')
save(S, 'tables', 'codep_summary_per_gene.tsv')

# ------------------------------------------------------------------------------- explorer neighbour lists and compact matrices
nb = {}
for dname, d in res.items():
    R = d['R']
    sig = np.zeros_like(R, dtype=bool)
    sig[iu] = d['call']
    sig = sig | sig.T
    A = np.where(sig, np.abs(np.nan_to_num(R)), -1.0)
    top = np.argsort(-A, axis=1)[:, :TOP_K]
    keep = np.take_along_axis(A, top, axis=1) > 0
    nb[dname] = np.where(keep, top, -1).astype(np.int32)
    nb[dname + '_degree'] = sig.sum(1).astype(np.int32)
os.makedirs(paths.INTERMEDIATE, exist_ok=True)
np.savez_compressed(os.path.join(paths.INTERMEDIATE, 'codep_neighbours.npz'), genes=np.array(fam), **nb)
for dname, d in res.items():
    call = np.zeros((F, F), bool)
    call[iu] = d['call']
    np.savez_compressed(os.path.join(paths.INTERMEDIATE, f'codep_{dname}.npz'), genes=np.array(fam), r=d['R'].astype(np.float16),
                        call_upper=np.packbits(call, axis=None), n_genes=F, n_screens=len(data.DATASETS[dname]()))
print('degree (number of calls) per gene - median NextGen %d, DepMap2D %d' % (np.median(nb['NextGen_degree']), np.median(nb['DepMap2D_degree'])))

# ------------------------------------------------------------------------------- headline numbers for RWE genes of interest
show = ['METTL3', 'METTL14', 'WTAP', 'VIRMA', 'YTHDC1', 'YTHDF2', 'PUS7', 'DKC1', 'TRUB1', 'METTL1', 'NSUN2', 'TRMT6']
print(S[S.symbol.isin(show)][['symbol', 'in_family', 'n_calls_NextGen', 'n_calls_DepMap2D', 'top_partners_NextGen']].to_string(index=False))
