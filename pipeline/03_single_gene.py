"""Step 3 - single-gene dependency profile of every RWE gene, NextGen next to DepMap 2D, with the paper's conventions.

  classes      paper NextGen / organoid classes; DepMap-2D class with the same rules (step 2)
  effect       mean gene effect and fraction of dependent screens (probability > 0.5) in NextGen, 3D organoids, CNS neurospheres, DepMap 2D
  organoid     3D organoids vs lineage-matched adherent 2D: the paper's table (Mann-Whitney U, BH over all genes; DropEssential flag)
  lineage      each lineage vs the rest (Mann-Whitney U, two-sided; BH over the RWE x lineage table of each dataset); hit = FDR < 0.05,
               mean difference < -0.3, lineage n >= 5 (thresholds of the paper's alteration / lineage analyses)
  biomarkers   expression addiction and paralog dependency: the paper's tables, FDR and calls exactly as in the paper (FDR < 0.005, |r| > 0.3)
  library      24 models screened with both Cas9 and Cas12a (Humagne-CD, the NextGen library): paired mean difference per gene
  sets         Fisher tests of dependency classes for RWE families / roles vs the genome (BH)
Outputs: results/tables/single_gene_rwe.tsv, lineage_tests_rwe.tsv, biomarkers_rwe.tsv, set_enrichment.tsv
"""
import os

import numpy as np
import pandas as pd
from _common import banner, data, paths, save, stats
from scipy import stats as ss

banner('3. single-gene dependency profiles (NextGen and DepMap 2D)')
ge, dep, sm = data.screen_effect(), data.screen_dependency(), data.screen_meta()
cls = pd.read_csv(os.path.join(paths.TABLES, 'dependency_classes.tsv'), sep='\t').set_index('key')
u = data.rwe_genes().dropna(subset=['key']).drop_duplicates('key')
keys = [k for k in u.key if k in ge.columns]
NG, TR = list(data.nextgen_screens()), list(data.traditional_lineage_matched_screens())
ORG = list(data.organoid_3d_screens())
CNS = [s for s in NG if sm.loc[s, 'ScreenType'] in ('3DN', '2DN')]
groups = {'NextGen': NG, 'Organoid3D': ORG, 'CNS': CNS, 'DepMap2D': TR}
T = u.set_index('key').loc[keys, ['symbol', 'family', 'role', 'evidence_tier', 'caveat_class', 'modifications', 'in_lab_list']].copy()
T['NextGenClass'] = cls.NextGenClass.reindex(keys).values
T['OrganoidsClass'] = cls.OrganoidsClass.reindex(keys).values
T['DepMap2DClass'] = cls.DepMap2DClass.reindex(keys).values
for g, scr in groups.items():
    T[f'{g}_n'] = len(scr)
    T[f'{g}_mean'] = ge.loc[scr, keys].mean().values
    T[f'{g}_sd'] = ge.loc[scr, keys].std().values
    d = dep.reindex(scr)[keys]
    T[f'{g}_frac_dep'] = ((d > 0.5).sum() / d.notna().sum()).values
T['NextGen_minus_2D_mean'] = T.NextGen_mean - T.DepMap2D_mean

# organoid vs adherent 2D: paper table
o2 = pd.read_csv(os.path.join(paths.FIGSHARE, 'organoids_vs_2d_gene_dependency_tests_full.csv')).set_index('Dependency')
for c in ('OrganoidsMean', 'AdherentMean', 'MeanDifference', 'MannWhitneyPValue', 'MannWhitneyFDR', 'DropEssential'):
    T[f'org_vs_2D_{c}'] = o2[c].reindex(keys).values
T['org_vs_2D_call'] = np.where((T.org_vs_2D_MannWhitneyFDR < 0.05) & (T.org_vs_2D_MeanDifference.abs() > 0.3) & ~T.org_vs_2D_DropEssential.fillna(False).astype(bool),
                               np.where(T.org_vs_2D_MeanDifference < 0, 'more dependent in organoids', 'less dependent in organoids'), '')

# library check: 24 models screened with both libraries (2D)
c12 = sm[sm.IsTraditional2D & sm.Library.eq('Humagne-CD')].index.intersection(ge.index)
c9 = sm[sm.IsTraditional2D & ~sm.Library.eq('Humagne-CD')].index.intersection(ge.index)
dual = sorted(set(sm.loc[c12, 'ModelID']) & set(sm.loc[c9, 'ModelID']))
E12 = ge.loc[c12, keys].groupby(sm.loc[c12, 'ModelID'].values).mean().reindex(dual)
E9 = ge.loc[c9, keys].groupby(sm.loc[c9, 'ModelID'].values).mean().reindex(dual)
dd = E12 - E9
T['library_dual_n'] = len(dual)
T['library_Cas12a_minus_Cas9'] = dd.mean().values
T['library_paired_p'] = ss.ttest_1samp(dd, 0, nan_policy='omit').pvalue
T['library_paired_FDR'] = stats.bh(T.library_paired_p.values)
T['library_sensitive'] = (T.library_paired_FDR < 0.05) & (T.library_Cas12a_minus_Cas9.abs() > 0.3)
print(f'library check on {len(dual)} dual-library models: library-sensitive RWE genes (FDR < 0.05, |paired difference| > 0.3): {int(T.library_sensitive.sum())}',
      T[T.library_sensitive].symbol.tolist())

# lineage vs rest, both datasets
rows = []
for dname, scr in (('NextGen', NG), ('DepMap2D', TR)):
    lin = sm.loc[scr, 'OncotreeLineage']
    X = ge.loc[scr, keys]
    for L, n in lin.value_counts().items():
        if n < 5:
            continue
        a, b = X[lin.values == L], X[lin.values != L]
        p = ss.mannwhitneyu(a, b, axis=0, nan_policy='omit', alternative='two-sided').pvalue
        for k, pm, ma, mb in zip(keys, np.atleast_1d(p), a.mean().values, b.mean().values):
            rows.append(dict(dataset=dname, lineage=L, key=k, symbol=data.symbol_of(k), ingroup_n=int(n), outgroup_n=len(scr) - int(n), ingroup_mean=ma, outgroup_mean=mb,
                             mean_difference=ma - mb, p=pm))
LT = pd.DataFrame(rows)
LT['FDR'] = np.nan
for dname in LT.dataset.unique():
    m = LT.dataset == dname
    LT.loc[m, 'FDR'] = stats.bh(LT.loc[m, 'p'].values)
LT['hit'] = (LT.FDR < 0.05) & (LT.mean_difference < -0.3)
w = LT.pivot_table(index=['key', 'lineage'], columns='dataset', values=['mean_difference', 'FDR', 'hit'], aggfunc='first')
w.columns = [f'{b}_{a}' for a, b in w.columns]
w = w.reset_index()
w['symbol'] = w.key.map(data.symbol_of)
w['class'] = np.select([w.NextGen_hit.fillna(False).astype(bool) & w.DepMap2D_hit.fillna(False).astype(bool), w.NextGen_hit.fillna(False).astype(bool),
                        w.DepMap2D_hit.fillna(False).astype(bool)], ['both', 'NextGen only', 'DepMap 2D only'], '')
print(f'lineage-selective RWE dependencies (FDR < 0.05, mean difference < -0.3): NextGen {int(LT[(LT.dataset == "NextGen")].hit.sum())}, DepMap 2D {int(LT[(LT.dataset == "DepMap2D")].hit.sum())}')
print(w[w['class'] != ''].sort_values(['class', 'NextGen_mean_difference'])[['symbol', 'lineage', 'NextGen_mean_difference', 'NextGen_FDR', 'DepMap2D_mean_difference', 'DepMap2D_FDR', 'class']]
      .to_string(index=False, float_format=lambda v: f'{v:.3g}'))
best = LT[LT.dataset == 'NextGen'].sort_values('mean_difference').groupby('key').head(1).set_index('key')
T['NextGen_top_lineage'] = best.lineage.reindex(keys).values
T['NextGen_top_lineage_diff'] = best.mean_difference.reindex(keys).values
T['NextGen_top_lineage_FDR'] = best.FDR.reindex(keys).values
best2 = LT[LT.dataset == 'DepMap2D'].sort_values('mean_difference').groupby('key').head(1).set_index('key')
T['DepMap2D_top_lineage'] = best2.lineage.reindex(keys).values
T['DepMap2D_top_lineage_diff'] = best2.mean_difference.reindex(keys).values
T['DepMap2D_top_lineage_FDR'] = best2.FDR.reindex(keys).values

# biomarkers: the paper's expression-addiction and paralog tables, FDR and calls as in the paper
DEP = ['Strongly Selective', 'High Variance', 'Pan-dependency']
ea = pd.read_csv(os.path.join(paths.FIGSHARE, 'expression_addiction_tests_full.csv'))
ea = ea[ea.NextGenDepClass.isin(DEP)].dropna(subset=['NextGenPValue']).copy()
ea['FDR'] = stats.bh(ea.NextGenPValue.values)
ea['Significant'] = (ea.FDR < 0.005) & (ea.NextGenPearsonR < -0.3)
ea['type'] = 'expression addiction (own expression)'
pa = pd.read_csv(os.path.join(paths.FIGSHARE, 'paralog_dependency_tests_full.csv'))
pa = pa[pa.NextGenDepClass.isin(DEP) & (pa.BiomarkerVariance > 1)].dropna(subset=['NextGenPValue']).copy()
pa['FDR'] = stats.bh(pa.NextGenPValue.values)
pa['Significant'] = (pa.NextGenPearsonR > 0.3) & (pa.FDR < 0.005)
pa['type'] = 'paralog expression'
BM = pd.concat([ea, pa], ignore_index=True)
BM = BM[BM.Dependency.isin(keys)].copy()
BM['symbol'] = BM.Dependency.map(data.symbol_of)
BM['biomarker_symbol'] = BM.Biomarker.map(data.symbol_of)
print(f'biomarker calls among RWE genes (paper tables / thresholds): expression addiction {int(BM[(BM.type.str.startswith("expression")) & BM.Significant].shape[0])},'
      f' paralog {int(BM[(BM.type.str.startswith("paralog")) & BM.Significant].shape[0])}')
print(BM[BM.Significant][['symbol', 'type', 'biomarker_symbol', 'NextGenPearsonR', 'FDR', 'TraditionalPearsonR']].to_string(index=False, float_format=lambda v: f'{v:.3g}'))
sig = BM[BM.Significant].groupby('Dependency').apply(lambda d: '; '.join(f'{t.split(" ")[0]}:{b} r={r:+.2f}' for t, b, r in zip(d.type, d.biomarker_symbol, d.NextGenPearsonR)))
T['biomarkers'] = sig.reindex(keys).fillna('').values

# set enrichment (Fisher) of dependency classes vs the genome
rows = []
screened = cls.index[cls.NextGenClass.notna() & cls.NextGenClass.ne('Unscreened')]
ctrl = set(data.control_genes().key.dropna())
sets = {'RWE: all': set(T.index)}
for col in ('family', 'role', 'caveat_class'):
    for v in T[col].fillna('none').unique():
        sets[f'{col}: {v}'] = set(T.index[T[col].fillna('none') == v])
sets['control genes (non-RNA MTases / IGFBPs / ALKBH)'] = ctrl
for dname, col in (('NextGen', 'NextGenClass'), ('DepMap2D', 'DepMap2DClass')):
    c_ = cls[col].reindex(screened)
    for outcome, flag in (('Pan-dependency', c_.eq('Pan-dependency')), ('any dependency class', c_.isin(DEP))):
        for name, ks in sets.items():
            ks = [k for k in ks if k in flag.index]
            if len(ks) < 3:
                continue
            a = int(flag.loc[ks].sum())
            rest = flag.drop(ks)
            tab = [[a, len(ks) - a], [int(rest.sum()), len(rest) - int(rest.sum())]]
            orr, p = ss.fisher_exact(tab)
            rows.append(dict(dataset=dname, outcome=outcome, set=name, n=len(ks), n_hit=a, frac=a / len(ks), genome_frac=rest.mean(), odds_ratio=orr, p=p))
SE = pd.DataFrame(rows)
SE['FDR'] = np.nan
for dname in SE.dataset.unique():
    m = SE.dataset == dname
    SE.loc[m, 'FDR'] = stats.bh(SE.loc[m, 'p'].values)
print(SE[SE.set.isin(['RWE: all', 'role: writer', 'role: writer-cofactor', 'role: reader', 'role: eraser', 'control genes (non-RNA MTases / IGFBPs / ALKBH)'])
         & SE.outcome.eq('Pan-dependency')].to_string(index=False, float_format=lambda v: f'{v:.3g}'))

T = T.reset_index()
save(T, 'tables', 'single_gene_rwe.tsv')
save(w, 'tables', 'lineage_tests_rwe.tsv')
save(BM, 'tables', 'biomarkers_rwe.tsv')
save(SE, 'tables', 'set_enrichment.tsv')
print('\nNextGen class vs DepMap-2D class (RWE genes):')
print(pd.crosstab(T.NextGenClass, T.DepMap2DClass).to_string())
print('\norganoid vs adherent 2D calls (paper table; FDR < 0.05, |diff| > 0.3, not common-essential):', T.org_vs_2D_call.value_counts().to_dict())
print(T[T.org_vs_2D_call != ''].sort_values('org_vs_2D_MeanDifference')[['symbol', 'org_vs_2D_MeanDifference', 'org_vs_2D_MannWhitneyFDR', 'library_Cas12a_minus_Cas9', 'library_sensitive']]
      .to_string(index=False, float_format=lambda v: f'{v:.3g}'))
