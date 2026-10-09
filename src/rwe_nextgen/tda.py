"""Gene dependency classes (Pan-dependency / Strongly Selective / High Variance / Weakly Selective / Non-dependency).

Re-implementation of the classification in the NextGen paper's code (src/generate.py: tda_classify_common_essentials,
tda_classify_strongly_selective, tda_classify_high_variance, partition_genes_by_tda_results with min_n_dep = 3 and the
priority order Non-dependency < Weakly Selective < High Variance < Strongly Selective < Pan-dependency).  Used to classify the
DepMap 2D comparator screens with exactly the rules the paper applied to NextGen; `pipeline/02_dependency_classes.py` first
checks that it reproduces the paper's published NextGen classes.
"""
import numpy as np
import pandas as pd
from scipy import stats


def common_essentials(ge, quantile=0.9):
    rank = ge.rank(axis=1, method='min', na_option='bottom')
    rank = rank.div(ge.notna().sum(axis=1), axis=0)
    pct = rank.quantile(q=quantile)
    kde = stats.gaussian_kde(pct.dropna(), bw_method=0.3)
    pts = np.arange(0.1, 0.9, 0.001)
    thr = pts[np.argmin(kde(pts))]
    return pd.DataFrame({'CE_percentile': pct, 'Common_Essential': pct <= thr}), float(thr)


def strongly_selective(ge, dep, prob=0.5, sel=-0.86):
    skew, kurt = ge.skew(axis=0), ge.kurtosis(axis=0)
    out = pd.DataFrame({'Mean': ge.mean(), 'Variance': ge.var(), 'Skewness': skew, 'Kurtosis': kurt})
    out['dep_lines'] = (dep > prob).sum()
    out['lines_with_data'] = dep.notna().sum()
    out['skew_x_kurtosis'] = skew * kurt
    return out


def high_variance(dep, expr):
    expr = expr.dropna(how='all')            # screens without expression data are dropped (inner merge in the paper's code)
    models = sorted(set(dep.index) & set(expr.index))
    genes = sorted(set(dep.columns) & set(expr.columns))
    d = dep.loc[models]
    e = expr.loc[models, genes]
    vals = e.values.ravel()
    vals = vals[~np.isnan(vals)]
    kde = stats.gaussian_kde(vals if len(vals) <= 3_000_000 else np.random.default_rng(1).choice(vals, 3_000_000, replace=False))
    pts = np.arange(1, 4, 0.005)
    cut = pts[np.argmin(kde(pts))]
    nonexp = (e < cut).mean()
    nonexp = set(nonexp[nonexp == 1].index)
    v = pd.DataFrame({'prob_variance': d.var()})
    v['unexpressed'] = v.index.isin(nonexp)
    thr = v.loc[v.unexpressed, 'prob_variance'].quantile(0.99)
    v['high_variance'] = (v.prob_variance > thr) & ~v.unexpressed
    return v, float(cut)


def classify(ge, dep, expr_screens, min_n_dep=3):
    """ge, dep: screens x genes for ONE screen set; expr_screens: screens x genes expression of the same screens"""
    ce, ce_thr = common_essentials(ge)
    ss = strongly_selective(ge, dep)
    hv, cut = high_variance(dep, expr_screens)
    res = ce.join(ss).join(hv)
    unscreened = res.index[res.lines_with_data == 0]
    calls = res.drop(unscreened)
    lab = pd.Series('Unscreened', index=res.index)
    order = [('Non-dependency', calls.index[calls.dep_lines == 0]),
             ('Weakly Selective', calls.index[calls.dep_lines > 0]),
             ('High Variance', calls.index[calls.high_variance.fillna(False).astype(bool)]),
             ('Strongly Selective', calls.index[(calls.skew_x_kurtosis < -0.86) & (calls.dep_lines >= min_n_dep)]),
             ('Pan-dependency', calls.index[calls.Common_Essential.fillna(False).astype(bool)])]
    for name, genes in order:
        lab.loc[genes] = name
    out = pd.DataFrame({'Class': lab, 'Variance': res.Variance, 'Selectivity': res.skew_x_kurtosis, 'DepLines': res.dep_lines, '90%ileLineRank': res.CE_percentile})
    return out, dict(common_essential_rank_threshold=ce_thr, expression_cutoff=cut)
