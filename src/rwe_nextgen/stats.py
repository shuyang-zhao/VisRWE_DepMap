"""Statistics used throughout, matching the paper's code: Pearson correlation on pairwise-complete screens (scipy.stats.pearsonr
semantics: two-sided p from a t distribution with n - 2 df) and one Benjamini-Hochberg FDR per results table (scipy.stats.false_discovery_control)."""
import numpy as np
from scipy import stats

FDR_THRESHOLD = 0.005      # paper: correlation-type biomarker analyses (expression addiction, paralogs)
R_THRESHOLD = 0.3


def pearson_vs_all(x, Y):
    """Pearson r, two-sided p and n of vector x against every column of matrix Y, pairwise-complete (NaN = missing)."""
    x = np.asarray(x, float)
    Y = np.asarray(Y, float)
    M = ~np.isnan(Y) & ~np.isnan(x)[:, None]
    n = M.sum(0).astype(float)
    xm = np.where(M, x[:, None], 0.0)
    ym = np.where(M, Y, 0.0)
    sx, sy = xm.sum(0), ym.sum(0)
    sxx, syy, sxy = (xm * xm).sum(0), (ym * ym).sum(0), (xm * ym).sum(0)
    with np.errstate(invalid='ignore', divide='ignore'):
        r = (n * sxy - sx * sy) / np.sqrt((n * sxx - sx ** 2) * (n * syy - sy ** 2))
        r = np.clip(r, -1.0, 1.0)
        p = r_to_p(r, n)
    r[n < 3] = np.nan
    p[n < 3] = np.nan
    return r, p, n.astype(int)


def pairwise_r(X, rows=None, block=400):
    """Pearson r and n between columns `rows` (default: all) and all columns of X (screens x genes), pairwise-complete.
    Returns float32 r (len(rows) x n_genes) and int16 n; computed in blocks with six matrix products each."""
    X = np.asarray(X, float)
    rows = np.arange(X.shape[1]) if rows is None else np.asarray(rows)
    M = (~np.isnan(X)).astype(np.float64)
    Z = np.where(M > 0, X, 0.0)
    Z2 = Z * Z
    R = np.empty((len(rows), X.shape[1]), np.float32)
    N = np.empty((len(rows), X.shape[1]), np.int16)
    for s in range(0, len(rows), block):
        b = rows[s:s + block]
        n = M[:, b].T @ M
        sx, sy = Z[:, b].T @ M, M[:, b].T @ Z
        sxx, syy = Z2[:, b].T @ M, M[:, b].T @ Z2
        sxy = Z[:, b].T @ Z
        with np.errstate(invalid='ignore', divide='ignore'):
            r = (n * sxy - sx * sy) / np.sqrt((n * sxx - sx ** 2) * (n * syy - sy ** 2))
        r[n < 3] = np.nan
        R[s:s + block] = np.clip(r, -1, 1)
        N[s:s + block] = n
    return R, N


def r_to_p(r, n):
    r = np.asarray(r, float)
    n = np.asarray(n, float)
    with np.errstate(invalid='ignore', divide='ignore'):
        t = r * np.sqrt((n - 2) / (1 - r ** 2))
        return 2 * stats.t.sf(np.abs(t), n - 2)


def bh(p):
    """Benjamini-Hochberg FDR (scipy.stats.false_discovery_control); NaN stay NaN."""
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        q[ok] = stats.false_discovery_control(p[ok])
    return q


def fisher_ci(r, n, level=0.95):
    z = np.arctanh(np.clip(r, -0.999999, 0.999999))
    se = 1 / np.sqrt(np.maximum(np.asarray(n, float) - 3, 1))
    k = stats.norm.isf((1 - level) / 2)
    return np.tanh(z - k * se), np.tanh(z + k * se)


def fisher_diff_p(r1, n1, r2, n2):
    """two-sided P that two independent correlations differ (Fisher z)"""
    z1, z2 = np.arctanh(np.clip(r1, -0.999999, 0.999999)), np.arctanh(np.clip(r2, -0.999999, 0.999999))
    se = np.sqrt(1 / np.maximum(np.asarray(n1, float) - 3, 1) + 1 / np.maximum(np.asarray(n2, float) - 3, 1))
    return 2 * stats.norm.sf(np.abs(z1 - z2) / se)


def within_group_r(x, y, groups):
    """Pearson r after removing the group (lineage) means from both vectors; pairwise-complete."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    g = np.asarray(groups)
    m = ~np.isnan(x) & ~np.isnan(y)
    x, y, g = x[m], y[m], g[m]
    xd, yd = x.copy(), y.copy()
    for G in np.unique(g):
        k = g == G
        xd[k] -= x[k].mean()
        yd[k] -= y[k].mean()
    den = np.sqrt((xd ** 2).sum() * (yd ** 2).sum())
    return float((xd * yd).sum() / den) if den > 0 else np.nan
