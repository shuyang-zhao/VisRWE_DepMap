"""Step 9 - data bundle for the interactive explorer (docs/index.html).

Genes   : every dependency gene of the co-dependency family plus the RWE genes outside it, with classes, categories and summary statistics.
Edges   : all called pairs with an RWE gene, plus the 20 strongest called partners of every family gene in each dataset (so that any node
          can be expanded); r in both datasets and call flags (bit 1 = NextGen, bit 2 = DepMap 2D, bit 4 = 2D library-sensitive gene).
Vectors : gene effect of every gene in every NextGen / DepMap-2D screen (int8, 0.025 resolution, gzip + base64) for the per-edge
          scatter plots and lineage decomposition computed in the browser; plus the traditional 2D screens made with the Cas12a library
          (library check: is a NextGen-only pair also present in 2D lines screened with the NextGen library?).
Layout  : coordinates of the RWE overview graph (union of both datasets), computed once here.
Output  : results/intermediate/explorer_data.json
"""
import base64
import datetime
import gzip
import json
import os

import igraph as ig
import numpy as np
import pandas as pd
from _common import banner, data, paths

TOP = 20
banner('9. explorer data bundle')
T = lambda *p: pd.read_csv(os.path.join(paths.TABLES, *p), sep='\t')  # noqa: E731
cls = T('dependency_classes.tsv').set_index('key')
cat = T('gene_categories.tsv').set_index('key')
sg = T('single_gene_rwe.tsv').set_index('key')
lib = T('library_sensitivity_2D.tsv').set_index('key')
lit = T('literature_rwe.tsv').set_index('key')
mcat, mpw = T('module_enrichment_categories.tsv'), T('module_enrichment_pathways.tsv')
calls = T('codep_rwe_calls.tsv')
Lr = T('lineage_restricted_pairs.tsv')
ge, sm = data.screen_effect(), data.screen_meta()
rwe = data.rwe_genes().dropna(subset=['key']).drop_duplicates('key').set_index('key')

Z = {d: np.load(os.path.join(paths.INTERMEDIATE, f'codep_{d}.npz'), allow_pickle=True) for d in data.DATASETS}
fam = [str(g) for g in Z['NextGen']['genes']]
F = len(fam)
genes = fam + [k for k in rwe.index if k not in set(fam) and k in ge.columns]
gi = {k: i for i, k in enumerate(genes)}
R = {d: Z[d]['r'] for d in Z}
CALL = {d: np.unpackbits(Z[d]['call_upper'])[:F * F].reshape(F, F).astype(bool) for d in Z}
for d in CALL:
    CALL[d] = CALL[d] | CALL[d].T
print(f'genes in the bundle: {len(genes):,} (family {F:,} + RWE outside the family {len(genes) - F})')

# ------------------------------------------------------------------------------------------------ edges
nb = np.load(os.path.join(paths.INTERMEDIATE, 'codep_neighbours.npz'))
pairs = set()
for a, b in zip(calls.gene_a, calls.gene_b):
    i, j = gi[a], gi[b]
    pairs.add((min(i, j), max(i, j)))
for d in data.DATASETS:
    top = nb[d][:, :TOP]
    for i in range(F):
        for j in top[i]:
            if j >= 0:
                pairs.add((min(i, j), max(i, j)))
for x, y in zip(Lr.gene_a, Lr.gene_b):                                                  # lineage-restricted pairs (step 5)
    i, j = gi[x], gi[y]
    pairs.add((min(i, j), max(i, j)))
pairs = np.array(sorted(pairs), dtype=np.int32)
a, b = pairs[:, 0], pairs[:, 1]
lib_flag = lib.DepMap2D_library_sensitive.reindex(genes).fillna(False).values.astype(bool)
flags = (CALL['NextGen'][a, b].astype(np.uint8) | (CALL['DepMap2D'][a, b].astype(np.uint8) << 1) | ((lib_flag[a] | lib_flag[b]).astype(np.uint8) << 2))
r_ng = np.nan_to_num(R['NextGen'][a, b].astype(float))
r_2d = np.nan_to_num(R['DepMap2D'][a, b].astype(float))
print(f'edges: {len(pairs):,} | called in NextGen {int((flags & 1).astype(bool).sum()):,}, DepMap 2D {int((flags & 2).astype(bool).sum()):,}')


def b64(arr):
    return base64.b64encode(np.ascontiguousarray(arr).tobytes()).decode('ascii')


def b64gz(arr):
    return base64.b64encode(gzip.compress(np.ascontiguousarray(arr).tobytes(), 9)).decode('ascii')


q8 = lambda v: np.clip(np.round(np.asarray(v) * 100), -100, 100).astype(np.int8)  # noqa: E731

# ------------------------------------------------------------------------------------------------ gene attributes
CLASSES = ['Pan-dependency', 'Strongly Selective', 'High Variance', 'Weakly Selective', 'Non-dependency', 'Unscreened']
CATS = sorted(cat.category.unique())
FAMS = sorted(rwe.family.dropna().unique())
ROLES = sorted(rwe.role.dropna().unique())
ix = lambda lst, v: lst.index(v) if isinstance(v, str) and v in lst else -1  # noqa: E731
deg = {d: np.concatenate([nb[f'{d}_degree'], np.zeros(len(genes) - F, int)]) for d in data.DATASETS}
mean = {}
fdep = {}
dep = data.screen_dependency()
for d, scr in data.DATASETS.items():
    s_ = list(scr())
    mean[d] = ge.loc[s_].reindex(columns=genes).mean().values
    dd = dep.reindex(s_)[genes]
    fdep[d] = ((dd > 0.5).sum() / dd.notna().sum()).values
G = dict(
    sym=[data.symbol_of(k) for k in genes], key=genes,
    rwe=[int(k in rwe.index) for k in genes],
    fam=[ix(FAMS, rwe.family.get(k)) if k in rwe.index else -1 for k in genes],
    role=[ix(ROLES, rwe.role.get(k)) if k in rwe.index else -1 for k in genes],
    cav=[int(isinstance(rwe.caveat_class.get(k), str) and rwe.caveat_class.get(k) != '') if k in rwe.index else 0 for k in genes],
    cat=[ix(CATS, cat.category.get(k)) for k in genes],
    cng=[ix(CLASSES, cls.NextGenClass.get(k)) for k in genes], c2d=[ix(CLASSES, cls.DepMap2DClass.get(k)) for k in genes],
    mng=[None if np.isnan(v) else round(float(v), 3) for v in mean['NextGen']], m2d=[None if np.isnan(v) else round(float(v), 3) for v in mean['DepMap2D']],
    fng=[None if np.isnan(v) else round(float(v), 3) for v in fdep['NextGen']], f2d=[None if np.isnan(v) else round(float(v), 3) for v in fdep['DepMap2D']],
    dng=[int(v) for v in deg['NextGen']], d2d=[int(v) for v in deg['DepMap2D']], lib=[int(v) for v in lib_flag], infam=[int(i < F) for i in range(len(genes))])

# ------------------------------------------------------------------------------------------------ RWE gene cards
cards = {}
for k in rwe.index:
    if k not in gi:
        continue
    s = sg.loc[k] if k in sg.index else None
    li = lit.loc[k] if k in lit.index else None
    mm = mpw[mpw.key == k]
    top_mod = {d: [(r.module, int(getattr(r, f'overlap_{d}')), float(f"{getattr(r, f'FDR_{d}'):.2g}")) for r in mm[mm[f'call_{d}']].sort_values(f'p_{d}').head(5).itertuples()]
               for d in data.DATASETS}
    cards[gi[k]] = dict(
        mods=rwe.modifications.get(k), subs=rwe.substrates.get(k), tier=rwe.evidence_tier.get(k), caveat=rwe.caveat_class.get(k) if isinstance(rwe.caveat_class.get(k), str) else '',
        org=None if s is None or pd.isna(s.org_vs_2D_MeanDifference) else [round(float(s.org_vs_2D_MeanDifference), 3), float(f'{s.org_vs_2D_MannWhitneyFDR:.2g}'), bool(s.org_vs_2D_DropEssential), s.org_vs_2D_call if isinstance(s.org_vs_2D_call, str) else ''],
        lin_ng=None if s is None or pd.isna(s.NextGen_top_lineage_diff) else [s.NextGen_top_lineage, round(float(s.NextGen_top_lineage_diff), 3), float(f'{s.NextGen_top_lineage_FDR:.2g}')],
        lin_2d=None if s is None or pd.isna(s.DepMap2D_top_lineage_diff) else [s.DepMap2D_top_lineage, round(float(s.DepMap2D_top_lineage_diff), 3), float(f'{s.DepMap2D_top_lineage_FDR:.2g}')],
        libdiff=None if s is None or pd.isna(s.library_Cas12a_minus_Cas9) else round(float(s.library_Cas12a_minus_Cas9), 3),
        bio=s.biomarkers if s is not None and isinstance(s.biomarkers, str) else '',
        stance=li.claim_stance if li is not None else '', conc=li.concordance_NextGen if li is not None else '', claim=li.claim_first_sentence if li is not None and isinstance(li.claim_first_sentence, str) else '',
        mod=top_mod)

# ------------------------------------------------------------------------------------------------ module enrichment (step 6): categories (all rows) and pathways (called rows)
def module_block(df):
    names = sorted(df.module.unique())
    mi = {n: i for i, n in enumerate(names)}
    rows = []
    for r in df.itertuples():
        if r.key not in gi:
            continue
        f = int(r.call_NextGen) | (int(r.call_DepMap2D) << 1) | (int(r.call_Cas12a2D) << 2)
        rows.append([gi[r.key], mi[r.module], float(f'{r.FDR_NextGen:.2g}'), float(f'{r.FDR_DepMap2D:.2g}'), float(f'{r.FDR_Cas12a2D:.2g}'),
                     int(r.overlap_NextGen), int(r.overlap_DepMap2D), f])
    return dict(names=names, sizes=[int(df[df.module == n].n_members.iloc[0]) for n in names], rows=rows)


modules = dict(category=module_block(mcat), pathway=module_block(mpw))

# ------------------------------------------------------------------------------------------------ lineage-restricted pairs (NextGen)
Lrows = [[gi[r.gene_a], gi[r.gene_b], r.lineage, round(float(r.r), 3), float(f'{r.FDR:.2g}'), round(float(r.r_other_lineages), 3), r.restriction] for r in Lr.itertuples()]

# ------------------------------------------------------------------------------------------------ overview layout (RWE call graph, union)
ov = calls[['gene_a', 'gene_b', 'r_NextGen', 'r_DepMap2D']].copy()
ov['w'] = np.fmax(ov.r_NextGen.abs().fillna(0), ov.r_DepMap2D.abs().fillna(0))
nodes = sorted(set(ov.gene_a) | set(ov.gene_b) | set(k for k in rwe.index if k in gi))
ni = {k: i for i, k in enumerate(nodes)}
g = ig.Graph(n=len(nodes), edges=[(ni[x], ni[y]) for x, y in zip(ov.gene_a, ov.gene_b)])
g.es['weight'] = list(ov.w)
lay = np.array(g.layout_drl(weights='weight', seed=np.random.RandomState(7).rand(len(nodes), 2) * 10).coords)
lay = (lay - lay.mean(0)) / np.abs(lay - lay.mean(0)).max() * 2400
layout = {'idx': [gi[k] for k in nodes], 'x': [int(v) for v in lay[:, 0]], 'y': [int(v) for v in lay[:, 1]]}
print(f'overview graph: {len(nodes):,} nodes, {len(ov):,} edges')

# ------------------------------------------------------------------------------------------------ screens and gene-effect vectors
LCOL = {'CNS/Brain': '#D55E00', 'Esophagus/Stomach': '#0072B2', 'Pancreas': '#009E73', 'Bowel': '#E69F00', 'Breast': '#CC79A7', 'Prostate': '#56B4E9', 'Ovary/Fallopian Tube': '#7F3C8D',
        'Uterus': '#A5AA99', 'Biliary Tract': '#8C564B', 'Ampulla of Vater': '#B8860B'}
LINEAGES = list(LCOL)
screens, vectors = {}, {}
cas12a_2d = lambda: [x for x in sm[sm.IsTraditional2D & sm.Library.eq('Humagne-CD')].index if x in ge.index]  # noqa: E731  library check
for d, scr in list(data.DATASETS.items()) + [('Cas12a2D', cas12a_2d)]:
    s_ = list(scr())
    m = sm.loc[s_]
    screens[d] = dict(id=s_, model=m.ModelID.tolist(), lin=[ix(LINEAGES, v) for v in m.OncotreeLineage], type=m.ScreenType.tolist(), lib=m.Library.tolist(),
                      sub=m.OncotreeSubtype.fillna('').tolist())
    V = ge.loc[s_].reindex(columns=genes).values.T                                      # genes x screens
    Q = np.where(np.isnan(V), -128, np.clip(np.round(V * 40), -127, 127)).astype(np.int8)
    vectors[d] = b64gz(Q)
    print(f'{d}: vectors {Q.shape} -> {len(vectors[d]) / 1e6:.1f} MB base64')

bundle = dict(
    meta=dict(version='3.0.0', built=datetime.date.today().isoformat(), fdr=0.005, rmin=0.3, n=dict(NextGen=len(data.nextgen_screens()), DepMap2D=len(data.traditional_lineage_matched_screens()), Cas12a2D=len(cas12a_2d())),
              n_family=F, n_calls=dict(NextGen=int(CALL['NextGen'].sum() // 2), DepMap2D=int(CALL['DepMap2D'].sum() // 2))),
    dict=dict(cls=CLASSES, cat=CATS, fam=FAMS, role=ROLES, lineage=LINEAGES, lcol=[LCOL[L] for L in LINEAGES]),
    genes=G, edges=dict(a=b64(a.astype(np.uint16)), b=b64(b.astype(np.uint16)), rng=b64(q8(r_ng)), r2d=b64(q8(r_2d)), f=b64(flags.astype(np.uint8)), n=len(pairs)),
    cards=cards, modules=modules, lineage_pairs=Lrows, layout=layout, screens=screens, vectors=vectors)
out = os.path.join(paths.INTERMEDIATE, 'explorer_data.json')
with open(out, 'w') as fh:
    json.dump(bundle, fh, separators=(',', ':'), ensure_ascii=False, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
print(f'wrote {os.path.relpath(out, paths.ROOT)} ({os.path.getsize(out) / 1e6:.1f} MB)')
