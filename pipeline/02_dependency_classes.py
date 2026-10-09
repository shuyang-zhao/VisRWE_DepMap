"""Step 2 - gene dependency classes in NextGen (paper) and in the DepMap 2D comparator (same rules).

The co-dependency family of the paper's correlation analyses is restricted to 'dependencies' (Strongly Selective, High Variance,
Pan-dependency).  NextGen classes come from the paper; the same classification is re-implemented (src/rwe_nextgen/tda.py), checked
against the paper's NextGen classes, and applied to the DepMap 2D screens so that both datasets are filtered the same way.
Output: results/tables/dependency_classes.tsv
"""
import pandas as pd
from _common import banner, data, save
from rwe_nextgen import tda

banner('2. dependency classes (NextGen: paper; DepMap 2D: same rules)')
ge, dep = data.screen_effect(), data.screen_dependency()
ex = data.expression()
paper = data.tda_classes()
out = {}
for name, scr in (('NextGen', data.nextgen_screens()), ('DepMap2D', data.traditional_lineage_matched_screens())):
    scr = [s for s in scr if s in dep.index]
    res, info = tda.classify(ge.loc[scr], dep.loc[scr], data.to_screens(ex, scr))
    out[name] = res
    print(f'{name}: {len(scr)} screens | common-essential rank threshold {info["common_essential_rank_threshold"]:.3f} | expression cut-off {info["expression_cutoff"]:.2f} | classes {res.Class.value_counts().to_dict()}')
cmp_ = pd.DataFrame({'paper': paper.NextGenClass, 'reimplementation': out['NextGen'].Class.reindex(paper.index)})
agree = (cmp_.paper == cmp_.reimplementation).mean()
print(f'\ncheck: re-implemented NextGen classes agree with the paper for {agree:.2%} of {len(cmp_):,} genes')
print(pd.crosstab(cmp_.paper, cmp_.reimplementation).to_string())
T = pd.DataFrame({'NextGenClass': paper.NextGenClass, 'OrganoidsClass': paper.OrganoidsClass, 'NextGenClass_reimplemented': cmp_.reimplementation})
T = T.join(out['DepMap2D'].add_prefix('DepMap2D'))
T.index.name = 'key'
DEP = ['Strongly Selective', 'High Variance', 'Pan-dependency']
T['dependency_NextGen'] = T.NextGenClass.isin(DEP)
T['dependency_DepMap2D'] = T.DepMap2DClass.isin(DEP)
T['in_codependency_family'] = T.dependency_NextGen | T.dependency_DepMap2D
print(f'\ndependencies: NextGen {int(T.dependency_NextGen.sum()):,} | DepMap 2D {int(T.dependency_DepMap2D.sum()):,} | union (co-dependency family) {int(T.in_codependency_family.sum()):,}')
print(pd.crosstab(T.NextGenClass, T.DepMap2DClass).to_string())
save(T.reset_index(), 'tables', 'dependency_classes.tsv')
