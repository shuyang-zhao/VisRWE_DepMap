"""Step 7 - published cancer-dependency claims vs the data.

Gene level : the claim text of the curated RWE list is coded into a stance (rule-based, with manual corrections after reading the
             sentences) and compared with the dependency class in NextGen and DepMap 2D.
Claim level: 86 PMID-verified claims (data/curated/literature_claim_audit.tsv) with the verdict of the claim-by-claim audit; the
             current NextGen / DepMap-2D dependency profile of each claimed gene is attached.
Outputs: results/tables/literature_rwe.tsv, literature_claims.tsv
"""
import os
import re

import numpy as np
import pandas as pd
from _common import banner, data, paths, save

banner('7. literature claims vs data')
NONE = re.compile(r'^\s*(nan|none\b|no\b)', re.I)
SUPPRESSOR = re.compile(r'not a cancer dependency|tumou?r[- ]suppress|growth suppressor|opposite of a dependency|rather than (as )?a dependency|not as a dependency|SUPPRESSOR|NOT essential', re.I)
DEPENDENCY = re.compile(r'required for|requirement for|dependency|depend(s|ent) on|knockdown (reduces|suppresses|inhibits|impairs|lowers|causes)|knockout (reduces|impairs|causes)|'
                        r'depletion (kills|impairs|reduces|blocks|prolongs)|loss (arrests|impairs)|critical for|essential', re.I)
EXPRESSION = re.compile(r'up-?regulated|over-?expressed|over-?expression|oncogen|poor (survival|outcome|prognosis)|prognos', re.I)
NEGATED = re.compile(r'not a cancer dependency|rather than (as )?a dependency|not as a dependency|opposite of a dependency', re.I)
# manual corrections after reading the curated sentences (negations, expression-only and non-viability claims the rules mislabel)
OVERRIDE = {
    **{g: 'none' for g in ['ADAT2', 'RPUSD4', 'RPUSD3', 'TRUB2', 'TRMT10C', 'SARS1', 'TRMT112', 'BUD23', 'PRRC2A', 'PUS7L', 'RBM15B', 'YTHDC2', 'APOBEC3G', 'THADA', 'METTL4', 'FUS']},
    **{g: 'expression_or_prognostic_only' for g in ['FBL', 'EIF3A', 'NOP2', 'TRMT2A', 'A1CF', 'RPUSD2', 'OSGEPL1', 'PUS1']},
    **{g: 'non_viability_role' for g in ['FTSJ3', 'CENPA', 'PCF11', 'MYBBP1A', 'PABPC1', 'SRSF2', 'CBLL1', 'G3BP1', 'G3BP2', 'SIRT7']},
    **{g: 'suppressor_or_not_dependency' for g in ['APOBEC1', 'ELP1', 'NSUN6', 'RBM47']},
}


def stance(text):
    head = ('' if not isinstance(text, str) else text).split('||')[0].strip()[:600]
    if NONE.match(head):
        return 'none'
    supp, dep, expr = bool(SUPPRESSOR.search(head)), bool(DEPENDENCY.search(head)), bool(EXPRESSION.search(head))
    if supp and dep and not NEGATED.search(head):
        return 'mixed'
    if supp:
        return 'suppressor_or_not_dependency'
    if dep:
        return 'dependency_claimed'
    return 'expression_or_prognostic_only' if expr else 'none'


def level(c):
    if c == 'Pan-dependency':
        return 'pan'
    if c in ('Strongly Selective', 'High Variance'):
        return 'selective_or_variable'
    if c == 'Weakly Selective':
        return 'weak'
    return 'none' if isinstance(c, str) and c != 'Unscreened' else 'no_data'


CONCORDANCE = {
    'none': {'pan': 'data only: pan-dependency without a cancer-dependency claim', 'selective_or_variable': 'data only: selective / variable dependency without a claim',
             'weak': 'consistent: no claim, weak or no dependency', 'none': 'consistent: no claim, no dependency'},
    'dependency_claimed': {'pan': 'claim consistent (dependent), but pan-dependency - claimed selectivity not supported', 'selective_or_variable': 'claim consistent (selective / variable dependency)',
                           'weak': 'claim not reproduced (weak dependency only)', 'none': 'claim not reproduced (no dependency)'},
    'suppressor_or_not_dependency': {'pan': 'discordant: suppressor / non-dependency claim but pan-dependency', 'selective_or_variable': 'discordant: suppressor claim but selective dependency',
                                     'weak': 'consistent: not a dependency', 'none': 'consistent: not a dependency'},
    'mixed': {k: f'mixed claims; {v} in data' for k, v in (('pan', 'pan-dependency'), ('selective_or_variable', 'selective dependency'), ('weak', 'weak / no dependency'), ('none', 'no dependency'))},
    'non_viability_role': {k: f'non-viability claim only; {v}' for k, v in (('pan', 'pan-dependency'), ('selective_or_variable', 'selective dependency'), ('weak', 'no real dependency'), ('none', 'no dependency'))},
    'expression_or_prognostic_only': {k: f'expression / prognostic claim only; {v}' for k, v in (('pan', 'pan-dependency'), ('selective_or_variable', 'selective dependency'),
                                                                                                    ('weak', 'no real dependency'), ('none', 'no dependency'))},
}

sg = pd.read_csv(os.path.join(paths.TABLES, 'single_gene_rwe.tsv'), sep='\t').set_index('key')
u = data.rwe_genes().dropna(subset=['key']).drop_duplicates('key').copy()
u['claim_stance'] = [OVERRIDE.get(s, stance(t)) for s, t in zip(u.symbol, u.cancer_dependency_claims)]
u['stance_source'] = np.where(u.symbol.isin(OVERRIDE), 'manual', 'rule')
u['claim_first_sentence'] = [(t.split('||')[0].strip()[:260] if isinstance(t, str) else '') for t in u.cancer_dependency_claims]
for d in ('NextGen', 'DepMap2D'):
    u[f'{d}Class'] = u.key.map(sg[f'{d}Class'])
    u[f'{d}_mean'] = u.key.map(sg[f'{d}_mean'])
    u[f'{d}_frac_dep'] = u.key.map(sg[f'{d}_frac_dep'])
    u[f'data_level_{d}'] = u[f'{d}Class'].map(level)
    u[f'concordance_{d}'] = [('no data' if lv == 'no_data' else CONCORDANCE[s][lv]) for s, lv in zip(u.claim_stance, u[f'data_level_{d}'])]
u['NextGen_top_lineage'] = u.key.map(sg.NextGen_top_lineage)
u['NextGen_top_lineage_diff'] = u.key.map(sg.NextGen_top_lineage_diff)
u['NextGen_top_lineage_FDR'] = u.key.map(sg.NextGen_top_lineage_FDR)
aud = pd.read_csv(os.path.join(paths.CURATED, 'literature_claim_audit.tsv'), sep='\t')
counts = aud.groupby('gene').verdict.value_counts().unstack(fill_value=0)
counts.columns = ['audit: ' + c for c in counts.columns]
u = u.merge(counts, left_on='symbol', right_index=True, how='left')
cols = ['symbol', 'key', 'family', 'role', 'evidence_tier', 'caveat_class', 'claim_stance', 'stance_source', 'NextGenClass', 'DepMap2DClass', 'concordance_NextGen', 'concordance_DepMap2D',
        'NextGen_mean', 'NextGen_frac_dep', 'DepMap2D_mean', 'DepMap2D_frac_dep', 'NextGen_top_lineage', 'NextGen_top_lineage_diff', 'NextGen_top_lineage_FDR', 'claim_first_sentence'] + list(counts.columns)
save(u[cols], 'tables', 'literature_rwe.tsv')
aud['NextGenClass'] = aud.gene_key.map(sg.NextGenClass)
aud['DepMap2DClass'] = aud.gene_key.map(sg.DepMap2DClass)
save(aud, 'tables', 'literature_claims.tsv')
print(pd.crosstab(u.claim_stance, u.data_level_NextGen).to_string())
print('\nclaim-level audit verdicts:', aud.verdict.value_counts().to_dict())
print('\n', u.concordance_NextGen.value_counts().to_string())
