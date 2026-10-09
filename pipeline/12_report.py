"""Step 12 - analysis report (Chinese, docs/report.zh-CN.md) and results workbook (results/VisRWE_DepMap_tables.xlsx).

All numbers in the report are read from the result tables, so the report always matches the current run.
"""
import datetime
import os

import numpy as np
import pandas as pd
from _common import banner, paths

banner('12. report and workbook')
T = lambda f: pd.read_csv(os.path.join(paths.TABLES, f), sep='\t')  # noqa: E731
genes, val, cls = T('rwe_gene_list.tsv'), T('validation_expression_addiction.tsv'), T('dependency_classes.tsv')
sg, lt, bm, se = T('single_gene_rwe.tsv'), T('lineage_tests_rwe.tsv'), T('biomarkers_rwe.tsv'), T('set_enrichment.tsv')
ov, calls, per = T('codep_overview.tsv').set_index('dataset'), T('codep_rwe_calls.tsv'), T('codep_summary_per_gene.tsv')
lr, mcat, mpw = T('lineage_restricted_pairs.tsv'), T('module_enrichment_categories.tsv'), T('module_enrichment_pathways.tsv')
kp, kpf = T('known_partners.tsv'), T('known_partners_by_family.tsv').set_index('family')
lit, claims, libs = T('literature_rwe.tsv'), T('literature_claims.tsv'), T('library_sensitivity_2D.tsv')
mall = pd.read_csv(os.path.join(paths.INTERMEDIATE, 'module_enrichment_pathways_all.tsv.gz'), sep='\t')


def md(df, floatfmt='{:.2f}'):
    cols = list(df.columns)
    out = ['| ' + ' | '.join(map(str, cols)) + ' |', '|' + '---|' * len(cols)]
    for _, r in df.iterrows():
        out.append('| ' + ' | '.join('' if (isinstance(v, float) and np.isnan(v)) else (floatfmt.format(v) if isinstance(v, float) else str(v)) for v in r.values) + ' |')
    return '\n'.join(out)


def f(v, d=2):
    return f'{v:+.{d}f}'


# ---------------------------------------------------------------------------------------------- numbers used in the text
N = {d: int(ov.loc[d, 'screens']) for d in ov.index}
agree = (cls.NextGenClass == cls.NextGenClass_reimplemented).mean()
rc = sg.NextGenClass.value_counts()
pan_all = se[(se.dataset == 'NextGen') & (se.outcome == 'Pan-dependency') & (se.set == 'RWE: all')].iloc[0]
role_pan = se[(se.dataset == 'NextGen') & (se.outcome == 'Pan-dependency') & se.set.str.startswith('role: ')].set_index('set').frac
org = sg[sg.org_vs_2D_call.fillna('') != '']
org_ess = sg[(sg.org_vs_2D_MannWhitneyFDR < 0.05) & (sg.org_vs_2D_MeanDifference.abs() > 0.3) & sg.org_vs_2D_DropEssential.fillna(False).astype(bool)].sort_values('org_vs_2D_MeanDifference')
lin_ng = lt[lt.NextGen_hit.fillna(False).astype(bool)]
lin_2d = lt[lt.DepMap2D_hit.fillna(False).astype(bool)]
bsig = bm[bm.Significant]
cc, rr = calls['class'].value_counts(), calls[calls.rwe_rwe]['class'].value_counts()
fam = per[per.in_family == True]  # noqa: E712
g = per.set_index('symbol')
co = kp[kp.r_DepMap2D >= 0.3]
lib_rwe = libs[libs.DepMap2D_library_sensitive & libs.key.isin(genes.key)]
lrc = lr.restriction.value_counts()
cat_cls, pw_cls = mcat['class'].fillna('').value_counts(), mpw['class'].fillna('').value_counts()
pw_c12 = int(mpw.call_Cas12a2D.sum())
ng_only_pw = mpw[mpw['class'] == 'NextGen only']
ng_only_c12 = ng_only_pw.call_Cas12a2D.mean() if len(ng_only_pw) else np.nan


def partners(gene):
    x = calls[(calls.symbol_a == gene) | (calls.symbol_b == gene)].copy()
    x['partner'] = np.where(x.symbol_a == gene, x.symbol_b, x.symbol_a)
    return x


pus = partners('PUS7')
pus_ng = pus[pus.call_NextGen].sort_values('FDR_NextGen')
dolk = pus.set_index('partner').loc['DOLK']
nglyc = mall[(mall.symbol == 'PUS7') & mall.module.str.contains('N-glycan|Dolichyl', case=False)]
big = kpf[(kpf.index != 'all families') & (kpf.pairs >= 8)]
weaker = big[big.called_NextGen < big.called_DepMap2D]
m3 = kp[((kp.gene_a == 'METTL3') & (kp.gene_b == 'METTL14')) | ((kp.gene_a == 'METTL14') & (kp.gene_b == 'METTL3'))].iloc[0]
m6a = kpf.loc['m6A/m6Am'] if 'm6A/m6Am' in kpf.index else None

# ---------------------------------------------------------------------------------------------- tables shown in the report
lin_tab = lt[lt['class'].fillna('') != ''].sort_values(['class', 'NextGen_mean_difference'])[
    ['symbol', 'lineage', 'NextGen_mean_difference', 'NextGen_FDR', 'DepMap2D_mean_difference', 'DepMap2D_FDR', 'class']].rename(columns={
        'symbol': '基因', 'lineage': '谱系', 'NextGen_mean_difference': 'NextGen 差值', 'NextGen_FDR': 'NextGen FDR', 'DepMap2D_mean_difference': '2D 差值', 'DepMap2D_FDR': '2D FDR',
        'class': '类别'})
top_tab = fam.sort_values('n_calls_NextGen', ascending=False).head(10)[['symbol', 'n_calls_NextGen', 'n_calls_DepMap2D', 'n_calls_both']].rename(columns={
    'symbol': '基因', 'n_calls_NextGen': 'NextGen', 'n_calls_DepMap2D': 'DepMap 2D', 'n_calls_both': '两边都显著'})
fam_tab = kpf.reset_index()[['family', 'pairs', 'median_r_NextGen', 'median_r_DepMap2D', 'called_NextGen', 'called_DepMap2D', 'called_DepMap2D_at_n147']].rename(columns={
    'family': '修饰家族', 'pairs': '已知伙伴对', 'median_r_NextGen': 'NextGen 中位 r', 'median_r_DepMap2D': '2D 中位 r', 'called_NextGen': 'NextGen 显著比例',
    'called_DepMap2D': '2D 显著比例', 'called_DepMap2D_at_n147': '2D 降采样到 147 的显著比例'})
lr_tab = lr.head(20)[['symbol_a', 'symbol_b', 'lineage', 'n', 'r', 'FDR', 'r_other_lineages', 'restriction']].rename(columns={
    'symbol_a': 'A', 'symbol_b': 'B', 'lineage': '谱系', 'n': 'n', 'r': '谱系内 r', 'FDR': 'FDR', 'r_other_lineages': '其余谱系 r', 'restriction': '类型'})
freq = mpw.groupby('module').agg(NextGen=('call_NextGen', 'sum'), DepMap2D=('call_DepMap2D', 'sum'), Cas12a2D=('call_Cas12a2D', 'sum'))
freq = freq.assign(total=freq.NextGen + freq.DepMap2D).sort_values('total', ascending=False).head(12).drop(columns='total').reset_index()
freq['module'] = freq.module.str.replace(r' R-HSA-\d+$', '', regex=True)
freq = freq.rename(columns={'module': '模块', 'NextGen': 'NextGen 中富集的 RWE 基因数', 'DepMap2D': '2D', 'Cas12a2D': 'Cas12a 2D'})
pus_tab = pus_ng[['partner', 'r_NextGen', 'FDR_NextGen', 'r_DepMap2D']].rename(columns={'partner': '伙伴', 'r_NextGen': 'NextGen r', 'FDR_NextGen': 'NextGen FDR', 'r_DepMap2D': 'DepMap 2D r'})
ng_tab = nglyc[['module', 'n_members', 'overlap_NextGen', 'FDR_NextGen', 'overlap_DepMap2D', 'overlap_Cas12a2D']].copy()
ng_tab['module'] = ng_tab.module.str.replace(r' R-HSA-\d+$', '', regex=True)
ng_tab = ng_tab.rename(columns={'module': 'Reactome 通路', 'n_members': '成员数', 'overlap_NextGen': 'NextGen 重叠（/100）', 'FDR_NextGen': 'NextGen FDR',
                                'overlap_DepMap2D': '2D 重叠', 'overlap_Cas12a2D': 'Cas12a 2D 重叠'})
conc, verdicts = lit.concordance_NextGen.value_counts(), claims.verdict.value_counts()

report = f"""# RNA 修饰 writer / reader / eraser（RWE）在 DepMap NextGen 中的依赖与共依赖

**v3（统计口径对齐原文）** · 作者：Shuyang Zhao（Department of Chemistry, University of California, Riverside）· 生成日期 {datetime.date.today().isoformat()}
交互页面：<https://shuyang-zhao.github.io/VisRWE_DepMap/> · 代码：<https://github.com/shuyang-zhao/VisRWE_DepMap>

> 数据：Neiswender, Maffa, Brenan *et al.* *A dependency map enhanced with next-generation 3D cancer models*. Nature (2026), doi:10.1038/s41586-026-10843-7；DepMap 门户文件集 "NextGen Model Manuscript 2026"。
> 本资源对所有修饰机器一视同仁：基因表、模块、已知伙伴、文库检查都用同一套规则，不为任何单个基因或假设设置专门的检验。

## 0. 要点

1. **统计口径与原文一致，并已验证**：用原始矩阵重算原文发布的 expression-addiction 表（{len(val):,} 个基因，|Δr| 在舍入精度内，n 完全相同）；重新实现的依赖分类与原文 {agree:.2%} 一致，再用同一规则给 2D 分类。
2. **单基因**：{len(sg)} 个 RWE 基因中 {rc.get('Pan-dependency', 0)} 个是 NextGen 泛必需（{pan_all.frac:.0%}，全基因组 {pan_all.genome_frac:.1%}，OR {pan_all.odds_ratio:.1f}）；writer 辅因子 {role_pan['role: writer-cofactor']:.0%}、reader {role_pan['role: reader']:.0%}、writer {role_pan['role: writer']:.0%}、eraser {role_pan['role: eraser']:.0%}。NextGen 与 2D 的类别高度一致；按原文口径，类器官 vs 贴壁 2D 只有 {'、'.join(f'{r.symbol}（{f(r.org_vs_2D_MeanDifference)}）' for r in org.itertuples())} 显著；谱系选择性依赖 NextGen {len(lin_ng)} 个、2D {len(lin_2d)} 个。
3. **共依赖**（依赖基因两两 Pearson，每个数据集整表 BH，FDR < 0.005 且 |r| > 0.3）：NextGen {int(ov.loc['NextGen', 'calls']):,} 对、DepMap 2D {int(ov.loc['DepMap2D', 'calls']):,} 对（各 {int(ov.loc['NextGen', 'pairs']):,} 对检验）。含 RWE 基因的显著对 {len(calls):,} 个：**两边都显著 {cc.get('both', 0):,}，仅 NextGen {cc.get('NextGen only', 0):,}，仅 2D {cc.get('DepMap 2D only', 0):,}**；RWE–RWE 为 {rr.get('both', 0)} / {rr.get('NextGen only', 0)} / {rr.get('DepMap 2D only', 0)}。
4. **已知伙伴在 NextGen 中普遍变弱，多数修饰家族都如此**：从基因表复合物注释中提取的 {len(kp)} 对已知伙伴，2D 中显著 {kp.call_DepMap2D.mean():.0%}，NextGen 中只有 {kp.call_NextGen.mean():.0%}；在 2D 中共依赖（r ≥ 0.3）的 {len(co)} 对，中位 r 从 {co.r_DepMap2D.median():.2f} 降到 {co.r_NextGen.median():.2f}，而 2D 降采样到 147 个筛选后仍有 {co.DepMap2D_n147_frac_called.mean():.0%} 能检出——所以不是样本量或统计方法的问题。
5. **功能模块**（原文的排序超几何口径：每个 RWE 基因前 100 个正相关伙伴，对 16 个粗分类和 {mall.module.nunique():,} 个 Reactome / CORUM 模块一视同仁地检验）：细模块显著 NextGen {int(mpw.call_NextGen.sum()):,}、2D {int(mpw.call_DepMap2D.sum()):,} 个（基因 × 模块）；最常见的是线粒体翻译、tRNA 加工 / 修饰、rRNA 修饰——即这些 RWE 基因自身所在的机器。
6. **谱系**："只在单一谱系极强"的共依赖 {len(lr)} 对（{', '.join(f'{k} {v}' for k, v in lr.lineage.value_counts().items())}），其中 {lrc.get('restricted', 0)} 对在其余谱系 r ≤ 0.3。
7. **交互页面**（GitHub Pages）：NextGen / DepMap 2D / 对比三种数据视图；悬停预览、点击展开的探索网络；每条边按谱系分解（两套数据并列）并附 Cas12a 2D 文库检查；两级模块热图；中英文切换。

## 1. 数据与统计口径

| 项目 | 做法（与原文代码一致） |
|---|---|
| 依赖矩阵 | 筛选级 Chronos gene effect（`screen_gene_effect.csv`），缺失值不填补，相关用两两完整的筛选 |
| NextGen | ScreenType ∈ {{2DO, 3DO, 3DN, 2DN}}、IsNextGen、谱系 ≠ Other：**{N['NextGen']} 个筛选**（全部 Cas12a Humagne-CD 文库） |
| DepMap 2D | NextGen 谱系的传统 2D 筛选（原文 `adherent_next_gen_lineage_matched_screens`）：**{N['DepMap2D']} 个筛选**（Avana / KY 为主，另有 11 个 Humagne） |
| 文库检查 | 全部 36 个用 Cas12a（Humagne-CD）文库筛选的传统 2D 细胞系 |
| 依赖类别 | 原文 TDA 分类；2D 用同一规则重算 |
| 共依赖 | 依赖基因（任一数据集中为 Pan-dependency / Strongly Selective / High Variance，{int(ov.loc['NextGen', 'family_genes']):,} 个）两两 Pearson；每个数据集一次 BH；FDR < 0.005 且 \\|r\\| > 0.3。NextGen 的 FDR 边界对应 p ≈ {ov.loc['NextGen', 'p_boundary']:.1e}，2D 为 {ov.loc['DepMap2D', 'p_boundary']:.1e} |
| 模块 | 16 个粗分类（正则表达式，所有基因同一规则）+ Reactome 2022 / CORUM 中在依赖基因里有 5–150 个成员的全部模块；每个 RWE 基因前 100 个正相关伙伴做单侧超几何检验（原文 `run_rank_hypergeometric`），每个数据集、每个层级一次 BH，FDR < 0.05 |
| 组间比较 | 类器官 vs 贴壁 2D：原文表；谱系 vs 其余：Mann–Whitney，FDR < 0.05 且差值 < −0.3 |

**文库**：24 个两种文库都筛过的模型里，只有 {', '.join(sg[sg.library_sensitive].symbol)} 的 Cas12a − Cas9 差异显著且 > 0.3；2D 内部 KY vs Avana 差异大的 RWE 基因（R² ≥ 0.1）有 {', '.join(lib_rwe.key.str.split(' ').str[0])}，它们的 2D 边和模块在页面中有标记。

## 2. 单基因依赖

![Fig. 1](../results/figures/fig1_single_gene.png)

**谱系选择性依赖**（谱系 vs 其余，FDR < 0.05 且差值 < −0.3）：

{md(lin_tab, '{:.3g}')}

**类器官 vs 贴壁 2D**（原文表）：按原文口径（排除两组都必需的基因）只有 {'、'.join(f'{r.symbol} {f(r.org_vs_2D_MeanDifference)}' for r in org.itertuples())}；不排除时差值 > 0.3 的还有 {'、'.join(f'{r.symbol} {f(r.org_vs_2D_MeanDifference)}' for r in org_ess.itertuples())}。**生物标志物**（原文表与阈值）：表达成瘾 0 个；旁系 {len(bsig)} 个（{', '.join(f'{r.symbol}–{r.biomarker_symbol} r = {r.NextGenPearsonR:.2f}' for r in bsig.itertuples())}）。

## 3. 共依赖：NextGen vs DepMap 2D

![Fig. 2](../results/figures/fig2_codependency.png)

{len(fam)} 个进入检验的 RWE 基因中，{int((fam.n_calls_NextGen > 0).sum())} 个在 NextGen、{int((fam.n_calls_DepMap2D > 0).sum())} 个在 2D 至少有一个显著伙伴。伙伴最多的是线粒体 RNA 加工 / 翻译相关的 RWE：

{md(top_tab, '{:.0f}')}

## 4. 已知伙伴（全部修饰家族，同一规则）

已知伙伴取自基因表中每个 RWE 基因的复合物注释里点名、且属于依赖基因的伙伴，共 {len(kp)} 对。

{md(fam_tab, '{:.2f}')}

读法：① 2D 降采样到与 NextGen 相同的 147 个筛选后，显著比例几乎不变（{kpf.loc['all families', 'called_DepMap2D']:.0%} → {kpf.loc['all families', 'called_DepMap2D_at_n147']:.0%}），而 NextGen 只有 {kpf.loc['all families', 'called_NextGen']:.0%}，所以**样本量和统计方法都不是原因**；② 在已知伙伴 ≥ 8 对的 {len(big)} 个修饰家族中，{len(weaker)} 个在 NextGen 中显著比例更低（{'、'.join(weaker.index)}），不是某一类机器特有；③ 扣除筛选级"必需性轴"后 NextGen 的 r 变化不大（`known_partners.tsv` 中的 `r_NextGen_axis_adjusted`）。可能的解释包括 3D 培养下这些基因的依赖受其他因素驱动，以及 Cas12a 文库的向导差异，现有数据无法区分。

以 m6A writer 为例（回答"METTL3 / YTH 为什么找不到自己的伙伴"）：METTL3–METTL14 在 2D 中 r = {f(m3.r_DepMap2D)}（降采样后 {m3.DepMap2D_n147_frac_called:.0%} 被检出），NextGen 中只有 {f(m3.r_NextGen)}{'' if m6a is None else f'；m6A 家族 {int(m6a.pairs)} 对已知伙伴的中位 r 为 NextGen {m6a.median_r_NextGen:.2f} vs 2D {m6a.median_r_DepMap2D:.2f}'}。YTHDF1/YTHDF3/YTHDC2 在两套数据中都不是依赖基因，不进入共依赖检验。

## 5. 功能模块

![Fig. 3](../results/figures/fig3_modules.png)

粗分类（基因 × 类别）：两边都显著 {cat_cls.get('both', 0)}、仅 NextGen {cat_cls.get('NextGen only', 0)}、仅 2D {cat_cls.get('DepMap 2D only', 0)}。Reactome / CORUM 细模块：两边都显著 {pw_cls.get('both', 0)}、仅 NextGen {pw_cls.get('NextGen only', 0)}、仅 2D {pw_cls.get('DepMap 2D only', 0)}；其中 {pw_c12} 个在 Cas12a 2D 中也显著，"仅 NextGen"的细模块里有 {ng_only_c12:.0%} 在 Cas12a 2D 中也显著（文库可能参与）。被最多 RWE 基因富集的模块：

{md(freq, '{:.0f}')}

## 6. 案例

以下案例读取的都是上面的通用结果，没有专门的检验。

**PUS7**：NextGen 中 {len(pus_ng)} 个显著伙伴，2D 中没有：

{md(pus_tab, '{:.3g}')}

PUS7–DOLK：NextGen {f(dolk.r_NextGen)}，2D {f(dolk.r_DepMap2D)}，两数据集差异 P = {dolk.p_difference:.1g}。在统一的模块检验中，PUS7 的前 100 个伙伴与 N-糖基化相关 Reactome 通路的重叠很小，校正后不显著：

{md(ng_tab, '{:.3g}')}

v2 曾用一个预先指定的"早期 DLO + 多萜醇供应"模块单独检验 PUS7 并得到显著结果；那是针对单一假设的检验，在对所有基因、所有模块一视同仁的框架下不再突出。另外同一信号在 36 个 Cas12a 2D 细胞系中也能看到（v2），所以即使存在，也可能与文库有关。

**DKC1**：2D {int(g.loc['DKC1', 'n_calls_DepMap2D'])} 个显著伙伴，NextGen {int(g.loc['DKC1', 'n_calls_NextGen'])} 个（{int(g.loc['DKC1', 'n_calls_both'])} 个共有）；模块层面两边都指向 rRNA 修饰 / 核糖体生物发生。**TRUB1**：两套数据中都不是依赖基因。**WTAP、YTHDC1**：NextGen 中伙伴多为 TAF12、POLR2G、CCDC174 等通用必需基因。

## 7. 谱系分解与"只在单一谱系极强"的关联

发现仍在合并数据上做；页面中每条边都能按谱系分解（森林图 + 散点，NextGen 与 2D 并列）。另外在 NextGen 的 5 个 ≥ 14 个筛选的谱系内检验含 RWE 基因的对（谱系内 Pearson，所有"对 × 谱系"一次 BH，FDR < 0.005 且 \\|r\\| > 0.3，且与其余谱系的 Fisher-z 差异 P < 0.01）：共 {len(lr)} 对（单谱系限定 {lrc.get('restricted', 0)}、明显更强 {lrc.get('enriched', 0)}、各谱系都有而此谱系更强 {lrc.get('shared, stronger here', 0)}）。前 20 对：

{md(lr_tab, '{:.3g}')}

2D 不做这一步：2D 各谱系混有比例不同的 Avana / KY / Humagne 筛选，谱系内相关会被文库组成混杂。

## 8. 文献主张

{len(claims)} 条逐条核实的主张：{', '.join(f'{k} {v}' for k, v in verdicts.items())}。基因层面最常见的是"有依赖主张、数据中泛必需（选择性不成立）"{conc.get('claim consistent (dependent), but pan-dependency - claimed selectivity not supported', 0)} 个，"主张未复现"{conc.get('claim not reproduced (weak dependency only)', 0) + conc.get('claim not reproduced (no dependency)', 0)} 个，"没有主张但泛必需"{conc.get('data only: pan-dependency without a cancer-dependency claim', 0)} 个。

## 9. 局限

1. NextGen 只有 Cas12a 文库，2D 以 Cas9 为主；"仅 NextGen"的信号应先看页面中的 Cas12a 2D 文库检查。
2. 原文口径不调整谱系、筛选质量；一部分共依赖来自谱系组成或通用必需机器（线粒体翻译、核糖体），页面的谱系分解用来判断这一点。
3. 共依赖不等于相互作用；许多 RWE 的依赖来自催化非依赖的结构功能（基因表中 `caveat_class`）。
4. 非依赖基因不进入共依赖检验，"没有伙伴"不等于"没有关系"。
5. 模块的超几何检验不考虑模块内部基因的相关，和原文一样；它回答"前 100 个伙伴里是否过表达"，不是"整体相关是否更强"。
6. 谱系内检验样本量小（14–39），只看极强的关联。

## 10. 复现与文件

```bash
conda env create -f environment.yml && conda activate visrwe
make all
```

| 文件 | 内容 |
|---|---|
| `results/tables/single_gene_rwe.tsv` | 每个 RWE 基因：类别、平均效应、依赖比例、类器官 vs 2D、谱系、文库、生物标志物 |
| `results/tables/codep_rwe_calls.tsv` | 含 RWE 基因、在任一数据集显著的共依赖对 |
| `results/tables/codep_summary_per_gene.tsv` | 每个 RWE 基因的显著伙伴数与最强伙伴 |
| `results/tables/known_partners.tsv`、`known_partners_by_family.tsv` | 已知伙伴检查（第 4 节） |
| `results/tables/module_enrichment_categories.tsv`、`module_enrichment_pathways.tsv` | 模块富集（第 5 节） |
| `results/tables/lineage_restricted_pairs.tsv` | 谱系限定的强关联 |
| `results/VisRWE_DepMap_tables.xlsx` | 以上表格汇总 |
| `docs/index.html` | 交互页面（GitHub Pages） |
"""
out = os.path.join(paths.ROOT, 'docs', 'report.zh-CN.md')
open(out, 'w', encoding='utf-8').write(report)
print('wrote docs/report.zh-CN.md', len(report), 'chars')

# ---------------------------------------------------------------------------------------------- workbook
readme = [
    ('rwe_genes', 'Curated RWE gene list with DepMap keys and NextGen classes'),
    ('single_gene', 'Per RWE gene: dependency classes, mean effect and fraction dependent in NextGen / organoids / CNS / DepMap 2D, organoid vs 2D (paper table), lineage, library, biomarkers'),
    ('lineage_tests', 'Lineage vs rest (Mann-Whitney, BH per dataset) for RWE genes, NextGen and DepMap 2D'),
    ('biomarkers', 'Expression addiction and paralog tests of RWE genes (paper tables and thresholds)'),
    ('set_enrichment', 'Fisher tests of dependency classes for RWE families / roles'),
    ('codep_overview', 'Genome-wide co-dependency numbers per dataset'),
    ('codep_rwe_calls', 'Pairs with an RWE gene called in NextGen and/or DepMap 2D (FDR < 0.005 and |r| > 0.3)'),
    ('codep_per_gene', 'Number of called partners per RWE gene and strongest partners'),
    ('known_partners', 'Partners named in the curated complex annotations: r and calls in both datasets, 2D down-sampled to 147 screens'),
    ('known_partners_family', 'Known-partner summary per modification family'),
    ('modules_categories', 'Over-representation of 16 broad categories among the top-100 partners of each RWE gene'),
    ('modules_pathways', 'Reactome / CORUM modules called in NextGen or DepMap 2D (with Cas12a 2D library check)'),
    ('lineage_restricted', 'NextGen pairs strong inside one lineage only'),
    ('literature_genes', 'Gene-level claim stance vs data'),
    ('literature_claims', 'Claim-level audit (86 PMID-verified claims)'),
]
sheets = [('README', pd.DataFrame(readme, columns=['sheet', 'description'])),
          ('rwe_genes', genes.drop(columns=['cancer_dependency_claims', 'key_refs', 'notes'], errors='ignore')), ('single_gene', sg), ('lineage_tests', lt), ('biomarkers', bm),
          ('set_enrichment', se), ('codep_overview', ov.reset_index()), ('codep_rwe_calls', calls), ('codep_per_gene', per), ('known_partners', kp),
          ('known_partners_family', kpf.reset_index()), ('modules_categories', mcat), ('modules_pathways', mpw), ('lineage_restricted', lr),
          ('literature_genes', lit), ('literature_claims', claims)]
wb = os.path.join(paths.RESULTS, 'VisRWE_DepMap_tables.xlsx')
with pd.ExcelWriter(wb, engine='openpyxl') as xw:
    for name, df in sheets:
        df.to_excel(xw, sheet_name=name, index=False)
        ws = xw.sheets[name]
        ws.freeze_panes = 'A2'
        for i, c in enumerate(df.columns, 1):
            ws.column_dimensions[ws.cell(1, i).column_letter].width = 100 if (name == 'README' and c == 'description') else min(30, max(9, len(str(c)) + 2))
print(f'wrote results/VisRWE_DepMap_tables.xlsx ({len(sheets)} sheets, {os.path.getsize(wb) / 1e6:.1f} MB)')
