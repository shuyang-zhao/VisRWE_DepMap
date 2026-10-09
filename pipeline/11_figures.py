"""Step 11 - figures (PDF + PNG, 183 mm wide).

Fig. 1  single-gene dependencies of RWE genes in NextGen vs DepMap 2D
Fig. 2  co-dependencies of RWE genes: NextGen vs DepMap 2D, and known partners by modification family
Fig. 3  functional modules among the top partners of RWE genes (Reactome / CORUM), NextGen vs DepMap 2D
Outputs: results/figures/fig{1,2,3}_*.pdf/png
"""
import os

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
from _common import banner, paths  # noqa: E402

banner('11. figures')
MM = 1 / 25.4
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7, 'axes.titlesize': 7.5, 'axes.labelsize': 7, 'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5,
                     'legend.fontsize': 6.5, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6, 'pdf.fonttype': 42, 'savefig.dpi': 300})
NG, D2, BOTH, GREY = '#0b8a6f', '#c25b12', '#6d4bc4', '#9aa3ad'
CLASS_COL = {'Pan-dependency': '#C71000', 'Strongly Selective': '#145A32', 'High Variance': '#E08A00', 'Weakly Selective': '#9a9a9a', 'Non-dependency': '#2a3f66', 'Unscreened': '#dddddd'}
T = lambda f: pd.read_csv(os.path.join(paths.TABLES, f), sep='\t')  # noqa: E731
os.makedirs(paths.FIGURES, exist_ok=True)


def panel(ax, letter):
    ax.text(-0.18, 1.06, letter, transform=ax.transAxes, fontsize=9, fontweight='bold', va='bottom')


def save(fig, name):
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(paths.FIGURES, f'{name}.{ext}'), bbox_inches='tight')
    plt.close(fig)
    print('  wrote results/figures/' + name + '.pdf/.png')


# ---------------------------------------------------------------------------------------------------------- Fig. 1
sg = T('single_gene_rwe.tsv')
lt = T('lineage_tests_rwe.tsv')
fig = plt.figure(figsize=(183 * MM, 72 * MM), layout='constrained')
gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1, 1.15])
ax = fig.add_subplot(gs[0])
roles = ['writer', 'writer-cofactor', 'reader', 'reader-cofactor', 'eraser', 'other']

classes = list(CLASS_COL)
x = 0
ticks, labels = [], []
for role in roles:
    sub = sg[sg.role == role]
    for k, (col, lab_) in enumerate((('NextGenClass', 'NG'), ('DepMap2DClass', '2D'))):
        frac = sub[col].value_counts(normalize=True).reindex(classes).fillna(0)
        bottom = 0
        for c in classes:
            ax.bar(x + k * 0.42, frac[c], 0.4, bottom=bottom, color=CLASS_COL[c], linewidth=0)
            bottom += frac[c]
    ticks.append(x + 0.21)
    labels.append(f'{role} ({len(sub)})')
    x += 1.15
ax.set_xticks(ticks)
ax.set_xticklabels(labels, fontsize=5.8, rotation=40, ha='right', rotation_mode='anchor')
ax.set_ylabel('fraction of RWE genes')
ax.set_title('Dependency class\n(bars: NextGen | DepMap 2D)', loc='left')
for c in classes:
    ax.bar(0, 0, color=CLASS_COL[c], label=c)
ax.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.5), ncol=2, fontsize=5.8, handlelength=1, columnspacing=0.8)
panel(ax, 'a')
ax = fig.add_subplot(gs[1])
ax.scatter(sg.DepMap2D_mean, sg.NextGen_mean, s=7, c=[CLASS_COL.get(c, GREY) for c in sg.NextGenClass], linewidths=0, alpha=0.85)
lim = [min(sg.DepMap2D_mean.min(), sg.NextGen_mean.min()) - 0.1, max(sg.DepMap2D_mean.max(), sg.NextGen_mean.max()) + 0.1]
ax.plot(lim, lim, color=GREY, lw=0.6, ls='--')
lab_genes = sg[(sg.org_vs_2D_call.fillna('') != '') | (sg.NextGen_minus_2D_mean.abs() > 0.35)]
for r in lab_genes.itertuples():
    ax.annotate(r.symbol, (r.DepMap2D_mean, r.NextGen_mean), fontsize=5.5, xytext=(2, 2), textcoords='offset points')
ax.set_xlabel('mean gene effect, DepMap 2D')
ax.set_ylabel('mean gene effect, NextGen')
ax.set_title('Mean gene effect per RWE gene', loc='left')
panel(ax, 'b')
ax = fig.add_subplot(gs[2])
hits = lt[lt['class'].fillna('') != ''].copy()
hits['lab'] = hits.symbol + ' · ' + hits.lineage.str.replace('Esophagus/Stomach', 'Eso/Stomach')
hits = hits.sort_values('NextGen_mean_difference')
y = np.arange(len(hits))
ax.scatter(hits.NextGen_mean_difference, y, s=12, color=NG, label='NextGen', zorder=3)
ax.scatter(hits.DepMap2D_mean_difference, y, s=12, facecolors='white', edgecolors=D2, label='DepMap 2D', zorder=3)
for yy, a_, b_ in zip(y, hits.NextGen_mean_difference, hits.DepMap2D_mean_difference):
    ax.plot([a_, b_], [yy, yy], color=GREY, lw=0.6, zorder=1)
ax.axvline(0, color=GREY, lw=0.5)
ax.axvline(-0.3, color=GREY, lw=0.5, ls=':')
ax.set_yticks(y)
ax.set_yticklabels(hits.lab, fontsize=5.8)
ax.set_xlabel('lineage mean − rest (gene effect)')
ax.set_title('Lineage-selective dependencies\n(FDR < 0.05, diff. < −0.3)', loc='left')
ax.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2)
panel(ax, 'c')
save(fig, 'fig1_single_gene')

# ---------------------------------------------------------------------------------------------------------- Fig. 2
calls = T('codep_rwe_calls.tsv')
fig = plt.figure(figsize=(183 * MM, 70 * MM), layout='constrained')
gs = fig.add_gridspec(1, 3, width_ratios=[0.8, 1.1, 1.2])
ax = fig.add_subplot(gs[0])
order = ['both', 'NextGen only', 'DepMap 2D only']
cols = {'both': BOTH, 'NextGen only': NG, 'DepMap 2D only': D2}
for k, (name, sub) in enumerate((('RWE–other', calls[~calls.rwe_rwe]), ('RWE–RWE', calls[calls.rwe_rwe]))):
    cnt = sub['class'].value_counts().reindex(order).fillna(0)
    left = 0
    for c in order:
        ax.barh(k, cnt[c], left=left, color=cols[c], height=0.6)
        left += cnt[c]
    ax.text(left, k, f' {int(left):,}', va='center', fontsize=6)
ax.set_yticks([0, 1])
ax.set_yticklabels(['RWE–other', 'RWE–RWE'])
ax.set_xscale('symlog', linthresh=100)
ax.set_xlabel('called pairs')
ax.set_title('Called RWE pairs', loc='left')
for c in order:
    ax.barh(0, 0, color=cols[c], label=c)
ax.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.45, -0.22), ncol=1)
panel(ax, 'a')
ax = fig.add_subplot(gs[1])
for c in ('DepMap 2D only', 'NextGen only', 'both'):
    sub = calls[calls['class'] == c]
    ax.scatter(sub.r_DepMap2D, sub.r_NextGen, s=1.2, color=cols[c], alpha=0.45, linewidths=0, rasterized=True)
ax.plot([-0.8, 1], [-0.8, 1], color=GREY, lw=0.6, ls='--')
for v in (-0.3, 0.3):
    ax.axhline(v, color=GREY, lw=0.4, ls=':')
    ax.axvline(v, color=GREY, lw=0.4, ls=':')
ax.set_xlabel('r, DepMap 2D')
ax.set_ylabel('r, NextGen')
ax.set_title('r of pairs called in either dataset\n(n = %s; dotted: |r| = 0.3)' % f'{len(calls):,}', loc='left')
panel(ax, 'b')
ax = fig.add_subplot(gs[2])
kp = T('known_partners_by_family.tsv')
kp = kp[(kp.family != 'all families') & (kp.pairs >= 3)].sort_values('median_r_DepMap2D')
y = np.arange(len(kp))
for yy, r in zip(y, kp.itertuples()):
    ax.plot([r.median_r_NextGen, r.median_r_DepMap2D], [yy, yy], color=GREY, lw=0.7)
ax.scatter(kp.median_r_NextGen, y, s=10 + 2 * kp.pairs, color=NG, zorder=3, label='NextGen')
ax.scatter(kp.median_r_DepMap2D, y, s=10 + 2 * kp.pairs, facecolors='white', edgecolors=D2, zorder=3, label='DepMap 2D')
ax.axvline(0, color=GREY, lw=0.5)
ax.set_yticks(y)
ax.set_yticklabels([f'{f} ({n})' for f, n in zip(kp.family, kp.pairs)], fontsize=5.8)
ax.set_xlabel('median r of known partner pairs')
ax.set_title('Known partners by modification family\n(curated complex annotations)', loc='left')
ax.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.4, -0.2), ncol=2)
panel(ax, 'c')
save(fig, 'fig2_codependency')

# ---------------------------------------------------------------------------------------------------------- Fig. 3
P = pd.read_csv(os.path.join(paths.INTERMEDIATE, 'module_enrichment_pathways_all.tsv.gz'), sep='\t')
fams = T('rwe_gene_list.tsv').set_index('key').family
fig = plt.figure(figsize=(183 * MM, 72 * MM), layout='constrained')
gs = fig.add_gridspec(1, 2, width_ratios=[0.9, 1.3])
ax = fig.add_subplot(gs[0])
per = P.groupby('key').agg(ng=('call_NextGen', 'sum'), d2=('call_DepMap2D', 'sum'), c12=('call_Cas12a2D', 'sum'))
per['family'] = fams.reindex(per.index).values
fam_list = sorted(per.family.dropna().unique())
pal = plt.get_cmap('tab20')
for k, f in enumerate(fam_list):
    sub = per[per.family == f]
    ax.scatter(sub.d2 + 0.5, sub.ng + 0.5, s=9, color=pal(k % 20), label=f, linewidths=0, alpha=0.85)
ax.set_xscale('log')
ax.set_yscale('log')
lim = [0.4, max(per.ng.max(), per.d2.max()) * 1.5]
ax.plot(lim, lim, color=GREY, lw=0.6, ls='--')
ax.set_xlabel('modules called, DepMap 2D (+0.5)')
ax.set_ylabel('modules called, NextGen (+0.5)')
ax.set_title('Reactome / CORUM modules per RWE gene', loc='left')
ax.legend(frameon=False, loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, fontsize=5.3, markerscale=1.2, handletextpad=0.2, columnspacing=0.6)
panel(ax, 'a')
ax = fig.add_subplot(gs[1])
cnt = P.groupby('module').agg(ng=('call_NextGen', 'sum'), d2=('call_DepMap2D', 'sum'))
top = cnt.assign(tot=cnt.ng + cnt.d2).sort_values('tot', ascending=False).head(18).iloc[::-1]
y = np.arange(len(top))
ax.barh(y + 0.2, top.ng, 0.38, color=NG, label='NextGen')
ax.barh(y - 0.2, top.d2, 0.38, color=D2, label='DepMap 2D')
labels = [m.replace('Reactome: ', '').split(' R-HSA-')[0].replace('CORUM: ', 'CORUM · ')[:52] for m in top.index]
ax.set_yticks(y)
ax.set_yticklabels(labels, fontsize=5.8)
ax.set_xlabel('RWE genes with the module among their top-100 partners (FDR < 0.05)')
ax.set_title('Most frequent modules', loc='left')
ax.legend(frameon=False, loc='lower right')
panel(ax, 'b')
save(fig, 'fig3_modules')
