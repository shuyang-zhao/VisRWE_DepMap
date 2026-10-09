"""Step 6 - functional modules of each RWE gene's co-dependency profile, NextGen vs DepMap 2D (and Cas12a 2D as a library check).

Two levels, defined the same way for every gene (no hand-made modules):
  categories : 16 broad functional categories assigned to each gene by regular expressions on its GO BP / GO CC / Reactome / KEGG / CORUM
               terms (first match in a fixed order). Used for node colours and the overview heat map.
  pathways   : every Reactome 2022 pathway and CORUM complex with 5-150 members in the dependency-gene family (step 4); identical member
               sets are kept once.
Test       : the paper's over-representation convention (DepMap-NextGen-Public run_rank_hypergeometric): the 100 most positively
             correlated partners of each RWE gene (universe = dependency-gene family) are tested against every module with a one-sided
             hypergeometric test; mean r with the module members is reported as the effect size.
FDR        : Benjamini-Hochberg per dataset and level over the (RWE gene x module) table; call FDR < 0.05.
Datasets   : NextGen, DepMap 2D, and the 36 traditional 2D screens made with the NextGen Cas12a library (library check).
Outputs: results/tables/gene_categories.tsv, module_enrichment_categories.tsv, module_enrichment_pathways.tsv (called in any dataset),
         results/intermediate/module_enrichment_pathways_all.tsv.gz
"""
import os
import re

import numpy as np
import pandas as pd
from _common import banner, data, paths, save, stats
from scipy import stats as ss

TOP_N = 100
MIN_SIZE, MAX_SIZE = 5, 150
UNANNOTATED = 'Unannotated / other'
CATEGORIES = [
    ('Ribosome biogenesis / rRNA processing', r'ribosom\w* (large|small)? ?(subunit )?biogenesis|rRNA (processing|metabolic|modification|methylation|pseudouridin)|maturation of (LSU|SSU|5\.8S)|'
                                              r'ribosome biogenesis|snoRNA|small-subunit processome|preribosome|ribonucleoprotein complex biogenesis|nucleolar|box [CH]/[DA]'),
    ('tRNA processing / modification / aminoacylation', r'tRNA|aminoacyl'),
    ('Mitochondrial translation / OXPHOS', r'mitochondrial (translation|gene expression|ribosom|respiratory|electron|ATP synthesis|complex I|large ribosomal|small ribosomal)|oxidative phosphorylation|'
                                           r'respiratory (electron|chain)|NADH dehydrogenase|cytochrome[- ]c oxidase|ATP synthase'),
    ('Mitochondria (other)', r'mitochondri'),
    ('Cytosolic translation / ribosome', r'translation(al)? (initiation|elongation|termination)|eukaryotic translation|cytosolic (large|small) ribosomal|ribosom'),
    ('mRNA splicing / processing / export', r'spliceosom|mRNA splicing|mRNA (processing|3.-end|export|surveillance)|snRNP|mRNA cleavage|polyadenylation|nuclear pore|nuclear mRNA|RNA export'),
    ('Transcription machinery', r'RNA polymerase|transcription (initiation|elongation|factor II|factor III)|Mediator|TFII|general transcription|SAGA|Integrator|transcription by|Paf1'),
    ('DNA replication / repair', r'DNA (replication|repair|damage)|homologous recombination|double-strand|nucleotide excision|base excision|mismatch repair|replication fork|Fanconi'),
    ('Cell cycle / mitosis / chromosome', r'mitotic|cell cycle|kinetochore|chromosome segregation|centromere|spindle|cohesin|condensin|centrosome|cytokinesis|anaphase|APC/C|anaphase-promoting'),
    ('Chromatin / epigenetic', r'chromatin|histone|nucleosome|SWI/SNF|NuRD|polycomb|PRC2|epigenetic'),
    ('Proteostasis (proteasome / ubiquitin / folding)', r'proteasome|ubiquitin|protein folding|chaperon|unfolded protein|ER-associated degradation|neddylation|SUMO'),
    ('ER / Golgi / vesicle trafficking', r'endoplasmic reticulum|Golgi|vesicle|COPI|COPII|SNARE|protein transport|GPI|endosome|lysosom|glycosylation'),
    ('Lipid / sterol metabolism', r'cholesterol|sterol|fatty acid|lipid|isoprenoid|steroid|sphingo|phospholipid'),
    ('Other metabolism', r'metabolic process|biosynthetic process|purine|pyrimidine|nucleotide|one-carbon|folate|glycolysis|tricarboxylic|amino acid'),
    ('Signaling / receptor / immune', r'signaling|receptor|cytokine|immune|T cell|G protein|GPCR|MAPK|kinase|phosphatase|Wnt|Notch|interferon|chemokine'),
    ('Cytoskeleton / adhesion / ECM', r'actin|microtubule|cytoskelet|adhesion|extracellular matrix|integrin|cilium|motility|myosin|kinesin'),
]

banner('6. functional modules')
ge, sm = data.screen_effect(), data.screen_meta()


def read_gmt(name):
    for line in open(os.path.join(paths.GENESETS, name + '.gmt')):
        f = line.rstrip('\n').split('\t')
        yield f[0], [g for g in f[2:] if g]


# ---------------------------------------------------------------------------------------------- level 1: categories for every gene
terms = {}
for lib in ('GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'Reactome_2022', 'KEGG_2021_Human', 'CORUM'):
    for name, genes in read_gmt(lib):
        if len(genes) <= 600:
            for g in genes:
                terms.setdefault(g, []).append(name)
rx = [(c, re.compile(p, re.I)) for c, p in CATEGORIES]
cat = pd.DataFrame({'key': ge.columns, 'symbol': [data.symbol_of(k) for k in ge.columns]})
cat['category'] = [next((c for c, r in rx if any(r.search(x) for x in terms.get(s, []))), UNANNOTATED) for s in cat.symbol]
cat['n_terms'] = cat.symbol.map(lambda s: len(terms.get(s, [])))
save(cat, 'tables', 'gene_categories.tsv')
print(cat.category.value_counts().to_string())

# ---------------------------------------------------------------------------------------------- module definitions on the family
z0 = np.load(os.path.join(paths.INTERMEDIATE, 'codep_NextGen.npz'), allow_pickle=True)
fam = [str(g) for g in z0['genes']]
fsym = np.array([data.symbol_of(k) for k in fam])
fpos = {s: i for i, s in enumerate(fsym)}
catmap = cat.set_index('key').category
modules = {}
for c in [c for c, _ in CATEGORIES]:
    modules[('category', c)] = np.where(catmap.reindex(fam).values == c)[0]
seen = set()
for lib, prefix in (('Reactome_2022', 'Reactome'), ('CORUM', 'CORUM')):
    for name, genes in read_gmt(lib):
        idx = np.array(sorted({fpos[g] for g in genes if g in fpos}))
        if MIN_SIZE <= len(idx) <= MAX_SIZE and tuple(idx) not in seen:
            seen.add(tuple(idx))
            modules[('pathway', f'{prefix}: {name}')] = idx
print(f'modules: {sum(1 for k in modules if k[0] == "category")} categories, {sum(1 for k in modules if k[0] == "pathway")} pathways / complexes '
      f'({MIN_SIZE}-{MAX_SIZE} family members)')

# ---------------------------------------------------------------------------------------------- full correlation matrices in three datasets
rwe = set(data.rwe_genes().key.dropna())
ti = np.array([i for i, k in enumerate(fam) if k in rwe])
cas12a_2d = [x for x in sm[sm.IsTraditional2D & sm.Library.eq('Humagne-CD')].index if x in ge.index]
matrices = {d: (np.load(os.path.join(paths.INTERMEDIATE, f'codep_{d}.npz'))['r'].astype(np.float32), list(data.DATASETS[d]())) for d in ('NextGen', 'DepMap2D')}
R12, _ = stats.pairwise_r(ge.loc[cas12a_2d, fam].values)
np.fill_diagonal(R12, np.nan)
matrices['Cas12a2D'] = (R12, cas12a_2d)


def module_ora(R):
    """top-N partners of every RWE gene vs every module: overlap, odds ratio, hypergeometric P (upper tail), mean r with members"""
    Rt = R[ti].astype(np.float64)
    valid = ~np.isnan(Rt).all(0)
    top = np.zeros_like(Rt, dtype=bool)
    order = np.argsort(-np.nan_to_num(Rt, nan=-np.inf), axis=1)[:, :TOP_N]
    np.put_along_axis(top, order, True, axis=1)
    universe = int(valid.sum())
    out = {}
    for key, mem in modules.items():
        mem = mem[valid[mem]]
        k = top[:, mem].sum(1)
        n = len(mem)
        p = ss.hypergeom.sf(k - 1, universe, n, TOP_N)
        with np.errstate(divide='ignore', invalid='ignore'):
            odds = (k * (universe - n - TOP_N + k)) / ((n - k) * (TOP_N - k))
        out[key] = (np.nanmean(Rt[:, mem], axis=1), k, odds, p)
    return out


res = {d: module_ora(matrices[d][0]) for d in matrices}
rows = []
for (level, name), mem in modules.items():
    for j, i in enumerate(ti):
        row = dict(level=level, module=name, n_members=len(mem), key=fam[i], symbol=fsym[i])
        for d in res:
            mr, k, odds, p = res[d][(level, name)]
            row[f'mean_r_{d}'], row[f'overlap_{d}'], row[f'odds_{d}'], row[f'p_{d}'] = mr[j], int(k[j]), odds[j], p[j]
        rows.append(row)
M = pd.DataFrame(rows)
for d in res:
    M[f'FDR_{d}'] = np.nan
    for level in ('category', 'pathway'):
        k = M.level == level
        M.loc[k, f'FDR_{d}'] = stats.bh(M.loc[k, f'p_{d}'].values)
    M[f'call_{d}'] = M[f'FDR_{d}'] < 0.05
M['class'] = np.select([M.call_NextGen & M.call_DepMap2D, M.call_NextGen, M.call_DepMap2D], ['both', 'NextGen only', 'DepMap 2D only'], '')
cols = ['level', 'module', 'n_members', 'symbol', 'key', 'class'] + [c for c in M.columns if c.endswith(('_NextGen', '_DepMap2D', '_Cas12a2D'))]
M = M[cols]
save(M[M.level == 'category'].drop(columns='level'), 'tables', 'module_enrichment_categories.tsv')
P = M[M.level == 'pathway'].drop(columns='level')
save(P, 'intermediate', 'module_enrichment_pathways_all.tsv.gz')
save(P[P.call_NextGen | P.call_DepMap2D], 'tables', 'module_enrichment_pathways.tsv')
for level in ('category', 'pathway'):
    x = M[M.level == level]
    print(f'{level}: RWE gene x module calls', x['class'].value_counts().to_dict(), '| called in Cas12a 2D:', int(x.call_Cas12a2D.sum()))
