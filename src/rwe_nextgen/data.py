"""Data access with the conventions of the NextGen paper's public code (broadinstitute/DepMap-NextGen-Public).

* gene effect / dependency are the SCREEN-level Chronos matrices (`screen_gene_effect.csv`, `screen_gene_dependency.csv`); one row per screen,
  missing values are kept and every correlation is computed on pairwise-complete screens.
* gene keys are DepMap "SYMBOL (ENTREZ)" strings; omics matrices are re-keyed by Entrez ID so that their columns match the gene-effect columns.
* screen sets follow `src/generate.py` of the paper repository (see `nextgen_screens`, `traditional_lineage_matched_screens`, ...).
"""
import functools
import os
import re

import numpy as np
import pandas as pd

from . import paths

NEXTGEN_TYPES = ['2DO', '3DO', '3DN', '2DN']


def _find(name):
    for d in (paths.PORTAL, paths.FIGSHARE):
        f = os.path.join(d, name)
        if os.path.exists(f):
            return f
    raise FileNotFoundError(f'{name} not found in {paths.PORTAL} or {paths.FIGSHARE} (see data/README.md)')


def entrez_of(key):
    m = re.search(r'\((\d+)\)\s*$', str(key))
    return m.group(1) if m else None


def symbol_of(key):
    return str(key).split(' (')[0]


@functools.lru_cache(maxsize=None)
def screen_effect():
    """screens x genes, Chronos gene effect"""
    return pd.read_csv(_find('screen_gene_effect.csv'), index_col=0)


@functools.lru_cache(maxsize=None)
def screen_dependency():
    """screens x genes, probability of dependency"""
    return pd.read_csv(_find('screen_gene_dependency.csv'), index_col=0)


@functools.lru_cache(maxsize=None)
def screen_meta():
    sm = pd.read_csv(_find('screen_metadata.csv'), index_col=0)
    for c in ('IsNextGen', 'IsTraditional2D', 'IsAdherent2D', 'PassesQC', 'HasExprData', 'HasCNData'):
        if c in sm.columns:
            sm[c] = sm[c].astype(str).str.upper().eq('TRUE')
    sm['Library2'] = np.where(sm.Library.eq('Humagne-CD'), 'Cas12a (Humagne-CD)', 'Cas9')
    return sm


@functools.lru_cache(maxsize=None)
def tda_classes():
    """the paper's gene dependency classes (NextGen and organoid screens)"""
    return pd.read_csv(os.path.join(paths.FIGSHARE, 'tda_dependency_classification.csv'), index_col=0)


@functools.lru_cache(maxsize=None)
def _reference_keys():
    k = pd.read_csv(os.path.join(paths.FIGSHARE, 'tda_dependency_classification.csv'), usecols=['Dependency']).Dependency
    s = pd.Series(k.values, index=k.map(entrez_of).values)
    return s[~s.index.duplicated()]


def rekey_by_entrez(df):
    ref = _reference_keys()
    new = [ref.get(entrez_of(c), c) if entrez_of(c) is not None else c for c in df.columns]
    out = df.copy()
    out.columns = new
    return out.loc[:, ~out.columns.duplicated()]


@functools.lru_cache(maxsize=None)
def expression():
    """models x genes, log2(TPM + 1), NextGen + traditional models"""
    parts = [rekey_by_entrez(pd.read_csv(_find(f), index_col=0)) for f in ('traditional_expression.csv', 'next_gen_expression.csv')]
    full = pd.concat(parts, axis=0)
    return full[~full.index.duplicated(keep='last')]


def to_screens(model_matrix, screens):
    """expand a model-indexed matrix to the given screens (rows = ScreenID)"""
    sm = screen_meta().reindex(screens)
    out = model_matrix.reindex(sm.ModelID.values)
    out.index = list(screens)
    return out


# ------------------------------------------------------------------------------------------- screen sets (paper definitions)
def _in_ge(idx):
    ge = set(screen_effect().index)
    return [s for s in idx if s in ge]


@functools.lru_cache(maxsize=None)
def nextgen_screens():
    """NextGen screens: ScreenType in {2DO, 3DO, 3DN, 2DN}, IsNextGen, lineage != 'Other'  (147)"""
    sm = screen_meta()
    return tuple(_in_ge(sm[sm.ScreenType.isin(NEXTGEN_TYPES) & sm.IsNextGen & sm.OncotreeLineage.ne('Other')].index))


@functools.lru_cache(maxsize=None)
def nextgen_lineages():
    sm = screen_meta()
    return tuple(sm[sm.IsNextGen].OncotreeLineage.value_counts().index)


@functools.lru_cache(maxsize=None)
def traditional_lineage_matched_screens():
    """'adherent_next_gen_lineage_matched_screens' of the paper: IsTraditional2D screens of the NextGen lineages (both libraries)"""
    sm = screen_meta()
    return tuple(_in_ge(sm[sm.IsTraditional2D & sm.OncotreeLineage.isin(nextgen_lineages())].index))


@functools.lru_cache(maxsize=None)
def organoid_3d_screens():
    sm = screen_meta()
    return tuple(_in_ge(sm[sm.IsNextGen & sm.ScreenType.eq('3DO') & sm.OncotreeLineage.ne('Other')].index))


@functools.lru_cache(maxsize=None)
def adherent_organoid_lineage_matched_screens():
    sm = screen_meta()
    lins = sm.loc[list(organoid_3d_screens()), 'OncotreeLineage'].unique()
    return tuple(_in_ge(sm[sm.IsAdherent2D & sm.OncotreeLineage.isin(lins)].index))


DATASETS = {'NextGen': nextgen_screens, 'DepMap2D': traditional_lineage_matched_screens}


# ------------------------------------------------------------------------------------------- curated inputs
@functools.lru_cache(maxsize=None)
def rwe_genes():
    """curated RWE gene list (one row per gene; `key` = DepMap gene key)"""
    return pd.read_csv(os.path.join(paths.CURATED, 'rwe_genes.tsv'), sep='\t')


@functools.lru_cache(maxsize=None)
def control_genes():
    return pd.read_csv(os.path.join(paths.CURATED, 'rwe_control_genes.tsv'), sep='\t')
