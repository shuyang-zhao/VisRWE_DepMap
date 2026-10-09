"""Shared bootstrap for pipeline steps: puts src/ on the import path and provides small I/O helpers."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'src'))

import pandas as pd  # noqa: E402

from rwe_nextgen import data, paths, stats  # noqa: E402,F401


def save(df, *parts, index=False, fmt='%.6g'):
    """write a TSV (or .tsv.gz) under results/ and report it"""
    p = paths.out(*parts)
    df.to_csv(p, sep='\t', index=index, float_format=fmt)
    print(f'  wrote {os.path.relpath(p, paths.ROOT)} ({len(df):,} rows)')
    return p


def load(*parts, **kw):
    return pd.read_csv(os.path.join(paths.RESULTS, *parts), sep='\t', **kw)


def banner(title):
    print('=' * 100)
    print(title)
    print('=' * 100)
