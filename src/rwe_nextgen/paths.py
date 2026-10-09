"""Project paths. Everything is resolved relative to the repository root, so the pipeline runs from any checkout."""
import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

RAW = os.path.join(ROOT, 'data', 'raw')
PORTAL = os.path.join(RAW, 'depmap_portal_nextgen')          # DepMap portal file set "NextGen Model Manuscript 2026"
FIGSHARE = os.path.join(RAW, 'figshare_nextgen_29472362')     # paper's processed tables (figshare 10.6084/m9.figshare.29472362)
GENESETS = os.path.join(RAW, 'genesets')                      # Enrichr GO BP/CC 2023, Reactome 2022, KEGG 2021, CORUM

CURATED = os.path.join(ROOT, 'data', 'curated')               # hand-curated inputs shipped with the repository
RESULTS = os.path.join(ROOT, 'results')
TABLES = os.path.join(RESULTS, 'tables')
FIGURES = os.path.join(RESULTS, 'figures')
INTERMEDIATE = os.path.join(RESULTS, 'intermediate')         # large machine-readable intermediates (re-creatable)
WEB = os.path.join(ROOT, 'web')


def out(*parts):
    """Path under results/, creating the parent directory."""
    p = os.path.join(RESULTS, *parts)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    return p
