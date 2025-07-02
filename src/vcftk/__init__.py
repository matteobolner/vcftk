from ._version import __version__
from .main import (
    VCFClass,
    setup,
    set_pandas_display_options,
    compute_all_allele_frequencies,
    extract_vep_annotations,
)

__all__ = [
    "__version__",
    "VCFClass",
    "setup",
    "set_pandas_display_options",
    "compute_all_allele_frequencies",
    "extract_vep_annotations",
]
