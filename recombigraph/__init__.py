"""recombigraph: pedigree-based recombination simulation with explicit homolog tracking."""

from .ancestry import Segment, Homolog
from .visualize import draw_pedigree_from_records

__version__ = "0.1.0"

__all__ = [
    "Segment",
    "Homolog",
    "draw_pedigree_from_records",
]
