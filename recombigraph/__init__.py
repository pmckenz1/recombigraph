"""recombigraph: pedigree-based recombination simulation with explicit homolog tracking."""

from typing import Any

from .ancestry import Segment, Homolog, Slot
from .pedigree import Pedigree, PedigreeRecord
from .genome import GenomeSpec, ChromosomeSpec
from .simulate import PedigreeModel, SimIndividual, SimulationResult
from .arg import LocalForest, LocalForestSequence
from .export import (
    forest_to_newicks,
    to_newick_records,
    to_dataframe,
    to_tskit,
)

# i'm using a lazy wrapper here to avoid auto-matplotlib import
def draw_pedigree_from_records(*args, **kwargs) -> Any:
    """lazy wrapper around the plotting helper"""
    from .visualize import draw_pedigree_from_records as _draw
    return _draw(*args, **kwargs)

__version__ = "0.1.0"

__all__ = [
    "Segment",
    "Homolog",
    "Slot",
    "PedigreeRecord",
    "Pedigree",
    "ChromosomeSpec",
    "GenomeSpec",
    "SimIndividual",
    "SimulationResult",
    "PedigreeModel",
    "draw_pedigree_from_records",
]
