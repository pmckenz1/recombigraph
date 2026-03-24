"""recombigraph: pedigree-based recombination simulation with explicit homolog tracking."""

from .ancestry import Segment, Homolog, Slot
from .pedigree import Pedigree, PedigreeRecord
from .genome import GenomeSpec, ChromosomeSpec
from .simulate import PedigreeModel, SimIndividual, SimulationResult
from .visualize import draw_pedigree_from_records

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