from dataclasses import dataclass, replace
from typing import Optional

@dataclass
class Segment:
    """ancestry segment on a homolog"""
    left: float
    right: float
    parent_homolog_id: Optional[int]
    founder_homolog_id: int

    def copy(self, **kwargs) -> "Segment":
        """return a copied segment with selected fields replaced"""
        return replace(self, **kwargs)

@dataclass
class Homolog:
    """homolog with explicit ancestry segments"""
    homolog_id: int
    chromosome: str
    individual_id: str
    time: int
    length: float
    segments: list[Segment]

    def to_slot(self, slot_id: int) -> "Slot":
        """convert a homolog into a meiosis slot"""
        return Slot(
            slot_id=slot_id,
            homolog_id=self.homolog_id,
            chromosome=self.chromosome,
            length=self.length,
            segments=[
                seg.copy(parent_homolog_id=self.homolog_id)
                for seg in self.segments
            ],
        )

@dataclass
class Slot:
    """temporary chromatid record used during meiosis"""
    slot_id: int
    homolog_id: int
    chromosome: str
    length: float
    segments: list[Segment]

    def copy(self, **kwargs) -> "Slot":
        """return a copied slot with selected fields replaced"""
        return replace(self, **kwargs)
