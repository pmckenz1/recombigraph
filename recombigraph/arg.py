from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Iterator
from typing import Any, Iterable, TYPE_CHECKING

if TYPE_CHECKING:
    from .ancestry import Homolog, Segment
    from .simulate import SimulationResult

@dataclass(frozen=True)
class LocalForest:
    """ancestral forest over one genomic interval"""
    chromosome: str
    left: float
    right: float
    edges: frozenset[tuple[int, int]]
    sample_homolog_ids: tuple[int, ...]

    def nodes(self) -> tuple[int, ...]:
        """return all node ids present in this local forest"""
        out = set(self.sample_homolog_ids)
        for parent, child in self.edges:
            out.add(parent)
            out.add(child)
        return tuple(sorted(out))

    def roots(self) -> tuple[int, ...]:
        """return nodes with no parent within this local forest"""
        nodes = set(self.nodes())
        children = {child for _, child in self.edges}
        return tuple(sorted(nodes - children))

    def children_map(self) -> dict[int, tuple[int, ...]]:
        """return a parent to children mapping"""
        out: dict[int, list[int]] = {}
        for parent, child in self.edges:
            out.setdefault(parent, []).append(child)
        return {k: tuple(sorted(v)) for k, v in out.items()}

    def parent_map(self) -> dict[int, int]:
        """return a child to parent mapping"""
        return {child: parent for parent, child in self.edges}

    def span(self) -> float:
        """return the genomic width of this forest"""
        return self.right - self.left

    def is_empty(self) -> bool:
        """return whether this interval is empty"""
        return self.right <= self.left


@dataclass(frozen=True)
class LocalForestSequence:
    """ordered local forests across one chromosome"""
    chromosome: str
    forests: tuple[LocalForest, ...]

    def __iter__(self) -> Iterator[LocalForest]:
        """iterate over local forests"""
        return iter(self.forests)

    def __len__(self) -> int:
        return len(self.forests)

    def __getitem__(self, idx: int) -> LocalForest:
        """return one local forest by position"""
        return self.forests[idx]

    @property
    def left(self) -> float:
        """return the leftmost coordinate in the sequence"""
        return self.forests[0].left if self.forests else 0.0

    @property
    def right(self) -> float:
        """return the rightmost coordinate in the sequence"""
        return self.forests[-1].right if self.forests else 0.0

    def breakpoints(self) -> tuple[float, ...]:
        """return the sequence breakpoint coordinates"""
        if not self.forests:
            return tuple()
        vals = [self.forests[0].left]
        vals.extend(f.right for f in self.forests)
        return tuple(vals)

    def sample_homolog_ids(self) -> tuple[int, ...]:
        """return the sampled homolog ids for this sequence"""
        if not self.forests:
            return tuple()
        return self.forests[0].sample_homolog_ids

    def nodes(self) -> tuple[int, ...]:
        """return all node ids seen across the sequence"""
        out = set()
        for forest in self.forests:
            out.update(forest.nodes())
        return tuple(sorted(out))

def _build_homolog_lookup(result: "SimulationResult") -> dict[int, Any]:
    """build a homolog lookup keyed by homolog id"""
    out: dict[int, Any] = {}
    for individual in result.individuals.values():
        for homologs in individual.homologs_by_chromosome.values():
            for homolog in homologs:
                out[homolog.homolog_id] = homolog
    return out


def _segment_covering(homolog: "Homolog", position: float) -> "Segment | None":
    """return the segment covering a genomic position if one exists"""
    for seg in homolog.segments:
        if seg.left <= position < seg.right:
            return seg
    return None


def _local_forest_at_position(
    result: "SimulationResult",
    chromosome: str,
    sample_homolog_ids: Iterable[int],
    position: float,
) -> set[tuple[int, int]]:
    """return edges for the local forest at one genomic position"""
    homolog_lookup = _build_homolog_lookup(result)
    seen = set()
    edges: set[tuple[int, int]] = set()

    def visit(hid: int) -> None:
        if hid in seen:
            return
        seen.add(hid)

        h = homolog_lookup[hid]
        if h.chromosome != chromosome:
            return

        seg = _segment_covering(h, position)
        if seg is None:
            return

        parent_id = seg.parent_homolog_id
        if parent_id is not None:
            edges.add((parent_id, hid))
            visit(parent_id)

    for hid in sample_homolog_ids:
        visit(hid)

    return edges


def _ancestral_breakpoints(
    result: "SimulationResult",
    chromosome: str,
    sample_homolog_ids: Iterable[int],
) -> list[float]:
    """collect all ancestry breakpoint coordinates for sampled homologs"""
    homolog_lookup = _build_homolog_lookup(result)
    seen = set()
    breaks: set[float] = set()

    def visit(hid: int) -> None:
        if hid in seen:
            return
        seen.add(hid)

        h = homolog_lookup[hid]
        if h.chromosome != chromosome:
            return

        for seg in h.segments:
            breaks.add(seg.left)
            breaks.add(seg.right)
            if seg.parent_homolog_id is not None:
                visit(seg.parent_homolog_id)

    for hid in sample_homolog_ids:
        visit(hid)

    return sorted(breaks)


def _local_forests(
    result: "SimulationResult",
    chromosome: str,
    sample_homolog_ids: Iterable[int],
) -> list[tuple[float, float, set[tuple[int, int]]]]:
    """return local forests across all ancestry defined intervals"""
    homolog_lookup = _build_homolog_lookup(result)

    chrom_lengths = {
        h.length for h in homolog_lookup.values() if h.chromosome == chromosome
    }
    if not chrom_lengths:
        raise ValueError(f"Unknown chromosome: {chromosome!r}")
    if len(chrom_lengths) != 1:
        raise ValueError(f"Inconsistent lengths found for chromosome: {chromosome!r}")
    chrom_length = chrom_lengths.pop()

    breaks = set(_ancestral_breakpoints(result, chromosome, sample_homolog_ids))
    breaks.add(0.0)
    breaks.add(chrom_length)
    breaks = sorted(breaks)

    out: list[tuple[float, float, set[tuple[int, int]]]] = []
    for left, right in zip(breaks[:-1], breaks[1:]):
        if right <= left:
            continue
        mid = 0.5 * (left + right)
        edges = _local_forest_at_position(result, chromosome, sample_homolog_ids, mid)
        out.append((left, right, edges))

    return out
