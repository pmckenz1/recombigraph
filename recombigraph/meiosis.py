from typing import Optional
import numpy as np
from .ancestry import Homolog, Segment, Slot

def get_segments_extent(
    h,
    start_loc: float,
    stop_loc: float,
    update_parent_id: bool = True,
) -> list[Segment]:
    """Return copies of segments from homolog h overlapping [start_loc, stop_loc)."""
    if start_loc >= stop_loc:
        raise ValueError("start_loc must be less than stop_loc")

    segs = []
    for seg in h.segments:
        if seg.left < stop_loc and start_loc < seg.right:
            trimmed_seg = seg.copy()
            trimmed_seg.left = max(trimmed_seg.left, start_loc)
            trimmed_seg.right = min(trimmed_seg.right, stop_loc)
            if update_parent_id:
                trimmed_seg.parent_homolog_id = h.homolog_id
            segs.append(trimmed_seg)

    return segs


def breakpoints_to_intervals(breakpoints: list[float], length: float) -> list[tuple[float, float]]:
    """Convert breakpoints into half-open intervals [left, right)."""
    bps = sorted(set(bp for bp in breakpoints if 0 < bp < length))
    edges = [0.0] + bps + [length]
    return [(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]


def merge_adjacent_segments(
    segments: list[Segment],
) -> list[Segment]:
    if not segments:
        return []

    segments = sorted(segments, key=lambda s: s.left)
    merged = [segments[0].copy()]
    for seg in segments[1:]:
        last = merged[-1]
        if (
            abs(last.right - seg.left) < 1e-12
            and last.parent_homolog_id == seg.parent_homolog_id
            and last.founder_homolog_id == seg.founder_homolog_id
        ):
            last.right = seg.right
        else:
            merged.append(seg.copy())

    return merged


def recombine_two_homologs(
    h0: Homolog,
    h1: Homolog,
    breakpoints: list[float],
    start_phase: int = 0,
    merge_adjacent: bool = True,
    homolog_id: Optional[int] = None,
    individual_id: Optional[str] = None,
    time: Optional[int] = None,
) -> Homolog:
    """DEPRECATED: Deterministically recombine two homologs and return a new recombinant Homolog."""
    if h0.length != h1.length:
        raise ValueError("homologs should be same length")
    if h0.chromosome != h1.chromosome:
        raise ValueError("homologs should be from the same chromosome")
    if start_phase not in (0, 1):
        raise ValueError("start_phase must be 0 or 1")

    intervals = breakpoints_to_intervals(breakpoints, h0.length)
    new_segments = []
    phase = start_phase

    for start, stop in intervals:
        if phase == 0:
            new_segments.extend(get_segments_extent(h0, start, stop))
        else:
            new_segments.extend(get_segments_extent(h1, start, stop))
        phase = 1 - phase

    if merge_adjacent:
        new_segments = merge_adjacent_segments(new_segments)

    return Homolog(
        homolog_id=homolog_id,
        chromosome=h0.chromosome,
        individual_id=individual_id,
        length=h0.length,
        time=time,
        segments=new_segments,
    )

def make_slots(h0: Homolog, h1: Homolog) -> list[Slot]:
    """produce the four chromatids a la PedigreeSim"""
    return [h0.to_slot(0),h0.to_slot(1),h1.to_slot(2),h1.to_slot(3)]

def sample_nonsister_pair(
    rng: np.random.Generator | None = None,
) -> tuple[int, int]:
    if rng is None:
        rng = np.random.default_rng()
    return int(rng.integers(0, 2)), int(rng.integers(2, 4))

def sample_breakpoints_haldane(length: float, rng: np.random.Generator | None = None) -> list[float]:
    """Sample crossover breakpoints under a simple Haldane model."""
    if rng is None:
        rng = np.random.default_rng()

    nxo = rng.poisson(length / 100.0)
    if nxo == 0:
        return []

    bps = sorted(rng.uniform(0.0, length, size=nxo).tolist())
    return bps
    
def crossover_slots(slot_a: Slot, slot_b: Slot, pos: float) -> tuple[Slot, Slot]:
    """Perform one crossover between two chromatids at position pos."""
    if slot_a.length != slot_b.length:
        raise ValueError("slot_a and slot_b must have the same length")
    if slot_a.chromosome != slot_b.chromosome:
        raise ValueError("slot_a and slot_b must be on the same chromosome")
    if not (0.0 < pos < slot_a.length):
        raise ValueError("pos must lie strictly within the chromosome")

    left_a = get_segments_extent(slot_a, 0.0, pos, update_parent_id=False)
    right_a = get_segments_extent(slot_a, pos, slot_a.length, update_parent_id=False)

    left_b = get_segments_extent(slot_b, 0.0, pos, update_parent_id=False)
    right_b = get_segments_extent(slot_b, pos, slot_b.length, update_parent_id=False)

    new_slot_a = Slot(
        slot_id=slot_a.slot_id,
        homolog_id=slot_a.homolog_id,
        chromosome=slot_a.chromosome,
        length=slot_a.length,
        segments=merge_adjacent_segments(left_a + right_b),
    )

    new_slot_b = Slot(
        slot_id=slot_b.slot_id,
        homolog_id=slot_b.homolog_id,
        chromosome=slot_b.chromosome,
        length=slot_b.length,
        segments=merge_adjacent_segments(left_b + right_a),
    )

    return new_slot_a, new_slot_b

def simulate_bivalent_meiosis(
    h0: Homolog,
    h1: Homolog,
    rng: np.random.Generator,
) -> list[Slot]:
    chromatids = make_slots(h0,h1)
    crossover_breakpoints = sample_breakpoints_haldane(h0.length, rng)
    events = [(brkpt, sample_nonsister_pair(rng)) for brkpt in crossover_breakpoints]
        
    for pos, (i, j) in events:
        sl0, sl1 = crossover_slots(chromatids[i], chromatids[j], pos)
        chromatids[sl0.slot_id] = sl0
        chromatids[sl1.slot_id] = sl1
    return chromatids

def make_gamete(
    h0: Homolog,
    h1: Homolog,
    rng: np.random.Generator | None = None,
) -> Slot:
    if rng is None:
        rng = np.random.default_rng()

    chromatids = simulate_bivalent_meiosis(h0, h1, rng)
    idx = rng.integers(0, len(chromatids))
    return chromatids[idx]

def slot_to_homolog(
    slot: Slot,
    homolog_id: int,
    individual_id: str,
    time: int,
) -> Homolog:
    return Homolog(
        homolog_id=homolog_id,
        chromosome=slot.chromosome,
        individual_id=individual_id,
        length=slot.length,
        time=time,
        segments=[seg.copy() for seg in slot.segments]
    )