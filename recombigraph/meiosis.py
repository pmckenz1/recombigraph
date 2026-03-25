import numpy as np
from dataclasses import dataclass
from .ancestry import Homolog, Segment, Slot

####
# segment utilities
####

def get_segments_extent(
    h,
    start_loc: float,
    stop_loc: float,
    update_parent_id: bool = True,
) -> list[Segment]:
    """return copies of segments from h overlapping [start_loc, stop_loc)"""
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


def merge_adjacent_segments(segments: list[Segment]) -> list[Segment]:
    """merge adjacent segments only if both parent and founder IDs match"""
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


####
# slot-level crossover recording
####

@dataclass(frozen=True)
class SlotCrossover:
    pos: float
    partner_slot: int


@dataclass
class SlotRecord:
    """
    fixed slot used during PedigreeSim-style meiosis
    """
    slot_id: int
    homolog_id: int
    chromosome: str
    length: float
    events: list[SlotCrossover]

    def copy(self):
        return SlotRecord(
            slot_id=self.slot_id,
            homolog_id=self.homolog_id,
            chromosome=self.chromosome,
            length=self.length,
            segments=[seg.copy() for seg in self.segments],
            events=list(self.events),
        )



def make_slots(h0: Homolog, h1: Homolog) -> list[SlotRecord]:
    """produce the four fixed chromatids/slots for a diploid bivalent"""
    s0 = h0.to_slot(0)
    s1 = h0.to_slot(1)
    s2 = h1.to_slot(2)
    s3 = h1.to_slot(3)
    return [
        SlotRecord(s0.slot_id, s0.homolog_id, s0.chromosome, s0.length, []),
        SlotRecord(s1.slot_id, s1.homolog_id, s1.chromosome, s1.length, []),
        SlotRecord(s2.slot_id, s2.homolog_id, s2.chromosome, s2.length, []),
        SlotRecord(s3.slot_id, s3.homolog_id, s3.chromosome, s3.length, []),
    ]


def sample_nonsister_pair(
    rng: np.random.Generator | None = None,
) -> tuple[int, int]:
    """sample two nonsister chromatids. one from slots 0/1 and one from 2/3"""
    if rng is None:
        rng = np.random.default_rng()
    return int(rng.integers(0, 2)), int(rng.integers(2, 4))


####
# PedigreeSim-style recombination dists
####

def ran_exp(mean: float, rng: np.random.Generator) -> float:
    return rng.exponential(mean)


def dist_to_first_recomb(
    meandist: float,
    rng: np.random.Generator,
    chiasma_interference: bool = False,
) -> float:
    pos = ran_exp(meandist, rng)
    if chiasma_interference:
        raise NotImplementedError("chiasma_interference=True not yet implemented")
    return pos


def dist_to_next_recomb(
    meandist: float,
    rng: np.random.Generator,
    chiasma_interference: bool = False,
) -> float:
    if chiasma_interference:
        raise NotImplementedError("chiasma_interference=True not yet implemented")
    return ran_exp(meandist, rng)


####
# doCrossingOver
####

def all_chrom_recomb(slots: list[SlotRecord]) -> bool:
    """
    diploid PedigreeSim check...
    for a bivalent, at least one chromatid of the first homolog pair must
    have a recombination. then the other homolog pair necessarily does too.
    """
    return len(slots[0].events) > 0 or len(slots[1].events) > 0


def do_crossing_over(
    h0: Homolog,
    h1: Homolog,
    rng: np.random.Generator,
    chiasma_interference: bool = False,
    allow_no_recomb: bool = True,
) -> list[SlotRecord]:
    """
    PedigreeSim rewrite of doCrossingOver for a diploid bivalent
    """
    if h0.length != h1.length:
        raise ValueError("homologs should be same length")
    if h0.chromosome != h1.chromosome:
        raise ValueError("homologs should be from the same chromosome")

    # PedigreeSim RECOMBDIST = 50 cM for diploid bivalent
    recomb_dist = 50.0

    while True:
        slots = make_slots(h0, h1)

        pos = dist_to_first_recomb(
            recomb_dist,
            rng=rng,
            chiasma_interference=chiasma_interference,
        )

        while pos < h0.length:
            i, j = sample_nonsister_pair(rng)
            slots[i].events.append(SlotCrossover(pos=pos, partner_slot=j))
            slots[j].events.append(SlotCrossover(pos=pos, partner_slot=i))
            pos += dist_to_next_recomb(
                recomb_dist,
                rng=rng,
                chiasma_interference=chiasma_interference,
            )

        for slot in slots:
            slot.events.sort(key=lambda ev: ev.pos)

        if allow_no_recomb or all_chrom_recomb(slots):
            return slots


####
# slotsToPatterns
####

def _find_matching_event_index(
    events: list[SlotCrossover],
    pos: float,
    tol: float = 1e-12,
) -> int | None:
    """
    find the event in `events` occurring exactly at `pos`
    """
    for idx, ev in enumerate(events):
        if abs(ev.pos - pos) < tol:
            return idx
    return None


def slots_to_patterns(
    slots: list[SlotRecord],
) -> list[list[tuple[float, float, int]]]:
    """
    convert recorded slot connections into chromatid patterns

    returns one pattern per starting slot.
    Each pattern is a list of (left, right, chrom_index), where chrom_index
    is slot_id // 2, matching PedigreeSim's meshing of slot to chromosome
    """
    patterns: list[list[tuple[float, float, int]]] = []
    length = slots[0].length

    # slot_id -> sorted event list
    slot_events = {slot.slot_id: slot.events for slot in slots}

    for start_slot in range(len(slots)):
        pattern: list[tuple[float, float, int]] = []

        current_slot = start_slot
        current_left = 0.0
        current_event_idx = 0

        while True:
            events = slot_events[current_slot]

            if current_event_idx >= len(events):
                pattern.append((current_left, length, current_slot // 2))
                break

            ev = events[current_event_idx]

            # add segment from current slot up to thiss crossover
            pattern.append((current_left, ev.pos, current_slot // 2))

            # jump to partner slot at the same crossover
            partner_slot = ev.partner_slot
            partner_events = slot_events[partner_slot]
            partner_idx = _find_matching_event_index(partner_events, ev.pos)

            if partner_idx is None:
                raise ValueError(
                    f"Could not find reciprocal crossover at pos={ev.pos} "
                    f"between slots {current_slot} and {partner_slot}"
                )

            current_slot = partner_slot
            current_left = ev.pos
            current_event_idx = partner_idx + 1

        patterns.append(pattern)

    return patterns


####
# patternsToGametes / fillGamete
####

def get_centromere_sort_order(rng: np.random.Generator) -> list[int]:
    """
    PedigreeSim rewrite of bivalent tetrad ordering
    """
    if rng.random() < 0.5:
        return [0, 0, 1, 1]
    return [1, 1, 0, 0]


def pattern_value_at(
    pattern: list[tuple[float, float, int]],
    pos: float,
    tol: float = 1e-12,
) -> int:
    """
    return the chromosome identity (0/1) carried by the pattern at position pos
    """
    for left, right, chrom_idx in pattern:
        if left <= pos < right:
            return chrom_idx
    if abs(pos - pattern[-1][1]) < tol:
        return pattern[-1][2]
    raise ValueError(f"Position {pos} not covered by pattern")


def fill_gamete_from_pattern(
    pattern: list[tuple[float, float, int]],
    source_homologs: list[Homolog],
    output_slot_id: int,
) -> Slot:
    if not pattern:
        raise ValueError("pattern must not be empty")

    chromosome = source_homologs[0].chromosome
    length = source_homologs[0].length

    new_segments: list[Segment] = []
    for left, right, chrom_index in pattern:
        src = source_homologs[chrom_index]
        new_segments.extend(
            get_segments_extent(src, left, right, update_parent_id=True)
        )

    # previously merged here... but best not to to match pedigreeSim output
    return Slot(
        slot_id=output_slot_id,
        homolog_id=-1,
        chromosome=chromosome,
        length=length,
        segments=new_segments,
    )


def patterns_to_gametes(
    patterns: list[list[tuple[float, float, int]]],
    h0: Homolog,
    h1: Homolog,
    rng: np.random.Generator,
    centromere_pos: float | None = None,
) -> list[Slot]:
    if len(patterns) != 4:
        raise ValueError("Diploid bivalent should yield 4 patterns")

    if centromere_pos is None:
        centromere_pos = h0.length / 2.0

    centro = get_centromere_sort_order(rng)

    ordered = list(patterns)

    # PedigreeSim-style in-place sorting
    for s in range(len(ordered) - 1):
        if pattern_value_at(ordered[s], centromere_pos) != centro[s]:
            t = s + 1
            while pattern_value_at(ordered[t], centromere_pos) != centro[s]:
                t += 1
            ordered[s], ordered[t] = ordered[t], ordered[s]

    source_homologs = [h0, h1]
    gametes = [
        fill_gamete_from_pattern(ordered[s], source_homologs, output_slot_id=s)
        for s in range(4)
    ]
    return gametes


####
# high-level meiosis
####

def simulate_bivalent_meiosis(
    h0: Homolog,
    h1: Homolog,
    rng: np.random.Generator,
    chiasma_interference: bool = False,
    allow_no_recomb: bool = True,
    centromere_pos: float | None = None,
) -> list[Slot]:
    """
    PedigreeSim-style meiosis:
    three steps. doCrossingOver -> slotsToPatterns -> patternsToGametes
    """
    slots = do_crossing_over(
        h0,
        h1,
        rng=rng,
        chiasma_interference=chiasma_interference,
        allow_no_recomb=allow_no_recomb,
    )
    patterns = slots_to_patterns(slots)
    gametes = patterns_to_gametes(
        patterns,
        h0,
        h1,
        rng,
        centromere_pos=centromere_pos,
    )
    return gametes


def make_gamete(
    h0: Homolog,
    h1: Homolog,
    rng: np.random.Generator | None = None,
    chiasma_interference: bool = False,
    allow_no_recomb: bool = True,
) -> Slot:
    if rng is None:
        rng = np.random.default_rng()

    chromatids = simulate_bivalent_meiosis(
        h0,
        h1,
        rng=rng,
        chiasma_interference=chiasma_interference,
        allow_no_recomb=allow_no_recomb,
    )

    # matching PedigreeSim behavior: conception takes gamete 0.
    # this is still a random gamete because patterns_to_gametes randomizes
    # first-division pole assignment via get_centromere_sort_order()
    return chromatids[0]


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
        segments=[seg.copy() for seg in slot.segments],
    )