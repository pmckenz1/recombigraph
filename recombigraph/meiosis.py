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
    """return copies of segments from h overlapping [start_loc, stop_loc)."""
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
    """merge adjacent segments only if both parent and founder ids match."""
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

    slot_id: 0,1 for first homolog sisters; 2,3 for second homolog sisters
    source_slot: original slot identity
    homolog_id/chromosome/length: copied from source homolog
    events: crossover positions + connected partner slot
    """
    slot_id: int
    homolog_id: int
    chromosome: str
    length: float
    segments: list[Segment]
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
        SlotRecord(s0.slot_id, s0.homolog_id, s0.chromosome, s0.length, [seg.copy() for seg in s0.segments], []),
        SlotRecord(s1.slot_id, s1.homolog_id, s1.chromosome, s1.length, [seg.copy() for seg in s1.segments], []),
        SlotRecord(s2.slot_id, s2.homolog_id, s2.chromosome, s2.length, [seg.copy() for seg in s2.segments], []),
        SlotRecord(s3.slot_id, s3.homolog_id, s3.chromosome, s3.length, [seg.copy() for seg in s3.segments], []),
    ]


def sample_nonsister_pair(
    rng: np.random.Generator | None = None,
) -> tuple[int, int]:
    """sample two nonsister chromatids: one from slots 0/1 and one from 2/3"""
    if rng is None:
        rng = np.random.default_rng()
    return int(rng.integers(0, 2)), int(rng.integers(2, 4))


####
# PedigreeSim-style recombination distances
####

def ran_exp(mean: float, rng: np.random.Generator) -> float:
    """Exponential draw with given mean."""
    return rng.exponential(mean)


def dist_to_first_recomb(
    meandist: float,
    rng: np.random.Generator,
    chiasma_interference: bool = False,
) -> float:
    """
    PedigreeSim-style distance to first recombination.

    Currently:
    - no interference: exponential
    - interference not yet implemented
    """
    pos = ran_exp(meandist, rng)
    if chiasma_interference:
        raise NotImplementedError("chiasma_interference=True not yet implemented")
    return pos


def dist_to_next_recomb(
    meandist: float,
    rng: np.random.Generator,
    chiasma_interference: bool = False,
) -> float:
    """
    PedigreeSim-style dist to next recombination

    for now:
    - no interference: exponential
    - interference not yet implemented
    """
    if chiasma_interference:
        raise NotImplementedError("chiasma_interference=True not yet implemented")
    return ran_exp(meandist, rng)


####
# doCrossingOver
####

def all_chrom_recomb(slots: list[SlotRecord]) -> bool:
    """
    replicating PedigreeSim's diploid check:
    for a bivalent, at least one chromatid from each homolog pair involved...
    since slots 0/1 are homolog 0 sisters and 2/3 are homolog 1 sisters,
    PedigreeSim effectively checks that one of the first pair has >=1 event
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
    PedigreeSim style doCrossingOver for diploid bivalent

    slots remain fixed. crossover events are recorded onto slot event lists
    """
    if h0.length != h1.length:
        raise ValueError("homologs should be same length")
    if h0.chromosome != h1.chromosome:
        raise ValueError("homologs should be from the same chromosome")

    # for a bivalent, PedigreeSim uses RECOMBDIST = 0.5 Morgan.
    # the length is in cM, so 0.5 Morgan = 50 cM
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

def slots_to_patterns(slots: list[SlotRecord]) -> list[list[tuple[float, float, int]]]:
    """
    convert recorded slot connections into chromatid patterns

    returns one pattern per starting slot
    
    each pattern is a list of (left, right, chromosome_index),
    where chromosome_index is slot_id // 2. this follows PedigreeSim's
    collapse from slot number to chromosome number
    """
    patterns = []

    # precompute for convenience
    slot_events = {slot.slot_id: slot.events for slot in slots}
    length = slots[0].length

    for start_slot in range(len(slots)):
        current_slot = start_slot
        current_pos = 0.0
        pattern = []

        # track how far along we are in each slot's event list
        next_idx = {sid: 0 for sid in slot_events}
        incoming_pos = None

        while True:
            events = slot_events[current_slot]
            idx = next_idx[current_slot]

            # if we just entered at a crossover position, skip that exact event
            if incoming_pos is not None:
                while idx < len(events) and abs(events[idx].pos - incoming_pos) > 1e-12:
                    idx += 1
                if idx < len(events):
                    idx += 1
                next_idx[current_slot] = idx
                incoming_pos = None

            if idx >= len(events):
                pattern.append((current_pos, length, current_slot // 2))
                break

            ev = events[idx]
            pattern.append((current_pos, ev.pos, current_slot // 2))
            current_pos = ev.pos

            prev_slot = current_slot
            current_slot = ev.partner_slot
            incoming_pos = ev.pos

            next_idx[prev_slot] += 1

        patterns.append(pattern)

    return patterns


####
# patternsToGametes / fillGamete
####

def get_centromere_sort_order(rng: np.random.Generator) -> list[int]:
    """
    PedigreeSim style bivalent tetrad ordering
    """
    if rng.random() < 0.5:
        return [0, 0, 1, 1]
    return [1, 1, 0, 0]


def fill_gamete_from_pattern(
    pattern: list[tuple[float, float, int]],
    source_homologs: list[Homolog],
    output_slot_id: int,
) -> Slot:
    """
    apply a pattern of chromosome indices (0/1) back to the original homologs
    """
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

    new_segments = merge_adjacent_segments(new_segments)

    # slot_id here is just a temporary tetrad index
    # homolog_id is assigned later in slot_to_homolog()
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
) -> list[Slot]:
    """
    sort patterns by centromeric chromosome identity and convert to gamete chromatids
    """
    if len(patterns) != 4:
        raise ValueError("Diploid bivalent should yield 4 patterns")

    centro = get_centromere_sort_order(rng)

    # sort patterns to match centromere order...
    ordered_patterns = list(patterns)
    for s in range(len(ordered_patterns) - 1):
        centromere_chrom = None
        for left, right, chrom_idx in ordered_patterns[s]:
            if left <= 0.0 < right or abs(left) < 1e-12:
                centromere_chrom = chrom_idx
                break
        if centromere_chrom is None:
            centromere_chrom = ordered_patterns[s][0][2]

        if centromere_chrom != centro[s]:
            t = s + 1
            while t < len(ordered_patterns):
                cent_t = ordered_patterns[t][0][2]
                if cent_t == centro[s]:
                    break
                t += 1
            if t == len(ordered_patterns):
                raise ValueError("Could not reorder patterns by centromere identity")
            ordered_patterns[s], ordered_patterns[t] = ordered_patterns[t], ordered_patterns[s]

    source_homologs = [h0, h1]
    gametes = [
        fill_gamete_from_pattern(ordered_patterns[s], source_homologs, output_slot_id=s)
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
) -> list[Slot]:
    """
    PedigreeSim-style meiosis.
    the full thing: doCrossingOver -> slotsToPatterns -> patternsToGametes
    """
    slots = do_crossing_over(
        h0,
        h1,
        rng=rng,
        chiasma_interference=chiasma_interference,
        allow_no_recomb=allow_no_recomb,
    )
    patterns = slots_to_patterns(slots)
    gametes = patterns_to_gametes(patterns, h0, h1, rng)
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

    # PedigreeSim treats gamete 0 as already random, but sampling one at random
    # is fine for now and keeps the existing interface.
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
        segments=[seg.copy() for seg in slot.segments],
    )