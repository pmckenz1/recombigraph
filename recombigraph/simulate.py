from dataclasses import dataclass
import numpy as np

from .ancestry import Segment, Homolog
from .meiosis import make_gamete, slot_to_homolog


@dataclass
class SimIndividual:
    individual_id: str
    time: int
    homologs_by_chromosome: dict[str, list[Homolog]]


def make_founder_individual(
    individual_id: str,
    chromosome: str,
    length: float,
    time: int,
    next_homolog_id: int,
) -> tuple[SimIndividual, int]:
    h0 = Homolog(
        homolog_id=next_homolog_id,
        chromosome=chromosome,
        individual_id=individual_id,
        time=time,
        length=length,
        segments=[Segment(0.0, length, None, next_homolog_id)],
    )
    next_homolog_id += 1

    h1 = Homolog(
        homolog_id=next_homolog_id,
        chromosome=chromosome,
        individual_id=individual_id,
        time=time,
        length=length,
        segments=[Segment(0.0, length, None, next_homolog_id)],
    )
    next_homolog_id += 1

    indiv = SimIndividual(
        individual_id=individual_id,
        time=time,
        homologs_by_chromosome={chromosome: [h0, h1]},
    )
    return indiv, next_homolog_id


def make_offspring_individual(
    child_id: str,
    parent1: SimIndividual,
    parent2: SimIndividual,
    chromosome: str,
    time: int,
    next_homolog_id: int,
    rng: np.random.Generator,
) -> tuple[SimIndividual, int]:
    p1_h0, p1_h1 = parent1.homologs_by_chromosome[chromosome]
    p2_h0, p2_h1 = parent2.homologs_by_chromosome[chromosome]

    gam1 = make_gamete(p1_h0, p1_h1, rng)
    gam2 = make_gamete(p2_h0, p2_h1, rng)

    child_h0 = slot_to_homolog(gam1, next_homolog_id, child_id, time)
    next_homolog_id += 1
    child_h1 = slot_to_homolog(gam2, next_homolog_id, child_id, time)
    next_homolog_id += 1

    indiv = SimIndividual(
        individual_id=child_id,
        time=time,
        homologs_by_chromosome={chromosome: [child_h0, child_h1]},
    )
    return indiv, next_homolog_id


def infer_generation(
    parent1: str,
    parent2: str,
    individuals: dict[str, SimIndividual],
) -> int:
    if parent1 == "NA" and parent2 == "NA":
        return 0
    return max(individuals[parent1].time, individuals[parent2].time) + 1


def simulate_pedigree(
    records: list[list[str]],
    chromosome: str = "A",
    length: float = 100.0,
    seed: int | None = None,
) -> dict[str, SimIndividual]:
    rng = np.random.default_rng(seed)
    individuals = {}
    next_homolog_id = 0

    for name, parent1, parent2 in records:
        time = infer_generation(parent1, parent2, individuals)

        if parent1 == "NA" and parent2 == "NA":
            indiv, next_homolog_id = make_founder_individual(
                name, chromosome, length, time, next_homolog_id
            )
        else:
            indiv, next_homolog_id = make_offspring_individual(
                name,
                individuals[parent1],
                individuals[parent2],
                chromosome,
                time,
                next_homolog_id,
                rng,
            )

        individuals[name] = indiv

    return individuals