from dataclasses import dataclass
import numpy as np

from .ancestry import Segment, Homolog
from .meiosis import make_gamete, slot_to_homolog
from .pedigree import Pedigree
from .genome import GenomeSpec
from .visualize import draw_pedigree_from_records


@dataclass
class SimIndividual:
    individual_id: str
    time: int
    homologs_by_chromosome: dict[str, list[Homolog]]


@dataclass
class SimulationResult:
    individuals: dict[str, SimIndividual]
    pedigree: Pedigree
    genome: GenomeSpec


def make_founder_individual(
    individual_id: str,
    genome: GenomeSpec,
    time: int,
    next_homolog_id: int,
) -> tuple[SimIndividual, int]:
    homologs_by_chromosome = {}

    for chrom in genome:
        h0 = Homolog(
            homolog_id=next_homolog_id,
            chromosome=chrom.name,
            individual_id=individual_id,
            time=time,
            length=chrom.length,
            segments=[Segment(0.0, chrom.length, None, next_homolog_id)],
        )
        next_homolog_id += 1

        h1 = Homolog(
            homolog_id=next_homolog_id,
            chromosome=chrom.name,
            individual_id=individual_id,
            time=time,
            length=chrom.length,
            segments=[Segment(0.0, chrom.length, None, next_homolog_id)],
        )
        next_homolog_id += 1

        homologs_by_chromosome[chrom.name] = [h0, h1]

    indiv = SimIndividual(
        individual_id=individual_id,
        time=time,
        homologs_by_chromosome=homologs_by_chromosome,
    )
    return indiv, next_homolog_id


def make_offspring_individual(
    child_id: str,
    parent1: SimIndividual,
    parent2: SimIndividual,
    genome: GenomeSpec,
    time: int,
    next_homolog_id: int,
    rng: np.random.Generator,
) -> tuple[SimIndividual, int]:
    homologs_by_chromosome = {}

    for chrom in genome:
        p1_h0, p1_h1 = parent1.homologs_by_chromosome[chrom.name]
        p2_h0, p2_h1 = parent2.homologs_by_chromosome[chrom.name]

        gam1 = make_gamete(p1_h0, p1_h1, rng)
        gam2 = make_gamete(p2_h0, p2_h1, rng)

        child_h0 = slot_to_homolog(gam1, next_homolog_id, child_id, time)
        next_homolog_id += 1

        child_h1 = slot_to_homolog(gam2, next_homolog_id, child_id, time)
        next_homolog_id += 1

        homologs_by_chromosome[chrom.name] = [child_h0, child_h1]

    indiv = SimIndividual(
        individual_id=child_id,
        time=time,
        homologs_by_chromosome=homologs_by_chromosome,
    )
    return indiv, next_homolog_id


def simulate_pedigree(
    pedigree: Pedigree,
    genome: GenomeSpec,
    seed: int | None = None,
) -> SimulationResult:
    rng = np.random.default_rng(seed)
    individuals = {}
    next_homolog_id = 0

    for rec in pedigree:
        time = pedigree.generation_map[rec.name]

        if rec.is_founder:
            indiv, next_homolog_id = make_founder_individual(
                rec.name,
                genome,
                time,
                next_homolog_id,
            )
        else:
            indiv, next_homolog_id = make_offspring_individual(
                rec.name,
                individuals[rec.parent1],
                individuals[rec.parent2],
                genome,
                time,
                next_homolog_id,
                rng,
            )

        individuals[rec.name] = indiv

    return SimulationResult(
        individuals=individuals,
        pedigree=pedigree,
        genome=genome,
    )


class PedigreeModel:
    def __init__(self, pedigree, chromosomes, seed: int | None = None):
        self.pedigree = pedigree if isinstance(pedigree, Pedigree) else Pedigree(pedigree)
        self.genome = chromosomes if isinstance(chromosomes, GenomeSpec) else GenomeSpec(chromosomes)
        self.seed = seed

    def draw_pedigree(self, **kwargs):
        records = [
            [rec.name, rec.parent1, rec.parent2]
            for rec in self.pedigree.records
        ]
        return draw_pedigree_from_records(records, **kwargs)

    def simulate(self) -> SimulationResult:
        return simulate_pedigree(
            pedigree=self.pedigree,
            genome=self.genome,
            seed=self.seed,
        )