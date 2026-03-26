from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, Iterable, Optional, Tuple, Union
import numpy as np

from .ancestry import Segment, Homolog
from .meiosis import make_gamete, slot_to_homolog
from .pedigree import Pedigree
from .genome import GenomeSpec
from .arg import LocalForest, LocalForestSequence, _local_forests


@dataclass
class SimIndividual:
    """simulated diploid individual with homologs by chromosome"""
    individual_id: str
    time: int
    homologs_by_chromosome: dict[str, list[Homolog]]


@dataclass
class SimulationResult:
    """simulation output with individuals pedigree and genome"""
    individuals: dict[str, SimIndividual]
    pedigree: Pedigree
    genome: GenomeSpec

    def __repr__(self) -> str:
        chromosome_names = [chrom.name for chrom in self.genome]
        times = sorted({ind.time for ind in self.individuals.values()})
        final_ids = self.final_generation_ids()

        if times:
            gen_text = f"{times[0]}..{times[-1]}"
        else:
            gen_text = "NA"

        return (
            f"SimulationResult("
            f"n_individuals={len(self.individuals)}, "
            f"generations={gen_text}, "
            f"chromosomes={chromosome_names}, "
            f"final_generation_ids={final_ids}"
            f")"
        )

    def summary(self) -> str:
        chromosome_names = [chrom.name for chrom in self.genome]
        times = sorted({ind.time for ind in self.individuals.values()})
        n_homologs = sum(
            len(homologs)
            for ind in self.individuals.values()
            for homologs in ind.homologs_by_chromosome.values()
        )
        n_segments = sum(
            len(h.segments)
            for ind in self.individuals.values()
            for homologs in ind.homologs_by_chromosome.values()
            for h in homologs
        )

        lines = [
            "SimulationResult summary",
            f"  individuals           : {len(self.individuals)}",
            f"  generations           : {times[0]}..{times[-1]}",
            f"  chromosomes           : {chromosome_names}",
            f"  final generation ids  : {self.final_generation_ids()}",
            f"  total homologs        : {n_homologs}",
            f"  total segments        : {n_segments}",
        ]
        return "\n".join(lines)

    def individuals_dataframe(self):
        import pandas as pd

        rows = []
        for ind in self.individuals.values():
            rows.append(
                {
                    "individual_id": ind.individual_id,
                    "time": ind.time,
                    "n_chromosomes": len(ind.homologs_by_chromosome),
                    "n_homologs": sum(len(v) for v in ind.homologs_by_chromosome.values()),
                }
            )
        return pd.DataFrame(rows).sort_values(["time", "individual_id"]).reset_index(drop=True)

    def homologs_dataframe(self):
        import pandas as pd

        rows = []
        for ind in self.individuals.values():
            for chrom, homologs in ind.homologs_by_chromosome.items():
                for h in homologs:
                    rows.append(
                        {
                            "individual_id": ind.individual_id,
                            "time": ind.time,
                            "chromosome": chrom,
                            "homolog_id": h.homolog_id,
                            "length": h.length,
                            "n_segments": len(h.segments),
                        }
                    )
        return pd.DataFrame(rows).sort_values(
            ["time", "individual_id", "chromosome", "homolog_id"]
        ).reset_index(drop=True)

    def segments_dataframe(self):
        import pandas as pd

        rows = []
        for ind in self.individuals.values():
            for chrom, homologs in ind.homologs_by_chromosome.items():
                for h in homologs:
                    for seg in h.segments:
                        rows.append(
                            {
                                "individual_id": ind.individual_id,
                                "time": ind.time,
                                "chromosome": chrom,
                                "homolog_id": h.homolog_id,
                                "left": seg.left,
                                "right": seg.right,
                                "span": seg.right - seg.left,
                                "parent_homolog_id": seg.parent_homolog_id,
                                "founder_homolog_id": seg.founder_homolog_id,
                            }
                        )
        return pd.DataFrame(rows).sort_values(
            ["time", "individual_id", "chromosome", "homolog_id", "left"]
        ).reset_index(drop=True)

    def final_generation_time(self) -> int:
        """Return the maximum generation time present in the simulation."""
        return max(ind.time for ind in self.individuals.values())

    def final_generation_individuals(self) -> list[SimIndividual]:
        """Return all individuals in the final generation."""
        t = self.final_generation_time()
        return [ind for ind in self.individuals.values() if ind.time == t]

    def final_generation_ids(self) -> list[str]:
        """Return individual_ids for the final generation."""
        return [ind.individual_id for ind in self.final_generation_individuals()]

    def final_generation_homolog_ids(self, chromosome: str | None = None) -> list[int]:
        """
        Return homolog_ids for all individuals in the final generation.

        If chromosome is provided, restrict to that chromosome.
        """
        inds = self.final_generation_individuals()

        if chromosome is not None:
            valid_chromosomes = {chrom.name for chrom in self.genome}
            if chromosome not in valid_chromosomes:
                raise ValueError(
                    f"Unknown chromosome {chromosome!r}. "
                    f"Valid chromosomes: {sorted(valid_chromosomes)}"
                )

        out = []
        for ind in inds:
            if chromosome is None:
                for homologs in ind.homologs_by_chromosome.values():
                    for h in homologs:
                        out.append(h.homolog_id)
            else:
                for h in ind.homologs_by_chromosome[chromosome]:
                    out.append(h.homolog_id)

        return out

    def local_forests(
        self,
        chromosome: str,
        sample_homolog_ids: Iterable[int],
    ) -> LocalForestSequence:
        """build local ancestry forests for sampled homologs on one chromosome"""
        sample_homolog_ids = tuple(sample_homolog_ids)

        # validate chromosome exists
        valid_chromosomes = {chrom.name for chrom in self.genome}
        if chromosome not in valid_chromosomes:
            raise ValueError(
                f"Unknown chromosome {chromosome!r}. "
                f"Valid chromosomes: {sorted(valid_chromosomes)}"
            )

        # build lookup
        homolog_lookup = {}
        for individual in self.individuals.values():
            for homologs in individual.homologs_by_chromosome.values():
                for homolog in homologs:
                    homolog_lookup[homolog.homolog_id] = homolog

        # validate sample IDs
        for hid in sample_homolog_ids:
            if hid not in homolog_lookup:
                raise ValueError(f"Unknown homolog_id: {hid}")
            if homolog_lookup[hid].chromosome != chromosome:
                raise ValueError(
                    f"Homolog {hid} belongs to chromosome "
                    f"{homolog_lookup[hid].chromosome!r}, not {chromosome!r}."
                )

        records = _local_forests(self, chromosome, sample_homolog_ids)
        forests = tuple(
            LocalForest(
                chromosome=chromosome,
                left=left,
                right=right,
                edges=frozenset(edges),
                sample_homolog_ids=sample_homolog_ids,
            )
            for left, right, edges in records
            if right > left
        )
        return LocalForestSequence(chromosome=chromosome, forests=forests)


def make_founder_individual(
    individual_id: str,
    genome: GenomeSpec,
    time: int,
    next_homolog_id: int,
) -> tuple[SimIndividual, int]:
    """create a founder with one pair of founder homologs per chromosome"""
    homologs_by_chromosome: dict[str, list[Homolog]] = {}

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
    """create an offspring by drawing one gamete from each parent"""
    homologs_by_chromosome: dict[str, list[Homolog]] = {}

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
    seed: Optional[int] = None,
) -> SimulationResult:
    """simulate homolog transmission through a pedigree"""
    rng = np.random.default_rng(seed)
    individuals: dict[str, SimIndividual] = {}
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
    """convenience wrapper around pedigree simulation"""

    def __init__(
        self,
        pedigree: Union[Pedigree, Iterable[Tuple[str, Optional[str], Optional[str]]]],
        chromosomes: Union[GenomeSpec, Dict[str, float], Iterable[Tuple[str, float]]],
        seed: Optional[int] = None,
    ) -> None:
        """normalize inputs and store simulation settings"""
        self.pedigree = pedigree if isinstance(pedigree, Pedigree) else Pedigree(pedigree)
        self.genome = chromosomes if isinstance(chromosomes, GenomeSpec) else GenomeSpec(chromosomes)
        self.seed = seed

    def draw_pedigree(self, **kwargs) -> Any:
        """draw the stored pedigree"""
        from .visualize import draw_pedigree_from_records

        records = [
            [rec.name, rec.parent1, rec.parent2]
            for rec in self.pedigree.records
        ]
        return draw_pedigree_from_records(records, **kwargs)

    def simulate(self) -> SimulationResult:
        """run the simulation with the stored seed"""
        return simulate_pedigree(
            pedigree=self.pedigree,
            genome=self.genome,
            seed=self.seed,
        )
