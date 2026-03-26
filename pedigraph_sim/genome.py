from dataclasses import dataclass
from collections.abc import Iterator
from typing import Dict, Iterable, Tuple, Union


@dataclass(frozen=True)
class ChromosomeSpec:
    """chromosome name and length"""
    name: str
    length: float


ChromosomeLike = Union[ChromosomeSpec, Tuple[str, float]]
ChromosomeCollection = Union[Dict[str, float], Iterable[ChromosomeLike]]


class GenomeSpec:
    """validated genome specification"""

    def __init__(self, chromosomes: ChromosomeCollection) -> None:
        """build a genome from chromosome records or a name to length mapping"""
        self.chromosomes = self._normalize_chromosomes(chromosomes)
        self._chrom_dict = {chrom.name: chrom for chrom in self.chromosomes}

        self._validate_unique_names()
        self._validate_lengths()

    def _normalize_chromosomes(
        self,
        chromosomes: ChromosomeCollection,
    ) -> list[ChromosomeSpec]:
        """normalize supported chromosome inputs into chromosome specs"""
        normed: list[ChromosomeSpec] = []

        if isinstance(chromosomes, dict):
            for name, length in chromosomes.items():
                normed.append(ChromosomeSpec(name=name, length=float(length)))
        else:
            for chrom in chromosomes:
                if isinstance(chrom, ChromosomeSpec):
                    normed.append(chrom)
                else:
                    name, length = chrom
                    normed.append(ChromosomeSpec(name=name, length=float(length)))

        return normed

    def _validate_unique_names(self) -> None:
        """ensure chromosome names are unique"""
        names = [chrom.name for chrom in self.chromosomes]
        if len(names) != len(set(names)):
            raise ValueError("Genome contains duplicate chromosome names")

    def _validate_lengths(self) -> None:
        """ensure chromosome lengths are positive"""
        for chrom in self.chromosomes:
            if chrom.length <= 0:
                raise ValueError(
                    f"Chromosome {chrom.name!r} has non-positive length {chrom.length}"
                )

    def __iter__(self) -> Iterator[ChromosomeSpec]:
        """iterate over chromosome specs"""
        return iter(self.chromosomes)

    def __len__(self) -> int:
        """return the number of chromosomes"""
        return len(self.chromosomes)

    def __getitem__(self, name: str) -> ChromosomeSpec:
        """look up a chromosome by name"""
        return self._chrom_dict[name]

    @property
    def chromosome_names(self) -> list[str]:
        """return chromosome names in order"""
        return [chrom.name for chrom in self.chromosomes]
