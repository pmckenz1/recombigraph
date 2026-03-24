from dataclasses import dataclass


@dataclass(frozen=True)
class ChromosomeSpec:
    name: str
    length: float


class GenomeSpec:
    def __init__(self, chromosomes):
        self.chromosomes = self._normalize_chromosomes(chromosomes)
        self._chrom_dict = {chrom.name: chrom for chrom in self.chromosomes}

        self._validate_unique_names()
        self._validate_lengths()

    def _normalize_chromosomes(self, chromosomes) -> list[ChromosomeSpec]:
        normed = []

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

    def _validate_unique_names(self):
        names = [chrom.name for chrom in self.chromosomes]
        if len(names) != len(set(names)):
            raise ValueError("Genome contains duplicate chromosome names")

    def _validate_lengths(self):
        for chrom in self.chromosomes:
            if chrom.length <= 0:
                raise ValueError(
                    f"Chromosome {chrom.name!r} has non-positive length {chrom.length}"
                )

    def __iter__(self):
        return iter(self.chromosomes)

    def __len__(self):
        return len(self.chromosomes)

    def __getitem__(self, name: str) -> ChromosomeSpec:
        return self._chrom_dict[name]

    @property
    def chromosome_names(self) -> list[str]:
        return [chrom.name for chrom in self.chromosomes]