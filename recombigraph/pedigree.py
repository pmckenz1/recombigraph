from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class PedigreeRecord:
    name: str
    parent1: Optional[str]
    parent2: Optional[str]

    @property
    def is_founder(self) -> bool:
        return self.parent1 is None and self.parent2 is None


class Pedigree:
    def __init__(self, records):
        self.records = self._normalize_records(records)
        self._record_dict = {rec.name: rec for rec in self.records}

        self._validate_unique_names()
        self._validate_parent_references()

        self.records = self._topological_sort()
        self._record_dict = {rec.name: rec for rec in self.records}

        self.generation_map = self._compute_generation_map()

    def _normalize_records(self, records) -> list[PedigreeRecord]:
        normed = []
        for rec in records:
            if isinstance(rec, PedigreeRecord):
                normed.append(rec)
            else:
                name, parent1, parent2 = rec
                parent1 = None if parent1 in ("NA", None) else parent1
                parent2 = None if parent2 in ("NA", None) else parent2
                normed.append(PedigreeRecord(name, parent1, parent2))
        return normed

    def _validate_unique_names(self):
        names = [rec.name for rec in self.records]
        if len(names) != len(set(names)):
            raise ValueError("Pedigree contains duplicate individual names")

    def _validate_parent_references(self):
        names = {rec.name for rec in self.records}
        for rec in self.records:
            for parent in (rec.parent1, rec.parent2):
                if parent is not None and parent not in names:
                    raise ValueError(f"Parent {parent!r} of {rec.name!r} not found in pedigree")

    def _topological_sort(self) -> list[PedigreeRecord]:
        remaining = {rec.name: rec for rec in self.records}
        resolved = set()
        ordered = []

        while remaining:
            progressed = False
            for name in list(remaining):
                rec = remaining[name]
                parents_ok = all(
                    parent is None or parent in resolved
                    for parent in (rec.parent1, rec.parent2)
                )
                if parents_ok:
                    ordered.append(rec)
                    resolved.add(name)
                    del remaining[name]
                    progressed = True

            if not progressed:
                unresolved = ", ".join(sorted(remaining))
                raise ValueError(
                    "Pedigree could not be topologically sorted. "
                    f"Possible cycle or missing dependency among: {unresolved}"
                )

        return ordered

    def _compute_generation_map(self) -> dict[str, int]:
        generation = {}
        for rec in self.records:
            if rec.is_founder:
                generation[rec.name] = 0
            else:
                generation[rec.name] = max(
                    generation[rec.parent1] if rec.parent1 is not None else 0,
                    generation[rec.parent2] if rec.parent2 is not None else 0,
                ) + 1
        return generation

    def __iter__(self):
        return iter(self.records)

    def __len__(self):
        return len(self.records)

    def __getitem__(self, name: str) -> PedigreeRecord:
        return self._record_dict[name]

    @property
    def founders(self) -> list[PedigreeRecord]:
        return [rec for rec in self.records if rec.is_founder]