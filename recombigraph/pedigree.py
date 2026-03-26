from dataclasses import dataclass
from collections.abc import Iterator
from typing import Iterable, Optional


@dataclass(frozen=True)
class PedigreeRecord:
    """pedigree row with zero or two parents"""
    name: str
    parent1: Optional[str]
    parent2: Optional[str]

    @property
    def is_founder(self) -> bool:
        """return whether this individual is a founder"""
        return self.parent1 is None and self.parent2 is None


PedigreeRecordLike = PedigreeRecord | tuple[str, Optional[str], Optional[str]]


class Pedigree:
    """validated pedigree sorted so parents precede children"""

    def __init__(self, records: Iterable[PedigreeRecordLike]) -> None:
        """build a pedigree from records or tuples"""
        self.records = self._normalize_records(records)
        self._record_dict = {rec.name: rec for rec in self.records}

        self._validate_unique_names()
        self._validate_parent_references()

        self.records = self._topological_sort()
        self._record_dict = {rec.name: rec for rec in self.records}

        self.generation_map = self._compute_generation_map()

    def _normalize_records(
        self,
        records: Iterable[PedigreeRecordLike],
    ) -> list[PedigreeRecord]:
        """normalize supported record inputs into pedigree records"""
        normed: list[PedigreeRecord] = []
        for rec in records:
            if isinstance(rec, PedigreeRecord):
                normed.append(rec)
            else:
                name, parent1, parent2 = rec
                parent1 = None if parent1 in ("NA", None) else parent1
                parent2 = None if parent2 in ("NA", None) else parent2
                normed.append(PedigreeRecord(name, parent1, parent2))
        return normed

    def _validate_unique_names(self) -> None:
        """ensure individual names are unique"""
        names = [rec.name for rec in self.records]
        if len(names) != len(set(names)):
            raise ValueError("Pedigree contains duplicate individual names")

    def _validate_parent_references(self) -> None:
        """ensure parent references are complete and valid"""
        names = {rec.name for rec in self.records}
        for rec in self.records:
            if (rec.parent1 is None) != (rec.parent2 is None):
                raise ValueError(
                    f"Individual {rec.name!r} must have either 0 parents or 2 parents, "
                    "not exactly 1."
                )
            for parent in (rec.parent1, rec.parent2):
                if parent is not None and parent not in names:
                    raise ValueError(f"Parent {parent!r} of {rec.name!r} not found in pedigree")
    def _topological_sort(self) -> list[PedigreeRecord]:
        """sort records so all parents appear before children"""
        remaining = {rec.name: rec for rec in self.records}
        resolved = set()
        ordered: list[PedigreeRecord] = []

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
        """assign generation indices from founders outward"""
        generation: dict[str, int] = {}
        for rec in self.records:
            if rec.is_founder:
                generation[rec.name] = 0
            else:
                generation[rec.name] = max(
                    generation[rec.parent1] if rec.parent1 is not None else 0,
                    generation[rec.parent2] if rec.parent2 is not None else 0,
                ) + 1
        return generation

    def __iter__(self) -> Iterator[PedigreeRecord]:
        """iterate over pedigree records in topological order"""
        return iter(self.records)

    def __len__(self) -> int:
        """return the number of pedigree records"""
        return len(self.records)

    def __getitem__(self, name: str) -> PedigreeRecord:
        """look up a pedigree record by name"""
        return self._record_dict[name]

    @property
    def founders(self) -> list[PedigreeRecord]:
        """return founder records in stored order"""
        return [rec for rec in self.records if rec.is_founder]
