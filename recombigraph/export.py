from __future__ import annotations

from typing import Any, Optional
from typing import TYPE_CHECKING

from recombigraph.arg import LocalForest, LocalForestSequence

if TYPE_CHECKING:
    from .simulate import SimulationResult

def _build_homolog_lookup(result: "SimulationResult") -> dict[int, Any]:
    """build a homolog lookup keyed by homolog id"""
    out: dict[int, Any] = {}
    for individual in result.individuals.values():
        for homologs in individual.homologs_by_chromosome.values():
            for homolog in homologs:
                out[homolog.homolog_id] = homolog
    return out


def _sample_label(result: "SimulationResult", homolog_id: int) -> str:
    """build a stable label for a sampled homolog"""
    h = _build_homolog_lookup(result)[homolog_id]
    chrom = h.chromosome
    ind = h.individual_id

    # try to infer h0 / h1 from the individual's homolog ordering
    homologs = result.individuals[ind].homologs_by_chromosome[chrom]
    for idx, x in enumerate(homologs):
        if x.homolog_id == homolog_id:
            return f"{ind}_{chrom}_h{idx}"
    return f"{ind}_{chrom}_{homolog_id}"


def _branch_length(result, parent_id: int, child_id: int) -> float:
    """return the time difference between parent and child homologs"""
    lookup = _build_homolog_lookup(result)
    parent = lookup[parent_id]
    child = lookup[child_id]
    return float(child.time - parent.time)

def forest_to_newicks(
    result: "SimulationResult",
    forest: LocalForest,
    *,
    label_samples: bool = True,
    include_internal_labels: bool = False,
) -> list[str]:
    """serialize one local forest into one or more newick strings"""
    lookup = _build_homolog_lookup(result)
    children_map = forest.children_map()
    sample_set = set(forest.sample_homolog_ids)

    def node_label(node_id: int) -> str:
        if node_id in sample_set and label_samples:
            return _sample_label(result, node_id)
        if include_internal_labels:
            return str(node_id)
        return ""

    def rec(node_id: int, parent_id: Optional[int] = None) -> str:
        children = children_map.get(node_id, ())
        label = node_label(node_id)

        if not children:
            if parent_id is None:
                return f"{label or node_id};"
            blen = _branch_length(result, parent_id, node_id)
            return f"{label or node_id}:{blen}"

        sub = ",".join(rec(child, node_id) for child in children)
        if parent_id is None:
            suffix = label if label else ""
            return f"({sub}){suffix};"

        blen = _branch_length(result, parent_id, node_id)
        suffix = label if label else ""
        return f"({sub}){suffix}:{blen}"

    return [rec(root) for root in forest.roots()]

def to_newick_records(
    result: "SimulationResult",
    seq: LocalForestSequence,
    *,
    as_list: bool = True,
) -> list[dict[str, Any]]:
    """convert a local forest sequence into record dictionaries"""
    out: list[dict[str, Any]] = []
    for forest in seq:
        newicks = forest_to_newicks(result, forest)
        out.append(
            {
                "chromosome": forest.chromosome,
                "left": forest.left,
                "right": forest.right,
                "span": forest.span(),
                "n_roots": len(forest.roots()),
                "n_nodes": len(forest.nodes()),
                "n_edges": len(forest.edges),
                "newicks": newicks if as_list else "|".join(newicks),
            }
        )
    return out

def to_dataframe(result: "SimulationResult", seq: LocalForestSequence) -> Any:
    """convert newick records into a pandas dataframe"""
    import pandas as pd

    return pd.DataFrame(to_newick_records(result, seq, as_list=True))

def to_tskit(result: "SimulationResult", seq: LocalForestSequence) -> Any:
    """convert a local forest sequence into a tskit tree sequence"""
    try:
        import tskit
    except ImportError as e:
        raise ImportError("to_tskit() requires tskit to be installed.") from e

    lookup = _build_homolog_lookup(result)
    sample_set = set(seq.sample_homolog_ids())

    tables = tskit.TableCollection(sequence_length=seq.right)

    # store a little metadata so downstream tools can recover ids
    tables.individuals.metadata_schema = tskit.MetadataSchema.permissive_json()
    tables.nodes.metadata_schema = tskit.MetadataSchema.permissive_json()

    node_row_for_hid: dict[int, int] = {}
    individual_row_for_individual_id: dict[str, int] = {}

    def get_individual_row(individual_id: str) -> int:
        """return the tskit row id for an individual"""
        if individual_id not in individual_row_for_individual_id:
            individual_row_for_individual_id[individual_id] = tables.individuals.add_row(
                metadata={"individual_id": individual_id}
            )
        return individual_row_for_individual_id[individual_id]

    for hid in seq.nodes():
        h = lookup[hid]
        ind_row = get_individual_row(h.individual_id)
        flags = tskit.NODE_IS_SAMPLE if hid in sample_set else 0
        node_row = tables.nodes.add_row(
            flags=flags,
            time=float(h.time),
            individual=ind_row,
            metadata={
                "homolog_id": hid,
                "chromosome": h.chromosome,
                "individual_id": h.individual_id,
            },
        )
        node_row_for_hid[hid] = node_row

    for forest in seq:
        for parent_hid, child_hid in forest.edges:
            tables.edges.add_row(
                left=forest.left,
                right=forest.right,
                parent=node_row_for_hid[parent_hid],
                child=node_row_for_hid[child_hid],
            )

    tables.sort()
    return tables.tree_sequence()
