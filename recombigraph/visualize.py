from matplotlib import pyplot as plt
import networkx as nx


def _is_missing_parent(x: str | None) -> bool:
    """return whether a parent field should be treated as missing"""
    return x is None or x == "NA"


def compute_generation_map(
    records: list[tuple[str, str | None, str | None]],
) -> dict[str, int]:
    """compute generation indices from pedigree style records"""
    generation: dict[str, int] = {}
    for name, parent1, parent2 in records:
        if _is_missing_parent(parent1) and _is_missing_parent(parent2):
            generation[name] = 0
        else:
            generation[name] = max(generation[parent1], generation[parent2]) + 1
    return generation


def build_pedigree_graph(
    records: list[tuple[str, str | None, str | None]],
    arrows_to_parents: bool = False,
):
    """build a directed pedigree graph for plotting"""
    G = nx.DiGraph()
    generation = compute_generation_map(records)

    for name, parent1, parent2 in records:
        G.add_node(name, time=generation[name])

        if not _is_missing_parent(parent1):
            if arrows_to_parents:
                G.add_edge(name, parent1)
            else:
                G.add_edge(parent1, name)

        if not _is_missing_parent(parent2):
            if arrows_to_parents:
                G.add_edge(name, parent2)
            else:
                G.add_edge(parent2, name)

    return G

def draw_pedigree_from_records(
    records: list[tuple[str, str | None, str | None]],
    figsize: tuple[int, int] = (10, 6),
    arrows_to_parents: bool = False,
    node_alpha: float = 0.6,
    node_size: int = 1200,
    font_size: int = 8,
    cmap: str = "viridis",
    savepath: str | None = None,
    dpi: int = 300,
):
    """draw a simple pedigree graph from pedigree style records"""
    G = build_pedigree_graph(records, arrows_to_parents=arrows_to_parents)
    pos = nx.multipartite_layout(G, subset_key="time", align="vertical")

    plt.figure(figsize=figsize)
    
    # nodes
    nx.draw_networkx_nodes(
        G,
        pos,
        node_color=[G.nodes[n]["time"] for n in G.nodes],
        cmap=cmap,
        node_size=node_size,
        alpha=node_alpha,
    )

    # edges
    nx.draw_networkx_edges(
        G,
        pos,
        arrows=True,
        alpha=0.7,
    )

    # labels (no alpha applied)
    nx.draw_networkx_labels(
        G,
        pos,
        font_size=font_size,
        font_color="black",
    )
    plt.gca().invert_yaxis()
    plt.axis("off")

    if savepath is not None:
        plt.savefig(savepath, dpi=dpi, bbox_inches="tight")

    plt.show()
