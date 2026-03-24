from matplotlib import pyplot as plt
import networkx as nx


def compute_generation_map(records):
    generation = {}
    for name, parent1, parent2 in records:
        if parent1 == "NA" and parent2 == "NA":
            generation[name] = 0
        else:
            generation[name] = max(generation[parent1], generation[parent2]) + 1
    return generation


def build_pedigree_graph(records, arrows_to_parents=False):
    G = nx.DiGraph()
    generation = compute_generation_map(records)

    for name, parent1, parent2 in records:
        G.add_node(name, time=generation[name])

        if parent1 != "NA":
            if arrows_to_parents:
                G.add_edge(name, parent1)
            else:
                G.add_edge(parent1, name)

        if parent2 != "NA":
            if arrows_to_parents:
                G.add_edge(name, parent2)
            else:
                G.add_edge(parent2, name)

    return G


def draw_pedigree_from_records(
    records,
    figsize=(10, 6),
    arrows_to_parents=False,
    node_alpha=0.6,
    node_size=1200,
    font_size=8,
    cmap="viridis",
):
    G = build_pedigree_graph(records, arrows_to_parents=arrows_to_parents)
    pos = nx.multipartite_layout(G, subset_key="time", align="vertical")

    plt.figure(figsize=figsize)
    nx.draw_networkx(
        G,
        pos,
        with_labels=True,
        node_color=[G.nodes[n]["time"] for n in G.nodes],
        cmap=cmap,
        node_size=node_size,
        font_size=font_size,
        font_color="black",
        arrows=True,
        alpha=node_alpha,
    )
    plt.gca().invert_yaxis()
    plt.axis("off")
    plt.show()