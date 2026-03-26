import pedigraph_sim as pg


def test_basic_simulation_builds_expected_individuals():
    model = pg.PedigreeModel(
        pedigree=[
            ("child", "parent1", "parent2"),
            ("parent2", None, None),
            ("parent1", None, None),
        ],
        chromosomes={"chr1": 100.0},
        seed=1,
    )

    result = model.simulate()

    assert set(result.individuals) == {"parent1", "parent2", "child"}
    assert result.pedigree.generation_map == {
        "parent2": 0,
        "parent1": 0,
        "child": 1,
    }
    assert len(result.individuals["child"].homologs_by_chromosome["chr1"]) == 2
