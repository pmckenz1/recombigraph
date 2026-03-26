import pytest

import pedigraph_sim as pg


def make_three_generation_model():
    return pg.PedigreeModel(
        pedigree=[
            ("gp1", None, None),
            ("gp2", None, None),
            ("gp3", None, None),
            ("gp4", None, None),
            ("p1", "gp1", "gp2"),
            ("p2", "gp3", "gp4"),
            ("child", "p1", "p2"),
        ],
        chromosomes={"chr1": 100.0, "chr2": 50.0},
        seed=2,
    )


def test_local_forests_end_to_end_export():
    model = make_three_generation_model()
    result = model.simulate()

    sample_ids = [
        homolog.homolog_id
        for homolog in result.individuals["child"].homologs_by_chromosome["chr1"]
    ]

    seq = result.local_forests("chr1", sample_ids)
    records = pg.to_newick_records(result, seq)

    assert len(seq) >= 1
    assert seq.breakpoints()[0] == 0.0
    assert seq.breakpoints()[-1] == 100.0
    assert len(records) == len(seq)

    for forest, record in zip(seq, records):
        assert forest.chromosome == "chr1"
        assert forest.sample_homolog_ids == tuple(sample_ids)
        assert forest.left < forest.right
        assert forest.nodes()
        assert forest.roots()
        assert record["chromosome"] == "chr1"
        assert record["left"] == forest.left
        assert record["right"] == forest.right
        assert record["n_edges"] == len(forest.edges)
        assert len(record["newicks"]) == len(forest.roots())


def test_pedigree_rejects_exactly_one_parent():
    with pytest.raises(ValueError, match="either 0 parents or 2 parents"):
        pg.PedigreeModel(
            pedigree=[
                ("p1", None, None),
                ("child", "p1", None),
            ],
            chromosomes={"chr1": 10.0},
        )


def test_local_forests_rejects_mixed_chromosome_samples():
    model = make_three_generation_model()
    result = model.simulate()

    chr1_hid = result.individuals["child"].homologs_by_chromosome["chr1"][0].homolog_id
    chr2_hid = result.individuals["child"].homologs_by_chromosome["chr2"][0].homolog_id

    with pytest.raises(ValueError, match="belongs to chromosome"):
        result.local_forests("chr1", [chr1_hid, chr2_hid])


def test_local_forests_rejects_unknown_chromosome():
    model = pg.PedigreeModel(
        pedigree=[
            ("p1", None, None),
            ("p2", None, None),
            ("child", "p1", "p2"),
        ],
        chromosomes={"chr1": 10.0},
        seed=1,
    )
    result = model.simulate()
    sample_ids = [
        homolog.homolog_id
        for homolog in result.individuals["child"].homologs_by_chromosome["chr1"]
    ]

    with pytest.raises(ValueError, match="Unknown chromosome"):
        result.local_forests("missing", sample_ids)
