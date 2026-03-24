![License](https://img.shields.io/badge/license-MIT-green.svg)

# recombigraph

**Pedigree-based recombination simulation with explicit homolog tracking**

`recombigraph` simulates recombination along user-defined pedigrees while explicitly tracking the ancestry of chromosome segments through time. It's designed for applications where inheritance structure matters, like:

- studying pedigree-based genetic processes
- benchmarking inference methods

---

## Features

- Flexible diploid pedigree specification
- Explicit tracking of:
  - parent homolog IDs
  - founder homolog IDs
- Multiple chromosomes
- Built-in pedigree visualization
- Future (top priority): ARG reconstruction
- Future: different map functions
- Future: tetraploidy

---

## Installation

Install from source for now:

```bash
git clone https://github.com/pmckenz1/recombigraph.git
cd recombigraph
pip install -e .
```
---

## Quick Start

```
import recombigraph as rg

gen_list = [
    ['P0', 'NA', 'NA'],
    ['P1', 'NA', 'NA'],
    ['P2', 'NA', 'NA'],
    ['P3', 'NA', 'NA'],
    ['P4', 'NA', 'NA'],
    ['F1_0', 'P3', 'P1'],
    ['F1_1', 'P0', 'P0'],
    ['F1_2', 'P0', 'P3'],
    ['F1_3', 'P4', 'P3'],
    ['F1_4', 'P1', 'P2'],
    ['F2_0', 'F1_4', 'F1_3'],
    ['F2_1', 'F1_1', 'F1_3'],
    ['F2_2', 'F1_2', 'F1_4'],
    ['F2_3', 'F1_1', 'F1_1'],
    ['F2_4', 'F1_3', 'F1_4'],
    ['F3_0', 'F2_1', 'F2_2'],
    ['F3_1', 'F2_1', 'F2_2'],
    ['F3_2', 'F2_0', 'F2_4'],
    ['F3_3', 'F2_0', 'F2_0'],
    ['F3_4', 'F2_2', 'F2_2'],
]

model = rg.PedigreeModel(
    pedigree=gen_list,
    chromosomes={"A": 100.0, "B": 50.0},
    seed=123,
)

# visualize pedigree
model.draw_pedigree()
```

![Pedigree](./content/pedigree.svg)

```
# run simulation
result = model.simulate()

# inspect an individual
result.individuals["F3_0"]
```

```
# example output: result.individuals["F3_0"]

SimIndividual(
	individual_id='F3_0', 
	time=3, 
	homologs_by_chromosome={
		'A': [
			Homolog(
				homolog_id=60, 
				chromosome='A', 
				individual_id='F3_0', 
				time=3, 
				length=100.0, 
				segments=[
				Segment(left=0.0, right=36.03536993918074, parent_homolog_id=45, founder_homolog_id=16), 
				Segment(left=36.03536993918074, right=42.729877128295016, parent_homolog_id=45, founder_homolog_id=13), 
				Segment(left=42.729877128295016, right=64.05099368325804, parent_homolog_id=45, founder_homolog_id=12), 
				Segment(left=64.05099368325804, right=100.0, parent_homolog_id=45, founder_homolog_id=16)
				]
				), 
			Homolog(
				homolog_id=61, 
				chromosome='A', 
				individual_id='F3_0', 
				time=3, 
				length=100.0, 
				segments=[
				Segment(left=0.0, right=58.323841015105806, parent_homolog_id=48, founder_homolog_id=1), 
				Segment(left=58.323841015105806, right=61.2641671818151, parent_homolog_id=48, founder_homolog_id=0), 
				Segment(left=61.2641671818151, right=100.0, parent_homolog_id=48, founder_homolog_id=13)
				]
				)
			], 
		'B': [
			Homolog(
				homolog_id=62, 
				chromosome='B', 
				individual_id='F3_0', 
				time=3, 
				length=50.0, 
				segments=[
				Segment(left=0.0, right=40.55110556590415, parent_homolog_id=47, founder_homolog_id=19), 
				Segment(left=40.55110556590415, right=50.0, parent_homolog_id=47, founder_homolog_id=18)
				]
				), 
			Homolog(
				homolog_id=63, 
				chromosome='B', 
				individual_id='F3_0', 
				time=3, 
				length=50.0, 
				segments=[
				Segment(left=0.0, right=33.97120742047587, parent_homolog_id=50, founder_homolog_id=2), 
				Segment(left=33.97120742047587, right=50.0, parent_homolog_id=50, founder_homolog_id=14)
				]
				)
			]
		}
	)
```