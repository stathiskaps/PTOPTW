## Thesis context

This project studies the Team Orienteering Problem with Time Windows (TOPTW), an NP-hard extension of the Orienteering Problem. It implements an Iterated Local Search (ILS) approach and explores reducing runtime by partitioning the problem graph and applying local search per sub-graph, trading off solution score for faster execution on large instances.

# PTOPTW (Athens instance demo)

This repository contains a C++ implementation of an Iterated Local Search approach for an orienteering / routing problem with time windows, demonstrated on a custom **AthensTopology** instance.

The project includes:
- a solver (`ptoptw`) that reads a single instance file
- a small Leaflet-based visualization under `viz/topology/` to inspect the instance / routes

> Note: the included AthensTopology instance is **incomplete** (not all POIs / routes are present). This repo is mainly a reproducible demo of the solver + visualization pipeline.

---

## Requirements

- CMake (>= 3.16 recommended)
- A C++20 compiler (GCC/Clang)
- (Optional, for visualization) Node.js + npm

---

## Build

From the repository root:

```bash
cmake -S . -B build
cmake --build build -j
```

This produces:

```bash
./build/ptoptw
```

---

## Quickstart: run the Athens instance

### 1) Run

**Important:** run from the repo root (paths are relative).

```bash
./build/ptoptw -f Athens -i athens -m 4 -s 4 -c
```

- `-f` : dataset folder name under `./instances/` (here: `Athens`)
- `-i` : instance name (without `.txt`) (here: `athens`)
- `-m` : numeric parameter (internally used as number of routes / walks)
- `-s` : numeric parameter (internally used as number of sub-intervals / partitions)
- `-c` : use custom travel times embedded in the instance file (required for AthensTopology which includes `D` duration lines)

Optional flags:
- `-t <seconds>` : time-limited execution (e.g. `-t 10`)
- `-w` : enable writing solution output (see solver output logs / generated files)
- `-a` : run all instances under `./instances/Solomon/` and `./instances/Cordeau/` (legacy; not recommended for this Athens-only demo)

Help:

```bash
./build/ptoptw -h
```

---

## Visualization (Leaflet)

The visualization lives under:

`viz/topology/static/`

If you only want to open the static page:

```bash
cd viz/topology/static
python3 -m http.server 8000
```

Then open:

`http://localhost:8000`

If you want to regenerate or preprocess data (if applicable), check:

- `viz/topology/create.js`
- `viz/topology/package.json`

---

## Screenshot

![Unpartitioned (S=1) Athens instance with 7 routes (m = 7)](docs/athens_instance.png)

---

## Notes / limitations

- The solver currently constructs the instance file path internally as:
  `./instances/<folder>/<instance>.txt`
  (so `-f` and `-i` are identifiers, not full paths).
- The included Athens topology is incomplete; treat it as a demo dataset.

---
