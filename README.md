# PTOPTW (Athens instance demo)

C++ implementation of an ILS-based solver for the **Team Orienteering Problem with Time Windows (TOPTW)**, including an experimental **partitioned local search** variant aimed at reducing runtime on larger instances. This repo is set up as a runnable demo using a custom **AthensTopology** instance plus a small Leaflet visualization.

## Thesis context

The Team Orienteering Problem with Time Windows (TOPTW) is an NP-hard extension of the Orienteering Problem. Since exact methods do not scale to large instances, heuristic and approximate algorithms are typically used. This work builds on **Iterated Local Search (ILS)** and explores reducing execution time by **partitioning the problem graph**, applying local search on sub-graphs, and addressing the issues introduced by partitioning—trading some solution quality for faster execution on large inputs.

---

## Requirements

- CMake (>= 3.16 recommended)
- C++20 compiler (GCC/Clang)
- (Optional, for visualization) Python 3 (for a local static server) or any HTTP server

---

## Build

From the repository root:

```bash
cmake -S . -B build
cmake --build build -j
```

Binary:

```bash
./build/ptoptw
```

---

## Run (AthensTopology)

**Important:** run from the repository root (paths are relative).

The solver currently loads instances from:

`./instances/<folder>/<instance>.txt`

This repo includes the Athens instance at:

`./instances/Athens/AthensTopology.txt`

Run it like this:

```bash
./build/ptoptw -f Athens -i AthensTopology -m 4 -s 4 -c
```

Options:

- `-f` : dataset folder name under `./instances/` (here: `Athens`)
- `-i` : instance name (without `.txt`) (here: `AthensTopology`)
- `-m` : number of routes / teams
- `-s` : number of partitions / sub-intervals used by the partitioned local search variant
- `-c` : use custom travel times embedded in the instance file (required when the instance includes `D` duration lines)

Optional flags:

- `-t <seconds>` : time-limited execution (e.g. `-t 10`)
- `-w` : write solution output (if enabled in the solver)
- `-h` : show help

Example (short run):

```bash
./build/ptoptw -f Athens -i AthensTopology -m 4 -s 1 -c -t 5
```

---

## Visualization (Leaflet)

Static visualization lives under:

`viz/topology/static/`

Serve it locally:

```bash
cd viz/topology/static
python3 -m http.server 8000
```

Then open:

`http://localhost:8000`

If you need to regenerate or preprocess visualization assets, check:

- `viz/topology/create.js`
- `viz/topology/package.json`

---

## Screenshot

![Athens instance demo](docs/athens_instance.png)
*Figure 1 — Unpartitioned (s=1) Athens instance with 7 routes (m=7).*

---

## Notes / limitations

- This repo is currently documented for the **AthensTopology** instance format. Other benchmark formats (e.g., Cordeau/Solomon TOPTW datasets) are not supported by the current parser.
- The instance included in the repository may not match the “complete” dataset used during thesis experiments; treat it as a reproducible demo instance for running and visualization.

---

