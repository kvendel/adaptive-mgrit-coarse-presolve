# Adaptive MGRIT Coarse Presolve

This repository contains the source code, Komondor SLURM scripts, numerical data, and plotting files accompanying the manuscript **"A Coarse Presolve Strategy for Adaptive Time Grid Construction in Multigrid Reduction in Time"**.

The committed numerical data can be inspected and used to regenerate the manuscript figures without repeating the full HPC experiments. The C programs and SLURM scripts used to generate the new numerical results are also included.

## Repository contents

- `xbraid/`: XBraid, included as a Git submodule and pinned to a specific commit.
- `experiments/fig6/` through `experiments/fig13/`: material associated with the corresponding manuscript figures.
- `experiments/fig14-16/`: common source code, scripts, and timing data for Figures 14-16.
- `experiments/plotting/`: `figures.ipynb`, `plot_functions.py`, and the generated figure directory.
- `Makefile`: top-level build interface for the XBraid library and the experiment programs.

Each experiment directory contains the applicable subset of:

- `src/`: C source files; generated executables are placed next to the sources;
- `scripts/`: Komondor SLURM batch scripts;
- `results/` or figure-specific result directories: numerical data used by the plotting code;
- `Makefile`: compilation rules for that experiment.

## Scope of the reproducibility material

This repository contains the source codes, execution scripts, numerical results, and plotting routines for the new computational experiments reported in Figures 6-16.

Figures 2-4 revisit numerical results and refinement strategies from the authors' earlier work. The corresponding reproducibility material is available at:

https://github.com/kvendel/Computational-stability-in-adaptive-parallel-in-time-Runge-Kutta-methods

Figures 1 and 5 are schematic or analytical illustrations and do not require execution of the C experiment programs.

## Clone the repository

Clone the repository together with XBraid:

```bash
git clone --recurse-submodules https://github.com/kvendel/adaptive-mgrit-coarse-presolve.git
cd adaptive-mgrit-coarse-presolve
```

If the repository was cloned without submodules, initialize XBraid with:

```bash
git submodule update --init --recursive
```

Display the pinned XBraid commit with:

```bash
git -C xbraid rev-parse HEAD
```

## Software requirements

Building and running the full experiments requires:

- a C compiler;
- an MPI implementation;
- GNU Make;
- XBraid, supplied by the `xbraid/` submodule;
- SLURM for the supplied Komondor batch scripts.

Regenerating the figures from the committed data requires:

- Python 3;
- Jupyter Notebook;
- NumPy;
- Matplotlib.

The Python dependencies can be installed with:

```bash
python -m pip install -r experiments/plotting/requirements.txt
```

## Build on Komondor

The reported experiments were compiled on Komondor using the default Cray programming environment. The Cray compiler wrapper `cc` was used for C and MPI C compilation.

Build the XBraid library and all experiment programs from the repository root:

```bash
make
```

The top-level Makefile passes these defaults recursively:

```text
CC=cc
MPICC=cc
MPICXX=CC
```

The upstream `xbraid/makefile.inc` is not modified. Only the XBraid library `xbraid/braid/libbraid.a` is built; XBraid example programs, including examples with optional HYPRE dependencies, are not built.

Useful targets include:

```bash
make xbraid              # build xbraid/braid/libbraid.a only
make fig6                # build programs associated with Figure 6
make fig7                # build programs associated with Figure 7
make fig14-16            # build programs associated with Figures 14-16
make list                # list detected parallel and sequential sources
make clean-experiments   # remove experiment executables
make clean-xbraid        # remove the XBraid library build products
make clean               # remove all build products
```

On another system, override the compiler variables if needed, for example:

```bash
make CC=gcc MPICC=mpicc MPICXX=mpicxx
```

## Run experiments on Komondor

The C programs write to paths relative to the current working directory, such as `results/...`. Therefore, each SLURM job must use the corresponding experiment directory as its working directory. The batch scripts invoke executables through paths of the form `./src/<executable>`.

### Submit from an experiment directory

For example, to submit the Figure 7 experiment:

```bash
cd experiments/fig7
sbatch scripts/run_all.sh
```

### Submit from the repository root

Use `sbatch --chdir` to set the correct working directory:

```bash
sbatch \
  --chdir="$PWD/experiments/fig7" \
  experiments/fig7/scripts/run_all.sh
```

The same pattern applies to every experiment directory. Exact commands and script names are documented in the local `README.md` files.

The scripts are Komondor-specific. On another HPC system, the SLURM account, partition, resources, modules, and execution options may need to be adapted.

## Numerical data

The committed result files are the data used for the submitted manuscript. They allow the figures to be regenerated without rerunning the HPC experiments.

Wall-clock measurements are system dependent. New timings can differ because of hardware, MPI configuration, process placement, network behavior, software versions, and system load. Figures 14-16 use medians of the available repeated timing measurements as implemented in `experiments/plotting/plot_functions.py`.

## Regenerate the figures

Start Jupyter from the plotting directory:

```bash
cd experiments/plotting
jupyter notebook figures.ipynb
```

Run the notebook cells in order. Generated PDF files are written to:

```text
experiments/plotting/figures/
```

Access to Komondor is not required when regenerating the figures from the committed result files.

## Third-party software

XBraid is maintained separately and is included only as a Git submodule. XBraid remains subject to its own license; see the license information in `xbraid/`.
