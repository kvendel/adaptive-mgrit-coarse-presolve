# Figure 12

This directory contains the reproducibility material for local presolve step-ratio diagrams for the test problems.

## Contents

- `src/`: C source files and locally generated executables.
- `Makefile`: compilation rules for the programs in `src/`.
- `results/`: numerical data used by the plotting code.

## Build

From the repository root:

```bash
make fig12
```

Executables are created in `experiments/fig12/src/`.

## Experiment execution

No dedicated SLURM batch script is included in this directory. The committed data are sufficient to regenerate the figure. The plotting notebook reads the required step-sequence data directly from the repository.

## Regenerate the figure

From the repository root:

```bash
cd experiments/plotting
jupyter notebook figures.ipynb
```

Run the notebook cells in order. The corresponding plotting function is defined in `plot_functions.py`.

## Platform note

The supplied build and batch configuration was tested on Komondor. Other systems may require different compiler, MPI, SLURM, module, or resource settings.

