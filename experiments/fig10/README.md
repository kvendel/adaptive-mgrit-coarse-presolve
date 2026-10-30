# Figure 10

This directory contains the reproducibility material for terminal-error response as a function of the requested tolerance for the four controllers.

## Contents

- `src/`: C source files and locally generated executables.
- `Makefile`: compilation rules for the programs in `src/`.
- `scripts/`: Komondor SLURM batch scripts.
- `results/`: numerical data used by the plotting code.

## Build

From the repository root:

```bash
make fig10
```

Executables are created in `experiments/fig10/src/`.

## Run on Komondor

The job working directory must be this experiment directory because the executables use relative output paths.

### Submit from this experiment directory

```bash
cd experiments/fig10
sbatch scripts/prop.sh
```

### Submit from the repository root

```bash
sbatch \
  --chdir="$PWD/experiments/fig10" \
  experiments/fig10/scripts/prop.sh
```

The scripts launch executables from `./src/`. Output files are therefore written to the result location relative to this experiment directory.

## Regenerate the figure

From the repository root:

```bash
cd experiments/plotting
jupyter notebook figures.ipynb
```

Run the notebook cells in order. The corresponding plotting function is defined in `plot_functions.py`.

## Platform note

The supplied build and batch configuration was tested on Komondor. Other systems may require different compiler, MPI, SLURM, module, or resource settings.

