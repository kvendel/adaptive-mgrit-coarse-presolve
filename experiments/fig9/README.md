# Figure 9

This directory contains the reproducibility material for Van der Pol step-size sequences produced by the I, exponential, PI, and PID controllers.

## Contents

- `src/`: C source files and locally generated executables.
- `Makefile`: compilation rules for the programs in `src/`.
- `scripts/`: Komondor SLURM batch scripts.
- `results/`: numerical data used by the plotting code.

## Build

From the repository root:

```bash
make fig9
```

Executables are created in `experiments/fig9/src/`.

## Run on Komondor

The job working directory must be this experiment directory because the executables use relative output paths.

### Submit from this experiment directory

```bash
cd experiments/fig9
sbatch scripts/run_all.sh
```

### Submit from the repository root

```bash
sbatch \
  --chdir="$PWD/experiments/fig9" \
  experiments/fig9/scripts/run_all.sh
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

