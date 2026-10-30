# Figure 8

This directory contains the reproducibility material for terminal-error response as a function of the requested tolerance for the refinement-factor comparison.

## Contents

- `src/`: C source files and locally generated executables.
- `Makefile`: compilation rules for the programs in `src/`.
- `scripts/`: Komondor SLURM batch scripts.
- `exact_sols/`: reference solutions used to compute the reported errors.
- `results/`: numerical data used by the plotting code.

## Build

From the repository root:

```bash
make fig8
```

Executables are created in `experiments/fig8/src/`.

## Run on Komondor

The job working directory must be this experiment directory because the executables use relative output paths.

### Submit from this experiment directory

```bash
cd experiments/fig8
sbatch scripts/lin_BS.sh
```

```bash
cd experiments/fig8
sbatch scripts/lin_IE.sh
```

```bash
cd experiments/fig8
sbatch scripts/LV_BS.sh
```

```bash
cd experiments/fig8
sbatch scripts/LV_IE.sh
```

```bash
cd experiments/fig8
sbatch scripts/vdp_BS.sh
```

```bash
cd experiments/fig8
sbatch scripts/vdp_IE.sh
```

### Submit from the repository root

```bash
sbatch \
  --chdir="$PWD/experiments/fig8" \
  experiments/fig8/scripts/lin_BS.sh
```

```bash
sbatch \
  --chdir="$PWD/experiments/fig8" \
  experiments/fig8/scripts/lin_IE.sh
```

```bash
sbatch \
  --chdir="$PWD/experiments/fig8" \
  experiments/fig8/scripts/LV_BS.sh
```

```bash
sbatch \
  --chdir="$PWD/experiments/fig8" \
  experiments/fig8/scripts/LV_IE.sh
```

```bash
sbatch \
  --chdir="$PWD/experiments/fig8" \
  experiments/fig8/scripts/vdp_BS.sh
```

```bash
sbatch \
  --chdir="$PWD/experiments/fig8" \
  experiments/fig8/scripts/vdp_IE.sh
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

