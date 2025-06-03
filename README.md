# log-analysis-tools

## Overview

This repository scripts operators on `.mcap` files and generates three types of plots:

1. Position and velocity plots of the global EKF along with its inputs overlayed
2. Comparison plots between two `.mcap` files 
3. In detail GSOF plots

## Prerequisite

Install missing packages:
```
pip install mcap
pip install pymap3d
```

## Usage

```
Usage: generate_mcap_plots|generate_gsof_plots.py [OPTIONS]

Options:
  -l, --log-path TEXT  Path to .mcap
  -d, --dir-path TEXT  Path to a root directory containing .mcap files
  -m, --multiprocess   Enable multiprocessing
  --help               Show this message and exit.
```

```
Usage: generate_mcap_comparison_plots.py [OPTIONS]

Options:
  -d, --dir-path TEXT  Path to a root directory containing .mcap files
  -m, --multiprocess   Enable multiprocessing
  --help               Show this message and exit.
``