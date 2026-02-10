# SF6-Circuit-Breakers-Interruption-Data
Dataset and scripts for experiments in a 252 kV self-blast SF6 circuit breaker
## Overview
Electrode erosion caused by switching arcs is a key factor limiting electrical lifetime.
This dataset supports a coupled MHD–turbulence–phase-change study and experimental validation on a **252 kV self-blast SF6 circuit breaker**.

## Contents
- `data.csv/raw/`: raw measurements and exported data (Travel)
- `data/processed/`: pictures and records during the experiment process (Pictures)
- `data/data_dictionary.scm`: variable names, units, descriptions (Materials)
- `scripts/`: preprocessing and analysis scripts (Fitting)
- `simulation/inputs/`: simulation input files (UDF)

## How to use
1. See  `scripts/` for variable definitions and units.
2. Run preprocessing scripts in `scripts/preprocess/` to regenerate `data/processed/`.
3. Run analysis scripts in `scripts/analysis/` to reproduce key results/figures.

## File formats
- Tabular data (Excel)
- UDF (C)

## License
- Data: [e.g., CC BY 4.0]
- Code: [e.g., MIT]

## Citation
If you use this dataset, please cite the manuscript and the archived dataset DOI.

## Contact
Xinyuan Hong, Chongqing University
Email: x.hong@stu.cqu.edu.cn
