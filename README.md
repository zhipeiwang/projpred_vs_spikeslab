# Research Archive

This repository contains all relevant code and information to reproduce the results in the master's thesis:  
_Comparing Bayesian Post-estimation Variable Selection Methods: Projection Predictive Variable Selection Versus Stochastic Search Variable Selection_.

We compare the performance of two Bayesian variable selection methods Projection Predictive Variable Selection (PPVS) and Stochastic Search Variable Selection (SSVS) through two simulation studies: one in a low-dimensional setting and another in a high-dimensional setting. Two empirical datasets are also included to illustrate practical application. The full thesis is provided in this repository.

---
# Directory Structure
```
├── data
│   ├── p50
│   └── simdata
├── LICENSE
├── output
│   ├── figures
│   ├── output_study1
│   │   ├── ppvs
│   │   │   ├── cvvs4suggest_size_failed
│   │   │   ├── df_out
│   │   │   ├── prj_pred
│   │   │   ├── ppvs_out
│   │   │   ├── time_taken
│   │   │   └── summ_ref
│   │   └── ssvs
│   │   │   ├── df_out
│   │   │   ├── time_taken
│   │       └── ssvs_pred_results
│   ├── output_study2
│   │   ├── ppvs
│   │   │   ├── cvvs4suggest_size_failed
│   │   │   ├── df_out
│   │   │   ├── prj_pred
│   │   │   ├── ppvs_out
│   │   │   ├── time_taken
│   │   │   └── summ_ref
│   │   └── ssvs
│   │   │   ├── df_out
│   │   │   ├── time_taken
│   │       └── ssvs_pred_results
│   ├── result_objects_study1.RData
│   └── result_objects_study2.RData
├── master_thesis_5041333.pdf
├── projpred_vs_spikeslab.Rproj
├── README.md
├── renv
├── renv.lock
└── scripts
    ├── data_generation.R
    ├── empirical_examples.qmd
    ├── ssvs_run.R
    ├── study1Analysis.qmd
    ├── study2Analysis.qmd
    ├── workflow_functions.R
    └── workflow_run.R
```

---

## How to Reproduce the Results

### `scripts/`
- **`data_generation.R`**: Generates the simulation datasets.  
  - Running this without changes will create the datasets in `data/p50` (Study 1) and `data/simdata` (Study 2).
- **`workflow_run.R`**: Runs the full PPVS pipeline.
- **`ssvs_run.R`**: Runs the SSVS method.
- **`workflow_functions.R`**: Contains helper functions used in the scripts above.

> **Note:** These scripts contain clearly annotated blocks to indicate which parts to run for Study 1 vs. Study 2. You may need to manually comment/uncomment sections depending on the study you're running.

---

### `output/`
- Simulation results will be saved in `output/output_study1` or `output/output_study2`, depending on the run.
- The full result folders are excluded from the repository due to size limits, but can be available upon request, and key post-processed results are available:
  - `result_objects_study1.RData`
  - `result_objects_study2.RData`

These files contain aggregated results after post-processing and can be used to regenerate all figures without rerunning the simulations.

⚠️ Note on computational cost

- Running the full simulations, especially the PPVS method, can be extremely resource-intensive. Although personal laptops may have many CPU cores, PPVS is highly memory-demanding, which limits how many cores can be used effectively, and may cause RStudio to crash. As a result, the full pipeline can take more than two weeks to complete on a typical personal machine.

- Additionally, the resulting .RData files for each condition and replication can be very large. This is one reason why the raw simulation output is not included in this repository—only post-processed summaries are provided.
---

### Figures & Analysis
- **`study1Analysis.qmd`** and **`study2Analysis.qmd`**: These scripts combine post-processing and plotting. Running them will:
  - Load raw simulation results (if available),
  - Generate post-processed result objects, and
  - Create and save the corresponding plots.
- The code is modular and clearly commented. If you only want to regenerate plots from the result_objects_*.RData files, you can skip the post-processing section and run only the plotting code.
- **`empirical_examples.qmd`**: Reproduces the empirical data analysis and figures.
  - Datasets are publicly available and loaded directly within this script.

---

## Environment

The R environment and package versions used are managed via [renv](https://rstudio.github.io/renv/):
- `renv/`
- `renv.lock`

Running `renv::restore()` will install the correct package versions to reproduce the results.

---

# Ethical approval
Both simulation studies and empirical datasets used in the project are approved by the Ethical Review Board of the Faculty of Social and Behavioural Sciences of Utrecht University.
- Simulation studies: FETC 24-2042 
- Empirical datasets: FETC 25-1375 
  - The approval is valid through 12 May 2025.
    
# License
This repository is licensed under the MIT License.
You are free to use, modify, and distribute this work with appropriate attribution.

# Contact information
For any questions or issues, please contact: Zhipei (Kim) Wang via email: wzpskye@gmail.com.
