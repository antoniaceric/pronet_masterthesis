# Master's Thesis PRONET

This repository contains the analysis scripts for the PROENT project, investigating categorical speech perception and the distribution of
phonemic and prosodic representations across the bilateral frontotemporal language network in a combined TMS-fMRI approach using Representational Similarity Analysis (RSA).

This repository only contains the main scripts and custom function needed for the analyses for the thesis part of the project. Additional versions, variations and other analyses (such as other pipeline versions, searchlight(wip) or contrast analysis) can be found in the following repository: [git@github.com:seunggookim/pronet.git](https://github.com/seunggookim/pronet)


---

## Core Scripts and Functions of Pronet Master's Thesis

* **`pronet_main_may2026.m`**  
  The central analysis script. Computes multivariate neural distance matrices (Crossnobis), runs RSA modeling (NNLS / linear regression), and performs variance partitioning across experimental conditions.
* **`pronet_main_pretestbehav.m`**  
  Constructs theoretical and participant-specific model Representational Dissimilarity Matrices (RDMs) - ingredients for RSA

* **`pronet_main_may2026.m`**  
  The central analysis script. Computes multivariate neural distance matrices (Crossnobis), runs RSA modeling (NNLS / linear regression), and performs variance partitioning across experimental conditions + noise ceiling calculation.


> **Note:** Modular helper functions (e.g., first-level GLM estimation, Crossnobis distance calculations, permutation testing, and plotting) are called automatically by the main scripts and are documented directly within their respective files.
---
## Overview methodological pipelines
Pipelines are automatically run by `pronet_main_may2026.m`

<img width="1230" height="335" alt="pipeline_overview_082026" src="https://github.com/user-attachments/assets/2838e172-fdd4-4c0d-9393-e342cb0cd1e0" />


---

## Standard Workflow

1. **Model RDM Construction:** Run `pronet_main_pretestbehav.m` to generate predictor matrices from psychometric fits.
2. **Behavioral QA & Psychometrics:** Run `pronet_exclusion_behavioral.m`, `pronet_updatedsigmoidvslinear_290526.m`.
3. **Representational Similarity Analysis:** Execute `pronet_main_may2026.m` for neural RDM computation, RSA modeling, and statistical inference.

---

## Dependencies & Atlases

Add the following external resources to your MATLAB path:

* **RSA Toolbox for MATLAB**  
  [`rsagroup/rsatoolbox_matlab`](https://github.com/rsagroup/rsatoolbox_matlab)  
  *Core framework for representational similarity analysis and RDM construction.*
* **HCPex Atlas**  
  [`wayalan/HCPex`](https://github.com/wayalan/HCPex)  
  *Extended Human Connectome Project multimodal parcellation used for ROI definitions.*

---

## Requirements

* **MATLAB** (R2024b or later recommended)
* **SPM12** ([Statistical Parametric Mapping](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/))
