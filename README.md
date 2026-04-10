# Mapper-Based Localized Prediction with Data-Driven Cover Selection for High-Dimensional Data

## Overview
This repository implements a **Mapper-based localized prediction framework** for high-dimensional data using **Topological Data Analysis (TDA)**. The method leverages the geometric and connectivity structure of data to perform nonparametric prediction via graph-induced neighborhoods.

---

## Key Idea

We model prediction as estimation of conditional probabilities:

$$
P(Y \mid X = x)
$$

using **Mapper graph neighborhoods** instead of standard metric-based neighborhoods.

- Data is projected via a **filter function (OASDA)**
- A **data-driven overlapping cover** is constructed
- **Clustering + graph connectivity** defines local neighborhoods
- Prediction is performed via **localized weighted averaging**

---

## Features

- 📊 Topology-aware prediction  
- 🔬 Designed for high-dimensional biomedical data  
- ⚙️ Data-driven cover selection (bias–variance trade-off)  
- 📈 Works for:
  - Binary outcomes  
  - Nominal classification  
  - Ordinal outcomes  
- 🔍 Permutation-based variable importance  
- 📐 Theoretical guarantees:
  - Consistency  
  - Bayes-risk optimality  

---

## Repository Structure

```text
.
├── Simulation Study/
│   ├── Scenario 1/
│   └── Scenario 2/
├── Real Data Analysis/
│   ├── PPMI/
│   └── UCSC Xena/
├── Images/
│   ├── PPMI.png
│   └── Final_UCSC_Xena.png
└── README.md

```

---

## Method Summary

### 1. Mapper Construction
- Filter: OASDA projection  
- Cover: overlapping intervals  
- Clustering: hierarchical (Ward’s method)  
- Output: Mapper graph  

---

### 2. Prediction Rule

For a new point:

- Identify relevant intervals  
- Assign to nearest clusters  
- Aggregate over graph-connected nodes  
- Compute weighted class probabilities  

---

### 3. Cover Selection

Optimal scaling:

- Number of intervals: `l ~ n^(1/5)`  
- Observations per interval: `S ~ n^(4/5)`  

This balances:
- Bias (oversmoothing)  
- Variance (instability)  

---

## Simulation Study

We evaluate performance under:

- Nonlinear branching structures  
- Skewed distributions (log-normal)  
- Correlated predictors  
- Varying noise levels  
- Sample sizes: 250, 500, 1000  

### Key Findings

- Outperforms classical models under heterogeneous structure  
- Comparable performance under homogeneous settings  
- Strong robustness to noise  

---

## Applications

### 1. Parkinson’s Disease (PPMI)
Predict Hoehn–Yahr stage  

### 2. TCGA Glioma (UCSC Xena)
Classify tumor severity from RNA-seq data  

## Visualization

### UCSC Xena (TCGA Glioma Data)

<p align="center">
  <img src="Images/Final_UCSC_Xena.png" width="700"/>
</p>

This figure illustrates Mapper-derived structure in high-dimensional RNA-seq data.


---

## Requirements

- R (recommended)

---

