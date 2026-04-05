# Real Data Application: Proteogenomics

This folder contains the scripts required to reproduce the real-world data analysis (Section 5) from the manuscript. The application predicts Estrogen Receptor Alpha (ER-$\alpha$) protein expression using high-dimensional RNA-seq transcriptomic data from the TCGA Breast Cancer (BRCA) dataset.

## Files

* `run_real_data.R`: The master script for the application. It automatically downloads the matched RNA-seq and proteomic data, applies an unsupervised variance filter to retain $p=300$ genes, injects artificial concentrated leverage points into 15% of the patients, and benchmarks the algorithms.

## Dependencies

Because this script pulls standardized, high-dimensional multi-assay genomic data directly from The Cancer Genome Atlas (TCGA), it requires specific packages from **Bioconductor** in addition to the standard CRAN packages.

**1. Install CRAN packages:**

```R
install.packages("remotes")
remotes::install_version("robustbase", version = "0.99.6")
remotes::install_version("RobStatTM", version = "1.0.11")
```

**2. Install Bioconductor packages:**

```R
install.packages("remotes")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

remotes::install_version("BiocManager", version = "1.30.27")
remotes::install_version("curatedTCGAData", version = "1.32.1", repos = BiocManager::repositories())
remotes::install_version("TCGAutils", version = "1.30.2", repos = BiocManager::repositories())
```

**3. Plotting and Formatting Packages:**

```R
remotes::install_version("ggplot2", version = "4.0.1")
remotes::install_version("tidyr", version = "1.3.2")
remotes::install_version("xtable", version = "1.8.4")
remotes::install_version("scales", version = "1.4.0")
```

## Usage

Ensure your working directory is set to the root of the repository. You can source the script directly from your R console:

```R
source("application/run_real_data.R")
```

*Note: Upon first execution, the script will create a `data/` folder and download the TCGA datasets (this may take a few minutes depending on your internet connection). Subsequent runs will load the data locally. Results and plot data will be saved to a `results/` sub-directory.*
```








