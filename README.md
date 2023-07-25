# Delineating the Interplay Between Oncogenic Pathways and Immunity in Anaplastic Wilms Tumors: Implications for Prognosis and Therapeutic Vulnerability

This repository contains the analytical pipeline for our research project. We aimed to unravel the complex interplay between oncogenic pathways and immune responses in Anaplastic Wilms Tumors and derive implications for prognosis and potential therapeutic vulnerabilities.

## Prerequisites

The analysis was performed in the following environment:

- R version 4.2.2 (2022-10-31 ucrt)
- Platform: x86_64-w64-mingw32/x64 (64-bit)
- Running under: Windows 10 x64 (build 19044)

## Project Structure

### Main Script

- `main.R`: This is the central script containing the entire unstructured code used in this analysis.

### Custom Scripts Directory

The `/scripts` directory contains customized functions used in this project:

- `annTrackScale.R`: Scales and truncates continuous data for annotation in the heatmap.
- `batchPCA.R`: Checks if the batch effect has been properly removed.
- `calFPKM.R`: Calculates the FPKM value from raw count matrix.
- `creat.anntrack.R`: Generates corresponding colors for variables in factors that can be used in heatmap annotation.
- `display.progress.R`: Displays the progression in a 0-100 percentage manner in for loops.
- `generateInputFileForSubMap.R`: Generates input files for subclass mapping analysis.
- `LowExpFilter.R`: Filters out genes with low expression across the samples.
- `normData.R`: Normalizes raw count data.
- `plot.common.cluster.R`: Performs consensus hierarchical clustering with all the intermediate inputs saved.
- `standarize.fun.R`: Scales and truncates the data for the main part of the heatmap.
- `twoclassedgeR.R`: Performs differential expression analysis between two conditions using the edgeR algorithm based on raw count matrix.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

