# Ecological and Evolutionary Drivers of 2.3.4.4b H5Nx HPAI Spread in Europe

Sachin Subedi¹˒², M. H. M. Mubassir¹˒², Tanin Rajamand¹˒², Mohamed Bakheet²˒³,Ludy Registre Carmola²˒³, Leke Lyu⁵, Oluwatosin Babasola²˒³, Sihua Peng²˒³, Yangfan Liu⁶, Carsten Kirkeby⁶ and Justin Bahl¹˒²˒³˒⁴\*

¹ Institute of Bioinformatics, University of Georgia, Athens, GA, USA\
² Center for Ecology of Infectious Diseases, University of Georgia, Athens, GA, USA\
³ College of Veterinary Medicine, Department of Infectious Diseases, University of Georgia, Athens, GA, USA\
⁴ College of Public Health, Department of Epidemiology and Biostatistics, University of Georgia, Athens, GA, USA\
⁵ Rollins School of Public Health, Department of Biostatistics and Bioinformatics, Emory University, Atlanta, GA, USA\
⁶ Department of Veterinary and Animal Sciences, University of Copenhagen, Copenhagen, Denmark

\*Corresponding author: **Justin Bahl** — `justin.bahl@uga.edu`

------------------------------------------------------------------------

## Abstract

<small> Highly pathogenic avian influenza (HPAI) H5Nx clade 2.3.4.4b has been circulating in Europe since 2016, posing a threat to poultry, wildlife, and public health. Despite the ongoing economic and pandemic threat, the determinants of viral spillover, spread, and annual emergence remain poorly understood. To address this gap, we integrated genomic, ecological, and epidemiological data from 2016 to 2025 into a comparative genetic framework, allowing us to quantify viral spread among regions, and habitats. We analyzed 7,031 European H5Nx HA sequences, integrating habitat trait data from AVONET. Geographic strata were defined through a dual classification framework, hierarchical clustering of country centroids using Haversine distances cross-classified with European Environment Agency (EEA) biogeographic regions, yielding 10 ecologically and geographically grounded regional strata. Three subsampling strategies (stratified, equal, and epidemic-weighted proportional) yielded three phylogeographic datasets for robust inference. Our analyses show that the Atlantic region served as a persistent reservoir for the long-term circulation of H5Nx viruses. The Western Continental region was a key secondary hub that structured cross-regional transmission. Wetlands formed the backbone of transmission, maintaining long-term circulation before repeated jumps into farms, grasslands, and coastal areas. Anatidae migration and live-poultry trade were important predictors for inter-regional transmission. Together, these results show that European H5Nx dynamics are structured by regional transmission hubs and persistent wetland reservoirs, highlighting waterfowl migratory connectivity and trade corridors as priorities for targeted surveillance and control. </small> ---

## Repository structure

### `Pipelines/`

Data filtering, stratified subsampling, and the step-by-step workflow used to produce the inputs for downstream phylodynamic and GLM analyses.

-   `Pipeline.Rmd`\
    This repository follows a reproducible pipeline that turns raw sequence/metadata into analysis-ready inputs for BEAST discrete-trait phylogeography and GLM models. We retrieved all available European H5Nx hemagglutinin (HA) nucleotide sequences and metadata from GISAID, integrated wild-bird ecological information from AVONET, and assigned discrete traits (habitat and host type) using isolate/host metadata with context-specific rules (e.g., Farm for agricultural production species such as poultry, Urban for human cases and residential domestic mammals such as pet cats). To define geographically coherent regions, we built 10 data-driven Regions in R: Atlantic, Western Continental, Eastern Continental, Eastern Alpine, Central Alpine, Pannonian, Boreal Baltic, Scandinavian Highlands, Southeast Mediterranean, and Iberian.
-   `Subsampling.Rmd`\
    Generates three subsamples (equal, propotional and stratified) used across downstream BEAST analyses.

------------------------------------------------------------------------

### `Phylodynamics/`

#### 1) `Empirical Trees/`

Empirical-tree BEAST runs and utilities for combining/processing posterior trees.

-   `xml/`\
    BEAST XMLs for empirical-tree inference, one per subsample.
-   `Scripts/`\
    Utilities for running BEAST, combining runs/trees, and preparing tree sets for downstream analyses.

------------------------------------------------------------------------

#### 2) `Discrete Trait Analysis/`

Discrete-trait phylogeographic analyses (e.g., Region, Habitat, and Combined traits).

-   `xmls/`\
    Ready-to-run BEAST XML configurations.\
    Organized by `equal/`, `proportional/`, `stratified/` to reflect independently generated balanced datasets.

-   `Scripts/`\
    Post-processing and visualization scripts for BEAST outputs:

    -   `Beast/`\
        Shell scripts to summarize outputs, combine runs, and organize post-BEAST processing.
    -   `Markov Jumps/`\
        R scripts to summarize Markov jump counts across subsamples and trait models.
    -   `Rewards/`\
        Python scripts to compute and visualize Markov rewards across traits.
    -   `Transition Rates/`\
        R scripts to compute transition-rate summaries and Bayes factors; includes mapping/chord/summary plots.
    -   `Trees/`\
        Python scripts for reading posterior/MCC trees and extracting trait dynamics.

------------------------------------------------------------------------

#### 3) `Generalized Linear Modeling/`

Predictor-based modeling of transition rates using a GLM framework, including predictor data, preparation scripts and visualization.

-   `Scripts/Preparation/`\
    Scripts to build predictor layers used in GLM analyses (examples: Anatidae migration, climate variables, land cover, poultry density, road-transported goods, live poultry trade matrices, wild bird counts, sampling effort).
-   `Scripts/PCA and Correlation/`\
    Collinearity checks and dimensional summaries used to define non-collinear predictor sets.
-   `Scripts/Plots/`\
    R scripts that produce GLM result plots and Bayes factor summaries.
-   `xmls/`\
    BEAST-GLM XMLs organized by `equal,proportional and stratifed`, with multiple predictor sets.

------------------------------------------------------------------------
