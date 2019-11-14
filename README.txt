
Michalis Kassinopoulos, PhD candidate
Graduate Program in Biological and Biomedical Engineering, McGill University, Montreal, Canada
mkassinopoulos@gmail.com

Date: 11-Sept-2019

==================================================================

This repository provides scripts used for our work presented in:

Kassinopoulos, M., Mitsis, G.D., 2019. White Matter Denoising Improves the Identifiability of Large-Scale Networks and Reduces the Effects of Motion in fMRI Functional Connectivity. bioRxiv https://doi.org/10.1101/837609

Abstract: It is well established that confounding factors related to head motion and physiological processes (e.g. cardiac and breathing activity) should be taken into consideration when analyzing and interpreting results in fMRI studies. However, even though recent studies aimed to evaluate the performance of different preprocessing pipelines there is still no consensus on the optimal strategy. This may be partly because the quality control (QC) metrics used to evaluate differences in performance across pipelines often yielded contradictory results. Importantly, noise correction techniques based on physiological recordings or expansions of tissue-based techniques such as aCompCor have not received enough attention. Here, to address the aforementioned issues, we evaluate the performance of a large range of pipelines by using previously proposed and novel quality control (QC) metrics. Specifically, we examine the effect of three commonly used practices: 1) Removal of nuisance regressors from fMRI data, 2) discarding motion-contaminated volumes (i.e., scrubbing) before regression, and 3) low-pass filtering the data and the nuisance regressors before their removal. To this end, we propose a framework that summarizes the scores from eight QC metrics to a reduced set of two QC metrics that reflect the signal-to-noise ratio (SNR) and the reduction in motion artifacts and biases in the preprocessed fMRI data. Using resting-state fMRI data from the Human Connectome Project, we show that the best data quality, is achieved when the global signal (GS) and about 17% of principal components from white matter (WM) are removed from the data. In addition, while scrubbing does not yield any further improvement, low-pass filtering at 0.20 Hz leads to a small improvement.

====================================================================

The User guide will be available shortly. It will explain what each script is related to and how to reproduce some of the results presented in the study.

Please do not hesitate to contact me if you have any questions related to the use of these scripts.




