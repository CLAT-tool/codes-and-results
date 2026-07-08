# CLAT-tool: Codes and Results

This repository contains the experimental codes and results for the paper:

**A Generalized Framework for Conditional Linear Cryptanalysis and Its Application to AES-Like Ciphers**

The repository provides supporting materials for the proposed CLAT-based conditional linear cryptanalysis framework, including correlation evaluation, empirical correlation verification, and the experimental verification of the key recovery attack accelerated by the Fast Walsh-Hadamard Transform (FWHT).

## Repository Organization

This repository contains three folders. They are arranged according to the order in which they are referenced in the paper.


## Folder Description

### 1. FWHT Verification

This folder contains the experimental verification corresponding to the key recovery attack based on the FWHT technique described in Section 4.4. The specific details are provided in Appendix C.


### 2. Correlation Evaluation

This folder contains the codes and data used for evaluating the theoretical correlations of conditional linear trails.

The purpose of this folder is to compute or summarize the expected conditional linear correlations derived from the CLAT-based framework. These evaluations are used to estimate the quality of the discovered conditional linear distinguishers and to compare them with standard linear distinguishers.

The results in this folder support the theoretical evaluation of conditional linear trails for ciphers such as AES, Rijndael-256, ARIA, LED, Midori-128, and SKINNY-128.

### 3. Correlation Verification

This folder contains the experimental verification of the correlations of selected conditional linear trails.

The purpose of this folder is to empirically validate the correlations predicted by the theoretical evaluation. The experiments generate plaintext-ciphertext pairs, filter or partition the data according to the specified conditions, and compute the observed empirical correlations.


The verification results provide experimental evidence that the discovered conditional linear trails exhibit the expected statistical behavior.

## Requirements

The experiments are implemented for research and reproducibility purposes. The exact requirements may vary across folders, but a typical environment includes:

* C++;
* Gurobi;
* OpenMP
