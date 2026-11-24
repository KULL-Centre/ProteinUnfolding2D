ProteinUnfolding2D
------------------

Python library for thermodynamic analysis of protein unfolding using both temperature and chemical denaturant.

Hints for usage:
- Remember to use the 'processed' output from the NanoTemper and to put denaturant concentrations in the spreadsheet
- The first part of the notebook selects the sample and temperature range to use. Start by selecting all and look at the plots to identify outliers
- Some of the initial parameters are guessed, e.g. the heat capasity is guessed from the number of amino acids
- If the fit fails, try to obtain better initial parameters or remove outliers at extreem temperatures or samples that look wird. To improve the initial parameters, try to scan dH, comment out guessing of heat capacity and m-value or comment-in the code that plots the initial base-plane parameters (baselines in temperatuer and denaturent)

### Examples

The examples directory contains a number of use cases that shows how to use the library

- 20191202_ACBP.ipynb : Jupyter notebook used in the original description of the library [1]
- Sin3bPAH1_unfolding2d.ipynb : Kind of a negative example used to test output files from 2026 or newer
- PantaExample.ipynb : Test for NanoTemper Penta output. Credit to [Rasmus Norrild!](https://github.com/rnorrild/ProteinUnfolding2D)

### Dependencies

conda create -p ./nanotemper_env python==3.13 numpy xlrd openpyxl lmfit matplotlib jupyter

### References

1. Louise Hamborg, Emma Wenzel Horsted, Kristoffer Enøe Johansson, Martin Willemoës, Kresten Lindorff-Larsen & Kaare Teilum (2020) "Global analysis of protein stability by temperature and chemical denaturation" Analytical Biochemistry 605, 113863. https://doi.org/10.1016/j.ab.2020.113863
