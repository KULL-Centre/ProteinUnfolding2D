ProteinUnfolding2D
------------------

Python library for thermodynamic analysis of protein unfolding using both temperature and chemical denaturant.

Current version is made for data from NanoTemper Prometheus 

### Examples

The examples directory contains a number of use cases that shows how to use the library

- 20191202_ACBP.ipynb : Jupyter notebook used in the original description of the library [1]

### Dependencies

- Python version >= 3.6
- numpy
- xlrd (only version 1.2 reads .xlsx files. For newer versions, convert to .xls using Excel)
- lmfit version 1.0.1
- copy

### References

1. Louise Hamborg, Emma Wenzel Horsted, Kristoffer Enøe Johansson, Martin Willemoës, Kresten Lindorff-Larsen & Kaare Teilum (2020) "Global analysis of protein stability by temperature and chemical denaturation" Analytical Biochemistry 605, 113863. https://doi.org/10.1016/j.ab.2020.113863
