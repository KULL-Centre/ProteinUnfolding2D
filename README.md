ProteinUnfolding2D
------------------

Python library for thermodynamic analysis of protein unfolding using both temperature and chemical denaturant.

Current version is made for data from NanoTemper Prometheus, see fork by rnorrild for NanoTemper Panta version

Hints for usage:
- Remember to use the 'processed' output from the NanoTemper and to put denaturant concentrations in the spreadsheet
- The first part of the notebook selects the sample and temperature range to use. Start by selecting all and look at the plots to identify outliers
- Some of the initial parameters are guessed, e.g. the heat capasity is guessed from the number of amino acids
- If the fit fails, try to obtain better initial parameters or remove outliers at extreem temperatures or samples that look wird. To improve the initial parameters, try to scan dH, comment out guessing of heat capacity and m-value or comment-in the code that plots the initial base-plane parameters (baselines in temperatuer and denaturent)

### Examples

The examples directory contains a number of use cases that shows how to use the library

- 20191202_ACBP.ipynb : Jupyter notebook used in the original description of the library [1]

### Dependencies

- Python version 3.6
- numpy
- xlrd  version 1.2 (only version 1.2 reads .xlsx files. For newer versions, convert to .xls using Excel)
- lmfit version 0.9.x - 1.0.1
- copy
- matplotlib version 3.3 (other versions <3.5 may work as well)

Since the code now depends on older packages, it may be useful to setup a closed "python ecosystem" using e.g. Anaconda or Enthought

#### Python environment using Enthought Deployment Manager (EDM) command line (CLI)

Doc: http://docs.enthought.com/edm/user/cli.html

````bash
# install python 3.6 environment
edm environments create ProteinUnfolding2D --version 3.6
edm install "xlrd == 1.2.0-1" --environment ProteinUnfolding2D
edm install "lmfit == 0.9.2-14" --environment ProteinUnfolding2D
edm install "matplotlib == 3.3.4-3" --environment ProteinUnfolding2D
edm install jupyter --environment ProteinUnfolding2D

# activate environment in default shell
edm shell --environment ProteinUnfolding2D

# first time you may need to enable widgets to see the slider in denaturent plots
jupyter nbextension enable --py widgetsnbextension

# open notebook in default browser
jupyter notebook examples/ACBP_unfolding2d.ipynb
```


### References

1. Louise Hamborg, Emma Wenzel Horsted, Kristoffer Enøe Johansson, Martin Willemoës, Kresten Lindorff-Larsen & Kaare Teilum (2020) "Global analysis of protein stability by temperature and chemical denaturation" Analytical Biochemistry 605, 113863. https://doi.org/10.1016/j.ab.2020.113863
