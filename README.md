# Carnesecchi_et_al._2021
This repository supports Carnesecchi et al., 2021

Image series used consists of singular S2R+ culture cells, acquired on a NIKON A1R with NIS Elements. The Jython script accepts single cells in frame with fragments of other cells on the image borders, which will be ignored and excluded from analysis.

Prerequisites to be able to run the Jython script:

1. Have the slice alignment plugin from Tseng, Q. et al. (2012) installed on ImageJ, or adapt code to aligment tool of preference.
2. Have the stimulus ROI in the same directory as the timeseries.

Tips for running the R script:

1. In case the fitting of the model does not converge, e.g. "singular gradient matrix at initial parameter estimates", try different starting values.
