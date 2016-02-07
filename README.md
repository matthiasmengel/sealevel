Source code repository for

## Sea level constrained by observations and commitment

If you make use of the code, please cite

```
Future sea level rise constrained by observations and long-term commitment
M. Mengel et al. PNAS 2016
```

This repository contains the code to produce the results as presented in the paper.
There are two main parts: The calibration of the model and using the model for projections.

### Making projections
See the example.ipynb for projections using a sample global temperature mean timeseries.
If you do not have jupyter notebook, [install](http://) it or look at example.py.

### Recalibrating the model
See `src/run_calibration.py`.
Model calibration is based on data for past thermal expansion, glacier loss and ice sheet evolution, that is not completely openly available. Please contact the M. Mengel or the authors of the datasets applied for calibrition if you plan to recalibrate or extend the model.

### Installation

`git clone XX`

`cd sealevel`,
then with pip `pip install requirements.txt` or with conda `conda install requirements.txt`

The model is written in python. If this is new to you, you may have a look at anaconda python, which conveniently provides the packages to run this code. The model has not been tested under Windows.

### Dependencies

[scipy](http://), 
[numpy](http://), 
[dimarray](http://), 
[matplotlib](http://) for plotting, 
[jupyter](http://) for convenience.

### Issues and contributing

Please post issues on the [issue tracker](http://) or contact the author.

### Versions

v 1.0: code as used for the Mengel et al., PNAS 2016.

### License

This code is licensed under GPLv3, see the LICENSE.txt.


