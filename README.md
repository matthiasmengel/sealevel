[![DOI](https://zenodo.org/badge/50914429.svg)](https://zenodo.org/badge/latestdoi/50914429)

## Sea level constrained by observations and commitment

Project global sea level rise driven by global-mean-temperature change. Each component contributing to climate-driven sea level rise (thermal expansion, glaciers, the Greenland and the Antarctic ice sheets) is represented by indivdual functions and calibrated
separately.
Total sea level rise is the sum of the components.

### Versions

v.2.0.0: code as used for Mengel et al., Nature Communications 2017, forthcoming.

v 1.0.0: code as used for Mengel et al., PNAS 2016, [https://dx.doi.org/10.1073/pnas.1500515113](https://dx.doi.org/10.1073/pnas.1500515113).

If you make use of this code, please cite the respective paper.

### Making projections
See the [examples/projections.ipynb](examples/projection.ipynb) for projections using a sample global temperature mean timeseries.
If you do not have [jupyter](http://jupyter.org/) notebook, [install](http://jupyter.readthedocs.org/en/latest/install.html) it or look at [do_projections.py](do_projections.py).

### Recalibrating the model
See [do_calibration.py](do_calibration.py).
Model calibration is based on data for past thermal expansion, glaciers and ice caps loss and ice sheet evolution, which is not completely openly available. Please contact M. Mengel or the authors of the datasets applied for calibration if you plan to recalibrate or extend the model. The Antarctic ice sheet contribution following Deconto & Pollard 2016, as introduced in v.2.0.0 is not calibrated within this code, but in the [fast_ant_sid](https://github.com/matthiasmengel/fast_ant_sid) package.

### Installation

`git clone https://github.com/matthiasmengel/sealevel.git`

`cd sealevel`

then with pip

`pip install requirements.txt`

or with conda

`conda env create -f environment.yml`

The model is written in python. If this is new to you, you may have a look at [anaconda python](https://www.continuum.io/downloads). The code has not been tested for Windows. Adaptations to file path definitions would be needed.

### Dependencies

[scipy](http://www.scipy.org/),
[numpy](http://www.numpy.org/),
[dimarray](http://dimarray.readthedocs.org/en/latest/),
[pandas](http://pandas.pydata.org/)
[matplotlib](http://matplotlib.org/) for plotting,
[jupyter](http://jupyter.org/) for convenience.

See [requirements.txt](requirements.txt) for details.

### Issues and contributing

Please post issues on the [issue tracker](https://github.com/matthiasmengel/sealevel/issues) or contact the author. Contributions or propositions to extend the model are most welcome.

### Versions

v.2.0.0: code as used for Mengel et al., Nature Communications 2017, forthcoming.

v 1.0.0: code as used for Mengel et al., PNAS 2016, [https://dx.doi.org/10.1073/pnas.1500515113](https://dx.doi.org/10.1073/pnas.1500515113).

### License

This code is licensed under GPLv3, see the LICENSE.txt. See commit history for authors.
