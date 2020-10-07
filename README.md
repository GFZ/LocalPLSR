# LocalPLSR
**Summary of the local PLSR methodology**
Quantification of soil parameters e.g. organic carbon, in field soil samples using spectrally similar soil samples from a large scale soil spectral database. For this purpose, Laboratory spectra measured from the field soil samples are required and compared with the laboratory spectra from the soil spectral database. The most similar soil samples from the soil spectral database are then used to train a separate PLSR model for each field sample and apply this model to the respective field sample to quantify the soil parameter. The localPLSR thus represents a reproducible, simplified and thus somewhat less accurate alternative to the conventional quantification (in the chemical laboratory) of some (spectrally active) soil parameters. For this purpose, laboratory spectra of the field samples must be measured with laboratory spectrometers.

**Overview Script**
* Data preparation including soil spectral library (calibration) and field samples (validation)
    * Examples for smoothing and frist derivative
    * Examples of PCA and Mahalanobis distance
* LocalPLSR function
* Application example of localPLSR function