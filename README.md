<!-- Copyright (C) 2020
- Kathrin J. Ward (GFZ, kathrin.ward@gfz-potsdam.de), 
- Sabine Chabrillat (GFZ, sabine.chabrillat@gfz-potsdam.de),
- Saskia Foerster (GFZ, saskia.foerster@gfz-potsdam.de), 
- Helmholtz Centre Potsdam, German Research Centre for Geosciences (GFZ, https://www.gfz-potsdam.de/startseite/)

This program was developed within the context of the following publicly funded projects:
- EnMAP scientific preparation program under the DLR Space Administration, German Federal Ministry of Economic Affairs and Energy, 
  grant number 50EE1529 (https://www.enmap.org/)
 
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 2 of the License, complemented with 
the following provision: 
For the scientific transparency and verification of results obtained 
and communicated to the public after using a modified version of the 
work, You (as the recipient of the source code and author of this 
modified version, used to produce the published results in scientific 
communications) commit to make this modified source code available in 
a repository that is easily and freely accessible for a duration of 
five years after the communication of the obtained results.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>. -->


# LocalPLSR
**Summary of the local PLSR methodology:**
Quantification of soil parameters e.g. organic carbon, in field soil samples using spectrally similar soil samples from a large scale soil spectral database. For this purpose, Laboratory spectra measured from the field soil samples are required and compared with the laboratory spectra from the soil spectral database. The most similar soil samples from the soil spectral database are then used to train a separate PLSR model for each field sample. Each model is applied to the respective field sample to quantify the soil parameter. The localPLSR thus represents an alternative to the conventional quantification (in the chemical laboratory) of spectrally active soil parameters. For this purpose, laboratory spectra of the field samples need to be measured with laboratory spectrometers.

**Overview Script**
* Data preparation including soil spectral library (calibration) and field samples (validation)
    * Examples for smoothing and frist derivative
    * Examples of PCA and Mahalanobis distance
* LocalPLSR function
* Application example of localPLSR function

**How to run this script**
* Open and run this script in R
* Ensure to install the required R libraries prior to running the script
* Further user specific adaptations of the script might be necessary of which some are indicated by comments 
* The examples given for data preparation and application of the localPLSR function should be adapted to your data