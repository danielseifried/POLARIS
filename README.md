# POLARIS: POLArized RadIation Simulator

[![arXiv](https://img.shields.io/badge/arXiv-1604.05305-b31b1b)](https://arxiv.org/abs/1604.05305)
[![ascl](https://img.shields.io/badge/ascl-1807.001-262255)](https://ascl.net/1807.001)
[![bibcode](https://img.shields.io/badge/bibcode-2016A%26A...593A..87R-1c459b)](https://ui.adsabs.harvard.edu/abs/2016A&A...593A..87R)
[![doi](https://img.shields.io/badge/doi-10.1051%2F0004--6361%2F201424930-fab70c)](https://doi.org/10.1051/0004-6361/201424930)

is a 3D Monte Carlo radiative transfer code that

- allows to simulate intensity and polarization of light emerging from analytical astrophysical models as well as complex magneto-hydrodynamic simulations on various grids

- is capable to perform dust heating, -emission, -scattering, -grain alignment, line radiative transfer, and synchrotron simulations

- calculates synthetic intensity and polarization maps

- makes use of a full set of physical quantities (density, temperature, velocity, magnetic field distribution, and dust grain properties as well as different sources of radiation) as input


## Requirements

The following packages are required for the installation:

- gcc (preferred), icc, or clang++

- cmake (preferred), or ninja

- python3 (packages: *numpy*, *setuptools*)


## Installation (Linux)

Open a terminal/console and move into the POLARIS directory:
```bash
cd /YOUR/POLARIS/PATH/
```

To install POLARIS, run the installation script:
```bash
./compile.sh -f
```
For the first installation, the option `-f` is required to install the [CCfits](https://heasarc.gsfc.nasa.gov/fitsio/CCfits/) and [cfitsio](https://heasarc.gsfc.nasa.gov/fitsio/) libraries.
For more information, type:
```bash
./compile.sh -h
```
POLARIS can now be executed from any newly opened terminal/console.
However, to use it in already open terminals/consoles, execute the following command to update the environmental paths:
```bash
source ~/.bashrc
```

**HINT**: Please refer to the [manual](manual.pdf) for installation on **macOS**. An installer to use POLARIS with Windows is not available yet.


## Use POLARIS

For a guide how to run first simulations, please take a look in our [quickstart](quickstart.pdf).

For more information about POLARIS and its capabilities, please take a look in our [manual](manual.pdf).

If you use results from POLARIS in a publication, please cite [Reissl et al. (2016)](https://ui.adsabs.harvard.edu/abs/2016A%26A...593A..87R) or [Reissl et al. (2018)](https://ui.adsabs.harvard.edu/abs/2018ascl.soft07001R).
If line radiative transfer and/or Zeeman simulations are used, please cite [Brauer et al. (2017)](https://ui.adsabs.harvard.edu/abs/2017A%26A...601A..90B) as well.


## Copyright

POLARIS is licensed under [GPLv3](LICENSE.md).

Copyright (C) 2018 Stefan Reissl

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

For a list of corresponding authors, please see [AUTHORS](AUTHORS.md).

**contact**: polaris@astrophysik.uni-kiel.de
