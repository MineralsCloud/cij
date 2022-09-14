# <i>C<sub>ij</sub></i>: Semiemperical thermal elasticity

Calculate high temperature thermal elasticity in Python.

## Installation

The package can be installed with `pip` package manager.

### Install from PyPI (Python package index)

```shell
$ python3 -m pip install cij
```

### Manual install

At the command prompt, one should navigate to the directory that contains the
`setup.py` script and execute `pip install .`. Then, the package should be ready for use.

## Usage

### Command-line programs

After installation, the Cij program can be started by typing `cij` at your 
command prompt:

```
Usage: cij [OPTIONS] COMMAND [ARGS]...

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  extract           Create data table at specific P or T.
  extract-geotherm  Create data table at geotherm PT.
  fill              Fill non-zero Cij terms based on symmetry.
  modes             Plot interpolated mode frequency vs volume.
  plot              Plot SAM-Cij calculation results.
  run               Perform SAM-Cij calculation.
  run-static        Calculate elastic moduli and acoustic velocities.
```

### SAM-Cij calculations with `cij run`

#### Calculation settings file and example

The `settings.yaml` file is home to all calculation settings. One needs to specify calculation parameters, such as thermal EoS fitting parameters, phonon interpolation settings, input data location, and output variables to store in YAML format. The following is an example settings file.

```yml
qha:
  input: input01
  settings:
    # similar to settings in qha
    DT: 100
    P_MIN: 0
    DELTA_P: 0.5
    NTV: 81
    order: 3
    static_only: False
    T_MIN: 0
    NT: 31
    DT_SAMPLE: 100
    DELTA_P_SAMPLE: 5
    volume_ratio: 1.2
elast:
  input: elast.dat
  settings:
    mode_gamma:
      interpolator: spline
      order: 3
    symmetry:
      system: cubic
output:
  pressure_base:
    - cij
    - vs
    - vp
    - bm_V
    - bm_R
    - bm_VRH
    - G_V
    - G_R
    - G_VRH
    - v
  volume_base:
    - p
    # ...

```

#### Input data

Input data include a QHA input data file (`input01`) and a static elasticity input data (`elast.dat`). See the paper or detailed documentation for description and the [`examples`](./examples) folder for detailed example.

#### Command line arguments

```txt
Usage: cij run [OPTIONS] SETTINGS_FILENAME

  Perform SAM-Cij calculation.

Options:
  --version                       Show the version and exit.
  --debug [CRITICAL|ERROR|WARNING|INFO|DEBUG|NOTSET]
                                  Logging level
  --help                          Show this message and exit.
```


## Structure of the `cij` package

The cij package is written in Python 3. The Python source code is located in the `cij` sub-folder.
Input for three examples in the `examples` sub-folder, documentation in the `docs` sub-folder, and the installation script `setup.py`.

The Python code is organized into several modules:

- **`cij.core`**: Core functionalities
	- `calculator`: The calculator that controls the workflow.
	- `mode_gamma`: Interpolate phonon frequencies and calculate mode Grüneisen parameters.
	- `phonon_contribution`: Calculate phonon *c<sub>ij</sub>*<sup>ph</sup>.
	full_modulus – Interpolate *c<sub>ij</sub>*<sup>st</sup> vs. *V*, and calculate *c<sub>ij</sub><sup>S</sup>* and *c<sub>ij</sub><sup>T</sup>*.
tasks – Handles the ordering of *c<sub>ij</sub>*<sup>ph</sup> calculation.
- **`cij.util`**: Methods used in the main program
	- `voigt`: Convert between Voigt (*c<sub>ij</sub>*) and regular (*c<sub>ijkl</sub>*) notations of elastic coefficients.
	- `units`: Handle unit conversion.
- **`cij.io`**: Input output functions and classes.
- **`cij.plot`**: Plotting functionalities.
- **`cij.cli`**: Command-line programs
	- `cij run` (`main.py`) – Perform a SAM-Cij calculation.
	- `cij run-static` (`static.py`) – Calculate static elastic properties.
	- `cij extract` (`extract.py`) – Extract calculation results for a specific *T* or *P* to a table.
	- `cij extract-geotherm` (`geotherm.py`) – Extract calculation results along geotherm *PT* (given as input) to a table.
	- `cij plot` (`plot.py`) – Convert output data table to PNG plot.
	- `cij modes` (`modes.py`) – Plot phonon frequency interpolation results.
  - `cij fill` (`fill.py`) – Fill all the non-zero terms for elastic coefficients given the constraint of a crystal system.
- **`cij.data`**: Data distributed with the program, e.g., the relationship between *c<sub>ij</sub>*’s for different crystal systems, the naming scheme for output files, etc.
- **`cij.misc`**: Miscellaneous functionalities that are not used in the main program, e.g., methods that facilitate the preparation of input files.

## Author

The code in this repo is initially authored by [Chenxing Luo][1].

[1]: https://github.com/chazeon

## Documentation

See [GitHub pages][2].

[2]: https://mineralscloud.github.io/cij

## Build status

![GitHub Actions](https://github.com/MineralsCloud/cij/actions/workflows/main.yml/badge.svg)
[![codecov](https://codecov.io/gh/MineralsCloud/cij/branch/dev/graph/badge.svg?token=Ln1Fo4vNBE)](https://codecov.io/gh/MineralsCloud/cij)
[![pypi](https://img.shields.io/pypi/v/cij.svg)](https://pypi.org/project/cij/)
[![pypi](https://img.shields.io/pypi/dm/cij.svg)](https://pypi.org/project/cij/)

## How to cite

If you use this software in any publication, please cite:

Luo, C., Deng, X., Wang, W., Shukla, G., Wu, Z., & Wentzcovitch, R. M. (2021). cij: A Python code for quasiharmonic thermoelasticity. *Computer Physics Communications*, 108067. https://doi.org/10.1016/j.cpc.2021.108067

The paper is also available from arXiv: https://arxiv.org/abs/2101.12596

## Licence

Released under [GNU GPLv3](./LICENCE) license.
