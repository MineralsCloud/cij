# <i>C<sub>ij</sub></i>: Semiemperical thermal elasticity

Calculate high temperature thermal elasticity in Python.

## Installation

### Additional information

Pint need to be updated to the GitHub master branch for the code to work with atomic units, so it is installed from GitHub.

## Usage

### Configuration file

```yml
qha:
  input: input01
  settings:
    # QHA INPUT FILE HERE ...
elast:
  input: elast.dat
settings:
  mode_gamma:
    # Refer to cij.core.mode_gamma
    interpolator: spline
    order: 5
  cij_keys: [11, 22, 33, 12, 13, 23, 44, 55, 66]
  disable_phonon_contribution: False
```

### Command line arguments

```txt
usage: qha-cij [-h] [-V] CONF.yml

positional arguments:
  CONF.yml       name of the configuration file

optional arguments:
  -h, --help     show this help message and exit
  -V, --version  show program's version number and exit
```

## Documentation

## Author

The code in this repo is initially authored by Chenxing Luo.

## Licence

Released under GNU GPLv3 license.