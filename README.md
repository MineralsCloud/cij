# <i>C<sub>ij</sub></i>: Semiemperical thermal elasticity

Calculate high temperature thermal elasticity in Python.

## Installation

The package can be installed with `pip` package manager.

## Usage

### Command-line programs

The packge is shipped with the following command-line programs:

- **`qha-cij`**         - Perform SAM-Cij calculation.
- **`qha-cij-static`** - Calculate elastic moduli and acoustic velocities.
- **`qha-cij-extract`**     - Extract the value and create table for multiple variables at spcific P or T.
- **`qha-cij-fill`**        - Fill non-zero Cij terms based on symmetry.
- **`qha-cij-modes`**       - Plot interpolated mode frequency vs volume.
- **`qha-cij-plot`**        - Plot SAM-Cij calculation results.

And is avaliable as sub-commands of the `cij` program:


```
Usage: cij [OPTIONS] COMMAND [ARGS]...

Options:
  --version  Show the version and exit.
  --help     Show this message and exit.

Commands:
  extract     Extract the value and create table for multiple variables at...
  fill        Fill non-zero Cij terms based on symmetry.
  modes       Plot interpolated mode frequency vs volume.
  plot        Plot SAM-Cij calculation results.
  run         Perform SAM-Cij calculation.
  run-static  Calculate elastic moduli and acoustic velocities.
```

### SAM-Cij calculations with `qha-cij` or `cij run`

#### Configuration file

```yml
qha:
  input: input01
  settings:
    # QHA settings here ...
elast:
  input: elast.dat
  settings:
    mode_gamma:
      interpolator: spline
      order: 5
    symmetry: cubic
```

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

## Documentation

See [GitHub pages][1].

[1]: https://mineralscloud.github.io/qha-cij

## Author

The code in this repo is initially authored by Chenxing Luo.

## Licence

Released under GNU GPLv3 license.