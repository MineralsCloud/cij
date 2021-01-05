# <i>C<sub>ij</sub></i>: Semiemperical thermal elasticity

Calculate high temperature thermal elasticity in Python.

## Installation

The package can be installed with `pip` package manager.

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

And are avaliable as standalone commands:

- **`cij-run`**         - Perform SAM-Cij calculation.
- **`cij-run-static`**  - Calculate elastic moduli and acoustic velocities.
- **`cij-extract`**     - Extract the value and create data table for multiple variables at geotherm PT.
- **`cij-extract-geotherm`** - Extract the value and create data table for multiple variables at spcific P or T.
- **`cij-fill`**        - Fill non-zero Cij terms based on symmetry.
- **`cij-modes`**       - Plot interpolated mode frequency vs volume.
- **`cij-plot`**        - Plot SAM-Cij calculation results.


### SAM-Cij calculations with `cij run`

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

[1]: https://mineralscloud.github.io/cij

## Author

The code in this repo is initially authored by Chenxing Luo.

## Licence

Released under GNU GPLv3 license.