# Microscope

Microscope is a package containing a collection of functions which are run as part of the modular workflow to create the RShiny App [Microscope](https://cehrse-boost.datalabs.ceh.ac.uk/).

These functions enable the transformation and combining of spatial data with genomics data.

## Installing

The easiest way to install `microscope` is from GitHub using `devtools`.

In `R` run the following commands:

```
# Install devtools to be able to install packages from GitHub
install.packages("devtools")

# Load devtools and install microscope
library(devtools)
devtools::install_github("NERC-CEH/microscope")
```

## Running

The functions are designed to read in and write out `csv` files.  The documentation should make it clear as to what the structure of the `csv` files should be.  Some functions require additional files, and this is explained in the documentation.

## Running the workflow

The workflow which uses **microscope** can be found here (create repo for notebook) along with documentation on how to run the workflow.

# Authors

Briony Jones - Author / Creator

Contributors - Robin Long
Funder - 

