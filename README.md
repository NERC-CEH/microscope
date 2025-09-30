# Microscope

- Extention and reimplementation of [ID-TaXER](https://catalogue.ceh.ac.uk/documents/3da8144d-519f-4b93-b375-60dd49257cc9)
- The RShiny app which relies on this code is [here](https://cehrse-boost.datalabs.ceh.ac.uk/)
- The code to produce the RShiny app is located in [microscope/vignette/app.R](https://github.com/NERC-CEH/microscope/blob/main/vignette/app.R)
- The main code is in [microscope/R/analysis.R](https://github.com/NERC-CEH/microscope/blob/main/R/analysis.R) which is used to manipulate the data which goes into the [RShiny app](https://cehrse-boost.datalabs.ceh.ac.uk/)

## Description

Microscope is a package containing a collection of functions which are run as part of the modular workflow to create the RShiny App [Microscope](https://cehrse-boost.datalabs.ceh.ac.uk/).

These functions enable the transformation and combining of spatial data with genomics data.

## Installing

The easiest way to install `microscope` is from GitHub using `devtools`.

In `R` run the following commands:

```R
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

## Authors

Authors: Briony Jones

Contributors: Robin Long

## Funding Statement
This projected is part of the [Environmental Data Service (NERC EDS)](https://eds.ukri.org/environmental-data-service) pilot, funded under the UKRI Research Cloud Pilot Initiative, Digital Research Infrastructure programme.
