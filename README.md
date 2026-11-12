![Tests](https://github.com/dtc-ice-sheets/functional_time_series/actions/workflows/functional_time_series_test.yml/badge.svg)
![Deployment](https://github.com/dtc-ice-sheets/functional_time_series/actions/workflows/functional_time_series_deploy.yml/badge.svg)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/DTC-Ice-Sheets/functional_time_series/blob/main/LICENSE)

# DTC-IS _Functional Time Series_

The `functional_time_series` module is a sub-component of the
[Digital Twin Component for Ice Sheets (DTC-IS)](https://dtc-ice-sheets.org/), an ESA-funded, EO-driven digital twin project.

The `functional_time_series` module reduces a panel of space-time ice-sheet data to a few coefficients per cycle.
Each point's series is cut into curves - melt seasons, typically - and the whole panel is compressed, so that a season
at every location becomes a handful of numbers usable as a predictor alongside other variables. It usually runs on the
output of the [gaussian_random_field](https://github.com/DTC-Ice-Sheets/gaussian_random_field) module.

It was developed to provide hydropower insights for the Watson River catchment on the Greenland ice sheet, where it
compresses each melt season of a supraglacial lake record into a predictor of river discharge.

For an in-depth overview of the functional_time_series module please [read the docs](https://dtc-ice-sheets.github.io/functional_time_series/).

The documentation includes sections on:

1. [Getting started](https://dtc-ice-sheets.github.io/functional_time_series/quick-start.html)
1. [Scientific overview](https://dtc-ice-sheets.github.io/functional_time_series/overview.html)
1. [API reference](https://dtc-ice-sheets.github.io/functional_time_series/api.html)

For a more interactive way of exploring the system, please visit the digital twin's [dashboard](https://dashboards.dtc-ice-sheets.org/)
and/or explore the project [notebook repository](https://github.com/DTC-Ice-Sheets/dtc_is_notebooks).

## Installation for Use

You can install this package with pip:

``` bash
pip install functional_time_series
```

Unlike the `gaussian_random_field` module, this one needs nothing beyond Python - there is no R dependency.

## Installation for Development

``` bash
git clone https://github.com/dtc-ice-sheets/functional_time_series.git
cd functional_time_series
poetry install
poetry run pytest
```
