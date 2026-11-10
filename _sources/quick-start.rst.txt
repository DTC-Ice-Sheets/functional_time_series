Quick Start
===========

The ``functional_time_series`` package provides the
``functional_time_series.fts_spatial_entrypoint`` command-line entrypoint.
It reads a GeoZarr dataset containing a time series, fits the spatial
functional time series model, and writes the resulting curve coefficients as
a GeoZarr dataset.

Requirements
------------

Use Python 3.12 or newer, and install the project dependencies with Poetry:

.. code-block:: bash

   poetry install

Run the command from the repository root. For a list of arguments and their
defaults:

.. code-block:: bash

   poetry run python -m functional_time_series.fts_spatial_entrypoint --help

Run the spatial FTS model
-------------------------

The required arguments are the input data path and output path. The input
should be a GeoZarr time series with a ``spatial_ref`` coordinate containing a
valid ``crs_wkt`` attribute. Spatial data may be gridded (``x`` and ``y``
dimensions) or represented by a ``point_id`` dimension. For example, the
repository includes a GRF output GeoZarr that can be used as input:

.. code-block:: bash

   mkdir -p output
   poetry run python -m functional_time_series.fts_spatial_entrypoint \
      --data_path tests/test_data/test_inputs/watson_runoff_racmo55_grf_output.zarr \
      --output_path output/watson_runoff_racmo55_fts.zarr

The input path may be a local path or a URL for a publicly available Zarr store. The
output is written to a local GeoZarr store. The FTS output contains one
variable for each analyzed input variable, with dimensions ``curve_id`` and
``coef_id``, plus the spatial reference.

Options
-------

The following options are optional:

``--curve_time_points``
   Number of time steps in each curve (default: ``5``). Any trailing time
   steps that do not make a whole curve are omitted.
``--num_pca``
   Number of principal components to retain in the first dimension-reduction
   step (default: ``2``).
``--num_factors``
   Number of factors, and therefore output coefficients, to retain in the
   second dimension-reduction step (default: ``3``).

For example, to use 8 time steps per curve and retain 6 factors:

.. code-block:: bash

   poetry run python -m functional_time_series.fts_spatial_entrypoint \
      --data_path tests/test_data/test_inputs/watson_runoff_racmo55_grf_output.zarr \
      --output_path output/watson_runoff_racmo55_fts_custom.zarr \
      --curve_time_points 8 \
      --num_pca 4 \
      --num_factors 6
