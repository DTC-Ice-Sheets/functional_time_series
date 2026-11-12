"""Entrypoint for FTS module."""

import logging
from pathlib import Path

import cyclopts
import xarray as xr

from functional_time_series.fts_spatial_module.fts_spatial_module import FtsSpatialModule
from functional_time_series.helpers import ZarrType, save_geo_zarr
from functional_time_series.monitoring import setup_logging

logger = logging.getLogger(__name__)

app = cyclopts.App()


@app.default
def main(
    data_path: ZarrType,
    output_path: Path,
    curve_time_points: int = 5,
    num_pca: int = 2,
    num_factors: int = 3,
) -> None:
    """Run the Functional Time Series module.

    Parameters
    ----------
    data_path : ZarrType
        Path of the input GeoZarr file.
    output_path : Path
        Path of the output GeoZarr file.
    curve_time_points : int
        Number of time points to be considered as a single curve
    num_pca : int
        Number of principal components to retain. Only the first feeds the second step, so this does not
        change the result.
    num_factors: int
        Number of factors to retain in the second step dimension reduction
    """
    setup_logging()

    fts_module = FtsSpatialModule()
    ds = xr.load_dataset(data_path, engine="zarr", chunks={})

    output = fts_module.run_spatial_analysis(
        ds,
        curve_time_points=curve_time_points,
        num_pca=num_pca,
        num_factors=num_factors,
    )

    save_geo_zarr(output, output_path)


if __name__ == "__main__":
    app()
