"""Utility functions for use across functional_time_series modules."""

import logging
from pathlib import Path
from typing import Annotated

import numpy as np
import rioxarray  # ruff:ignore[unused-import]
import xarray as xr
from cyclopts import Parameter

logger = logging.getLogger(__name__)


def path_or_web_url_validator(type_: object, path_or_url: Path | str) -> str:  # ruff:ignore[unused-function-argument] - needed type_ to match cyclopts validator signature
    """Validate that the provided path is either a valid local path or a web URL.

    NOTE: Path("") resolves to the current directory (".") - which is considered a valid path.

    Parameters
    ----------
    type_ : object
        The type of the argument being validated - unused, purely required for cyclopts validator compatibility.
    path_or_url : Path | str
        The candidate path or URL to validate.

    Returns
    -------
    str
        The validated path or URL as a string.

    Raises
    ------
    ValueError
        If the provided path is neither a valid local path nor a web URL.
    """
    if not Path(path_or_url).exists() and not str(path_or_url).startswith(("http://", "https://")):
        raise ValueError(f"The provided path '{path_or_url}' is neither a valid local path nor a web URL.")
    return str(path_or_url)


# Alias for use in cyclopts Parameter definitions.
# Override the sphinx documentation entry in the lines below since cyclopts doesn't produce nice results for
# custom validators:
#: A local path or web URL, represented as a string.
#:
#: Validated by :func:`path_or_web_url_validator`.
#:
#: :meta hide-value:
type ZarrType = Annotated[str, Parameter(validator=path_or_web_url_validator)]


def save_geo_zarr(ds: xr.Dataset, zarr_path: Path | str, target_chunk_size: str = "1024MB") -> None:
    """Save dataset to GeoZarr adhering to CF conventions (https://cfconventions.org/).

    Parameters
    ----------
    ds: xr.Dataset
        The xarray dataset to save. This must follow the GeoZarr format,
        including having a `spatial_ref` attribute with a valid CRS.
    zarr_path: Path | str
        Full path of the output .zarr file to be created.
    target_chunk_size: str
        The target chunk size for the GeoZarr dataset, by default "1024MB".
    """
    ds = ds.rio.write_crs(ds["spatial_ref"].attrs["crs_wkt"], grid_mapping_name="spatial_ref")
    ds = ds.reset_coords("spatial_ref")
    if "time" in ds:
        ds["time"] = prepare_cf_time_data_array(ds["time"])
    # Prepare the spatial variables with the required attributes for GeoZarr.
    for name, axis in (("x", "X"), ("y", "Y")):
        if name in ds:
            ds[name] = prepare_cf_spatial_data_array(ds[name], name=name, axis=axis)
    ds = prepare_cf_data_variable_arrays(ds)

    ds.attrs = {"_ARRAY_DIMENSIONS": ["spatial_ref", *sorted(ds.dims)]}
    ds = ds.chunk(chunks=target_chunk_size)

    ds.to_zarr(
        str(zarr_path),
        mode="w",
        zarr_format=2,
        consolidated=True,
    )


def prepare_cf_data_variable_arrays(ds: xr.Dataset) -> xr.Dataset:
    """Assign CF grid-mapping and fill-value metadata to dataset variables."""
    for name in ds.variables:
        variable = ds[name]
        # Do not reuse source packing, time units, compression, or chunk encodings.
        variable.encoding = {"coordinates": None}
        if name != "spatial_ref":
            variable.attrs["grid_mapping"] = "spatial_ref"

        dtype = variable.dtype
        fill_value = None
        if variable.ndim:
            if np.issubdtype(dtype, np.unsignedinteger):
                fill_value = np.iinfo(dtype).max
            else:
                fill_value = {"f": -9999.0, "i": -9999, "b": 0, "U": ""}.get(dtype.kind)
        # We want to preserve NaNs so set the _FillValue in .attrs, rather than .encoding.
        variable.attrs["_FillValue"] = fill_value
    return ds


def prepare_cf_spatial_data_array(
    spatial_data_array: xr.DataArray,
    name: str,
    axis: str,
) -> xr.DataArray:
    """Assign the required GeoZarr attributes to a spatial coordinate."""
    return spatial_data_array.assign_attrs(
        variable_name_override=name,
        standard_name=f"projection_{name}_coordinate",
        long_name=f"projection_{name}_coordinate",
        units="m",
        axis=axis,
    )


def prepare_cf_time_data_array(time_data_array: xr.DataArray) -> xr.DataArray:
    """Convert a time coordinate to numeric epoch seconds with legacy CF metadata."""
    logger.info("Converting type %s to seconds since epoch", time_data_array.dtype)
    time_data_array = time_data_array.astype("datetime64[s]").astype("float64")
    return time_data_array.assign_attrs(
        {
            "standard_name": "time",
            "long_name": "time",
            "axis": "T",
            "bounds": "time_bnds",
            "units": "seconds since 1970-01-01T00:00:00",
            "calendar": "standard",
        }
    )
