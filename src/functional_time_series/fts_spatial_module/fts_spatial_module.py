"""Functional Time Series (FTS) module."""

import logging

import numpy as np
import xarray as xr

from functional_time_series.fts_spatial_module.hdfpca import hdfpca_factor_scores

logger = logging.getLogger(__name__)


class FtsSpatialModule:
    """FTS module: compresses each point's time series into curve coefficients."""

    @staticmethod
    def spatial_variables(ds: xr.Dataset) -> list[str]:
        """
        Names of the variables to analyse: those with exactly the dimensions time and point_id.

        Parameters
        ----------
        ds : xr.Dataset
            The input dataset.

        Returns
        -------
        list[str]
            Names of the variables with dimensions ("time", "point_id"), in dataset order.
        """
        return [str(v) for v in ds.data_vars if set(ds[v].dims) == {"time", "point_id"}]

    @staticmethod
    def as_point_cloud(ds: xr.Dataset) -> xr.Dataset:
        """
        Flatten a gridded dataset into the point cloud form this module analyses.

        The GRF module interpolates its fit onto a regular grid before saving, so its published output has
        the dimensions ("time", "x", "y") rather than ("time", "point_id"). Without this step
        :func:`spatial_variables` matches nothing, every variable is dropped, and the result is an empty
        dataset with no warning. Cells outside the fitted area are NaN and a curve cannot be built across a
        gap, so a point is kept only where every variable is finite at every timestep.

        Parameters
        ----------
        ds : xr.Dataset
            The input dataset, either already a point cloud carrying a "point_id" dimension, or a grid
            carrying x and y dimensions under those names or the old
            "projection_x_coordinate" and "projection_y_coordinate".

        Returns
        -------
        xr.Dataset
            The same data with the two spatial dimensions stacked into "point_id", any variable that is
            empty at a whole timestep removed, and the points that are not finite throughout the record
            removed. A dataset that is already a point cloud, or that is neither form, is returned
            unchanged.
        """
        if FtsSpatialModule.spatial_variables(ds):
            return ds
        if "projection_x_coordinate" in ds.dims:
            ds = ds.rename({"projection_x_coordinate": "x", "projection_y_coordinate": "y"})
        if "x" not in ds.dims or "y" not in ds.dims:
            return ds

        ds = ds.stack({"point_id": ["x", "y"]})
        names = FtsSpatialModule.spatial_variables(ds)
        if not names:
            return ds

        # A variable can be empty at a whole timestep while the others are not. The GRF
        # module fills a snapshot that has no spatial variation with the constant it
        # already holds and writes NaN for its uncertainty, since no posterior was
        # sampled - so the uncertainty field of a lake record arrives empty for the
        # months when the lakes have drained. There is no curve to fit through a
        # timestep with no values, and the mask below is shared, so leaving such a
        # variable in would drop every point in the catchment rather than analysing
        # what is there.
        empty = [v for v in names if not np.isfinite(ds[v].values).any(axis=1).all()]
        if empty:
            logger.warning("not analysing %s: no values at all in at least one timestep", ", ".join(empty))
            ds = ds.drop_vars(empty)
            names = FtsSpatialModule.spatial_variables(ds)
        if not names:
            return ds

        finite = np.all([np.isfinite(ds[v].values) for v in names], axis=(0, 1))
        if not finite.all():
            logger.warning(
                "dropping %d of %d points that are not finite throughout the record", int((~finite).sum()), finite.size
            )
        return ds.isel(point_id=finite)

    @staticmethod
    def whole_curves(data: np.ndarray, curve_time_points: int, name: str = "data") -> np.ndarray:
        """
        Truncate the time axis to a whole number of curves.

        Every curve needs exactly curve_time_points values, so any time steps left over after the last
        whole curve are dropped, with a warning.

        This matches hdftsa as called from R. There, ``matrix(Coeff[k, ], delta, nn)`` with a non-whole
        ``nn = tt/delta`` rounds the number of curves down and drops the same trailing values, so the
        result is identical. The differences are that R does it silently apart from one R warning per
        point, while this logs a single warning per variable.

        Parameters
        ----------
        data : np.ndarray
            Array of shape (n_points, n_times).
        curve_time_points : int
            Number of time points in one curve.
        name : str, optional
            Variable name, used only in the warning.

        Returns
        -------
        np.ndarray
            The input, truncated to the first whole number of curves along the time axis.
        """
        n_times = data.shape[1]
        n_whole = (n_times // curve_time_points) * curve_time_points
        if n_whole != n_times:
            logger.warning(
                "%s: %d time points is not a multiple of %d; using the first %d",
                name,
                n_times,
                curve_time_points,
                n_whole,
            )
        return data[:, :n_whole]

    def run_spatial_analysis(
        self,
        ds: xr.Dataset,
        curve_time_points: int = 24,
        num_pca: int = 2,
        num_factors: int = 3,
    ) -> xr.Dataset:
        """
        Run spatial analysis on the given dataset.

        Parameters
        ----------
        ds : xr.Dataset
            The input dataset for analysis.
            This must have the following dimensions:

            - time
            - point_id (flattened spatial coordinates)

            Each of the spatial variables with dimensions ("time", "point_id") will be analyzed.
            A gridded dataset with x and y dimensions in place of point_id is accepted too, and is
            flattened by :func:`as_point_cloud` first.
        curve_time_points : int, optional
            Number of time points to be considered as a single curve
        num_pca : int, optional
            Number of principal components to retain. Only the first feeds the second step, so this does
            not change the result.
        num_factors: int, optional
            Number of factors to retain in the second step dimension reduction

        Returns
        -------
        xr.Dataset
            The output dataset containing the model results, with dimensions:

            - curve_id (length equal to the number of time points divided by curve_time_points, rounded down)
            - coef_id (length equal to num_factors)

            Each of the input spatial variables will have a corresponding variable of the same name, with dimensions
            ("curve_id", "coef_id").
            The input spatial_ref variable and its CRS metadata are preserved as a coordinate, if present.
            Any other input variables will be discarded.

        Raises
        ------
        ValueError
            If the dataset holds no variable this module can analyse. Returning an empty dataset instead
            would record the run as a success and leave the caller with a store holding no data.
        """
        ds = self.as_point_cloud(ds)
        names = self.spatial_variables(ds)
        if not names:
            raise ValueError(
                "no variable has the dimensions time and point_id, and the dataset is not a grid that can "
                f"be flattened into them, so there is nothing to analyse: found {list(ds.data_vars)} over "
                f"dimensions {list(ds.dims)}",
            )

        data_vars = {}
        for var in names:
            data = self.whole_curves(ds[var].transpose("point_id", "time").values, curve_time_points, var)
            result = hdfpca_factor_scores(data, curve_time_points, num_pca, num_factors)
            data_vars[var] = (("curve_id", "coef_id"), result)

        output = xr.Dataset(
            data_vars=data_vars,
            coords={
                "curve_id": np.arange(ds.sizes["time"] // curve_time_points),
                "coef_id": np.arange(num_factors),
            },
        )
        if "spatial_ref" in ds:
            output = output.assign_coords(spatial_ref=ds["spatial_ref"].copy(deep=True))
        return output
