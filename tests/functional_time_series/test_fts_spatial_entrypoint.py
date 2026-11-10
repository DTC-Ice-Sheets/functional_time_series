"""Testing the FTS spatial entrypoint."""

import importlib
import re
import sys
from pathlib import Path
from types import ModuleType
from unittest.mock import MagicMock, patch

import cyclopts
import pytest
import xarray as xr

# The entrypoint's writer import belongs to the GRF package it was copied from.
grf_module = ModuleType("functional_time_series.grf_module")
grf_module.__path__ = []
geo_zarr_utils = ModuleType("functional_time_series.grf_module.geo_zarr_utils")
geo_zarr_utils.save_geo_zarr = MagicMock()
with patch.dict(
    sys.modules,
    {
        "functional_time_series.grf_module": grf_module,
        "functional_time_series.grf_module.geo_zarr_utils": geo_zarr_utils,
    },
):
    entrypoint = importlib.import_module("functional_time_series.fts_spatial_entrypoint")

app = entrypoint.app


def _invoke_app(
    input_path: str,
    output_path: Path,
    options: list[str] | None = None,
) -> None:
    argv = [
        "fts_spatial_entrypoint",
        "--data_path",
        input_path,
        "--output_path",
        str(output_path),
        *(options or []),
    ]
    with patch.object(sys, "argv", argv):
        app(result_action="return_value", exit_on_error=False)


def test_main_loads_analyzes_and_saves_with_default_parameters(
    test_dataset_1_path: Path, tmp_zarr_output_path: Path
) -> None:
    loaded_dataset = xr.Dataset()
    analyzed_dataset = xr.Dataset()

    with (
        patch.object(entrypoint, "setup_logging") as setup_logging,
        patch.object(entrypoint.xr, "load_dataset", return_value=loaded_dataset) as load_dataset,
        patch.object(entrypoint, "FtsSpatialModule") as fts_module_class,
        patch.object(entrypoint, "save_geo_zarr") as save_geo_zarr,
    ):
        fts_module_class.return_value.run_spatial_analysis.return_value = analyzed_dataset

        _invoke_app(str(test_dataset_1_path), tmp_zarr_output_path)

    setup_logging.assert_called_once_with()
    load_dataset.assert_called_once_with(str(test_dataset_1_path), engine="zarr", chunks={})
    fts_module_class.return_value.run_spatial_analysis.assert_called_once_with(
        loaded_dataset,
        curve_time_points=5,
        num_pca=2,
        num_factors=3,
    )
    save_geo_zarr.assert_called_once_with(analyzed_dataset, tmp_zarr_output_path)


def test_main_forwards_custom_analysis_parameters(
    test_dataset_1_path: Path, tmp_zarr_output_path: Path
) -> None:
    loaded_dataset = xr.Dataset()
    with (
        patch.object(entrypoint.xr, "load_dataset", return_value=loaded_dataset),
        patch.object(entrypoint, "FtsSpatialModule") as fts_module_class,
        patch.object(entrypoint, "save_geo_zarr"),
    ):
        _invoke_app(
            str(test_dataset_1_path),
            tmp_zarr_output_path,
            ["--curve_time_points", "8", "--num_pca", "4", "--num_factors", "6"],
        )

    fts_module_class.return_value.run_spatial_analysis.assert_called_once_with(
        loaded_dataset,
        curve_time_points=8,
        num_pca=4,
        num_factors=6,
    )


def test_main_accepts_a_web_url(tmp_zarr_output_path: Path) -> None:
    input_url = "https://example.com/input.zarr"
    loaded_dataset = xr.Dataset()
    with (
        patch.object(entrypoint.xr, "load_dataset", return_value=loaded_dataset) as load_dataset,
        patch.object(entrypoint, "FtsSpatialModule") as fts_module_class,
        patch.object(entrypoint, "save_geo_zarr"),
    ):
        _invoke_app(input_url, tmp_zarr_output_path)

    load_dataset.assert_called_once_with(input_url, engine="zarr", chunks={})
    fts_module_class.return_value.run_spatial_analysis.assert_called_once()


@pytest.mark.parametrize(
    "test_source",
    ["invalid_file_path", "invalid_path", "invalid_url"],
)
def test_main_rejects_invalid_input_paths(
    test_source: str, tmp_zarr_output_path: Path, request: pytest.FixtureRequest
) -> None:
    invalid_path = request.getfixturevalue(test_source)
    with (
        patch.object(entrypoint.xr, "load_dataset") as load_dataset,
        patch.object(entrypoint, "FtsSpatialModule") as fts_module_class,
        pytest.raises(
            cyclopts.ValidationError,
            match=re.escape(f"The provided path '{invalid_path}' is neither a valid local path nor a web URL."),
        ),
    ):
        _invoke_app(str(invalid_path), tmp_zarr_output_path)

    load_dataset.assert_not_called()
    fts_module_class.assert_not_called()
