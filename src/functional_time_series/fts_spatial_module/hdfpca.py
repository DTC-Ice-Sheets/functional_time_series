"""
NumPy port of the part of ``hdfpca`` from the R package hdftsa that the FTS module uses.

``hdfpca`` (Gao, Shang and Yang, 2019) reduces a panel of functional time series in two steps. First, each point's
series of curves is reduced by dynamic functional principal component analysis (``dfpca``), which takes the
eigenvectors of a long-run covariance rather than of the plain covariance. Second, the scores of the first principal
component from every point are stacked side by side and reduced again (``hdpca``), this time using the eigenvectors of
a sum of lagged autocovariance products. The FTS module keeps only the scores from that second step, so this port
stops there and leaves out the rest of ``hdfpca`` and the forecast that ``forecast.hdfpca`` builds on top of it.

The sources ported are ``hdfpca.R``, ``dfpca.R``, ``hdpca.R`` and ``wtflat.R`` from hdftsa 1.1.
"""

import numpy as np

# hdpca's h0: the number of lags summed in the second step.
LONG_RUN_LAGS = 5


def flat_top_weight(x: float, c: float = 0.5) -> float:
    """
    Flat-top kernel weight, ported from hdftsa's ``wtflat``.

    The weight is 1 on (-c, c) and falls linearly to 0 at -1 and 1. The R function uses a strict ``c < x`` on its
    rising edge, so ``x == c`` falls through to 0 while ``x == -c`` gives 1. The quirk is kept so that results match
    R exactly; it only matters when the number of curves is four times a square.

    Parameters
    ----------
    x : float
        Point at which to evaluate the kernel.
    c : float
        Half-width of the flat top.

    Returns
    -------
    float
        The kernel weight at `x`.
    """
    if -1 < x <= -c:
        return (x + 1) / (1 - c)
    if -c < x < c:
        return 1.0
    if c < x < 1:
        return (x - 1) / (c - 1)
    return 0.0


def descending_eigenvectors(matrix: np.ndarray, order: int) -> np.ndarray:
    """
    Eigenvectors of the leading `order` eigenvalues of one or more symmetric matrices, largest first.

    Like R's ``eigen(symmetric = TRUE)``, only the lower triangle is read.

    Parameters
    ----------
    matrix : np.ndarray
        Array of shape (..., n, n).
    order : int
        Number of eigenvectors to keep.

    Returns
    -------
    np.ndarray
        Array of shape (..., n, order), one eigenvector per column.
    """
    _, vectors = np.linalg.eigh(matrix, UPLO="L")
    return vectors[..., ::-1][..., :order]


def dynamic_fpca_scores(curves: np.ndarray, order: int) -> np.ndarray:
    """
    Principal component scores of each point's curves, ported from hdftsa's ``dfpca``.

    R loops over every lag and multiplies its autocovariance by the kernel weight, so it builds lags whose weight is
    zero. Those lags are skipped here, which leaves roughly the square root of the number of curves to compute.

    Parameters
    ----------
    curves : np.ndarray
        Array of shape (n_points, n_curves, curve_time_points).
    order : int
        Number of principal components to keep.

    Returns
    -------
    np.ndarray
        Array of shape (n_points, n_curves, order).
    """
    n_curves = curves.shape[1]
    bandwidth = np.sqrt(n_curves)
    centred = curves - curves.mean(axis=1, keepdims=True)

    long_run = np.einsum("pti,ptj->pij", centred, centred) / n_curves
    for lag in range(1, n_curves):
        # R stores lag `lag` at index h = lag + 1 and weights it by h / bandwidth, not lag / bandwidth.
        h = lag + 1
        forward = flat_top_weight(h / bandwidth)
        backward = flat_top_weight(-h / bandwidth)
        if forward == 0 and backward == 0:
            continue
        autocovariance = np.einsum("pti,ptj->pij", centred[:, lag:], centred[:, : n_curves - lag]) / n_curves
        long_run += forward * autocovariance + backward * np.swapaxes(autocovariance, 1, 2)

    return centred @ descending_eigenvectors(long_run, order)


def column_space_factor_directions(centred: np.ndarray, order: int) -> np.ndarray:
    """
    Leading eigenvectors of hdpca's ``f``, built directly as R builds it.

    ``f`` is the sum over lags 1 to 5 of ``C_h C_hᵀ``, where ``C_h`` is the lag-h autocovariance between columns, so
    it is n_columns square. That is cheap when there are few columns.

    Parameters
    ----------
    centred : np.ndarray
        Array of shape (n_rows, n_columns), with each column's mean removed.
    order : int
        Number of eigenvectors to keep.

    Returns
    -------
    np.ndarray
        Array of shape (n_columns, order), one eigenvector per column.
    """
    n_rows = centred.shape[0]
    f = np.zeros((centred.shape[1], centred.shape[1]))
    for h in range(1, LONG_RUN_LAGS + 1):
        lagged = centred[h:].T @ centred[: n_rows - h] / (n_rows - h)
        f += lagged @ lagged.T
    return descending_eigenvectors(f, order)


def row_space_factor_directions(centred: np.ndarray, order: int) -> np.ndarray:
    """
    Leading eigenvectors of hdpca's ``f``, found without building ``f``.

    Each ``C_h C_hᵀ`` in ``f`` equals ``centredᵀ M_h centred`` for a matrix ``M_h`` that is only n_rows square, so
    ``f`` is ``Zᵀ Z`` for a ``Z`` with n_rows rows, and its eigenvectors are the right singular vectors of ``Z``. That
    is cheap when there are few rows, which for a grid of points is the usual case.

    Parameters
    ----------
    centred : np.ndarray
        Array of shape (n_rows, n_columns), with each column's mean removed.
    order : int
        Number of eigenvectors to keep.

    Returns
    -------
    np.ndarray
        Array of shape (n_columns, order), one eigenvector per column.
    """
    n_rows = centred.shape[0]
    weight = np.zeros((n_rows, n_rows))
    for h in range(1, LONG_RUN_LAGS + 1):
        earlier = centred[: n_rows - h]
        weight[h:, h:] += earlier @ earlier.T / (n_rows - h) ** 2

    values, vectors = np.linalg.eigh(weight)
    z = (vectors * np.sqrt(np.clip(values, 0, None))).T @ centred
    _, _, right = np.linalg.svd(z, full_matrices=False)
    return right[:order].T


def factor_scores(x: np.ndarray, order: int) -> np.ndarray:
    """
    Scores of the leading factors of a panel of series, ported from the ``coef`` that hdftsa's ``hdpca`` returns.

    R takes the eigenvectors of ``f``, a matrix that is n_columns square. A column here is a spatial point, so for a
    grid ``f`` can have tens of thousands of rows, while a row is a curve, of which there are usually few. The
    eigenvectors are therefore found from whichever of the columns or the rows is smaller; both give the same result.

    An eigenvector's sign is arbitrary, and R and NumPy may pick differently, so each column of scores is flipped to
    make its largest-magnitude entry positive. R's result therefore agrees with this one up to the sign of each column.

    Parameters
    ----------
    x : np.ndarray
        Array of shape (n_rows, n_columns): one row per time step, one column per series.
    order : int
        Number of factors to keep.

    Returns
    -------
    np.ndarray
        Array of shape (n_rows, order).

    Raises
    ------
    ValueError
        If there are too few rows for the lags summed, or more factors are asked for than the data can supply.
        hdftsa stops with "subscript out of bounds" for too few rows, or for more factors than columns. Given more
        factors than rows but no more than columns, it returns scores along eigenvectors whose eigenvalue is zero,
        which are arbitrary; this refuses that case too.
    """
    n_rows, n_columns = x.shape
    if n_rows <= LONG_RUN_LAGS:
        raise ValueError(
            f"{n_rows} curves is too few: the second step sums autocovariances over lags 1 to {LONG_RUN_LAGS}, so it "
            f"needs at least {LONG_RUN_LAGS + 1}"
        )
    if order > min(n_rows, n_columns):
        raise ValueError(
            f"cannot keep {order} factors from {n_rows} curves at {n_columns} points; at most "
            f"{min(n_rows, n_columns)} are available"
        )

    centred = x - x.mean(axis=0)
    if n_columns <= n_rows:
        scores = centred @ column_space_factor_directions(centred, order)
    else:
        scores = centred @ row_space_factor_directions(centred, order)

    largest = scores[np.abs(scores).argmax(axis=0), np.arange(order)]
    return scores * np.where(largest < 0, -1.0, 1.0)


def hdfpca_factor_scores(
    data: np.ndarray,
    curve_time_points: int,
    num_pca: int,
    num_factors: int,
) -> np.ndarray:
    """
    Compute the factor scores of the first principal component, as ``hdfpca(...)$model$model2[[1]]$coef`` in hdftsa.

    Only the first principal component of each point feeds the second step, so `num_pca` does not change the result
    beyond needing to be at least 1. The data is converted to float64 first, as R works in double precision whatever
    the input.

    Parameters
    ----------
    data : np.ndarray
        Array of shape (n_points, n_times), with n_times a whole multiple of `curve_time_points`.
    curve_time_points : int
        Number of time points in one curve.
    num_pca : int
        Number of principal components to keep in the first step.
    num_factors : int
        Number of factors to keep in the second step.

    Returns
    -------
    np.ndarray
        Array of shape (n_curves, num_factors).
    """
    n_points, n_times = data.shape
    curves = np.asarray(data, dtype=np.float64).reshape(n_points, n_times // curve_time_points, curve_time_points)
    first_component = dynamic_fpca_scores(curves, num_pca)[:, :, 0]
    return factor_scores(first_component.T, num_factors)
