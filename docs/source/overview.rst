Overview
========

The Functional Time Series (FTS) module compresses a panel of space-time
ice-sheet data into a few coefficients per cycle. Much of what an ice sheet
does is seasonal — melt begins, peaks and ends; lakes fill and drain; runoff
rises and falls — and a list of monthly values describes that clumsily,
treating each month as an independent number when what distinguishes one season
from another is largely the shape of the curve those months trace out.

The module cuts each point's series into curves of a fixed length and reduces
the whole panel to a small matrix per variable, which can then be used as a
predictor alongside other quantities.

Scientific background
---------------------

The reduction follows ``hdfpca`` from the R package ``hdftsa``, after Gao,
Shang and Yang (2019), *High-dimensional functional time series forecasting: an
application to age-specific mortality rates*, Journal of Multivariate Analysis
170, 232–243. It has two stages.

The first is per location. Each point's sequence of curves is reduced by
**dynamic functional principal component analysis**. The word *dynamic* is the
important one: rather than the ordinary covariance between curves, the
eigenvectors are taken of a **long-run covariance**, which sums autocovariances
at lags one to five under a flat-top weighting. What this picks out is the
structure carried across successive curves, not simply the directions of
greatest variance within them.

The second stage works across space. The leading score from every point is
stacked side by side and reduced again, using the eigenvectors of a sum of
lagged autocovariance products across the whole panel. Serial dependence drives
this stage too; it is not a purely spatial decomposition. The module keeps the
scores from this second stage and stops there — the rest of ``hdfpca``, and the
forecasting the R package builds on top of it, are not ported.

The resulting **coefficients are weights on shared patterns, not physical
quantities**. A curve with a large leading coefficient resembles the dominant
pattern across the record more strongly than one with a small coefficient. The
sign carries no meaning: each column is flipped so that its largest-magnitude
entry is positive, a convention imposed for reproducibility. The eigenvectors
are not returned, so the input cannot be reconstructed from the output.

Repository structure
--------------------

``fts_spatial_entrypoint.py``
    Command-line entry point. Takes an input path and an output path, with the
    curve length and the number of components and factors to retain.

``fts_spatial_module/fts_spatial_module.py``
    Selects the variables to analyse, flattens gridded input to a point cloud,
    cuts the series into whole curves, and runs the reduction.

``fts_spatial_module/hdfpca.py``
    The two-stage reduction itself, ported to NumPy from ``hdfpca.R``,
    ``dfpca.R``, ``hdpca.R`` and ``wtflat.R`` in ``hdftsa`` 1.1. The reduction
    was first applied in this project through an R implementation by Israel
    Martinez Hernandez, which called ``hdfpca`` directly; the port replaces
    that dependency on R, reproduces its edge-case conventions deliberately,
    and is checked against stored reference output.

Data flow
---------

1. Read the input and select the variables defined over time and space.
   Gridded input is flattened to a point cloud, keeping only points where every
   variable is finite at every timestep, so that a curve can be built without
   gaps and the coefficients stay comparable between variables. Input that
   already carries a point dimension is taken as given and is not filtered.
2. Cut each point's series into whole curves, discarding trailing timesteps
   that do not complete one. A variable empty across a whole timestep is
   dropped rather than being allowed to empty the shared mask.
3. Run the two-stage reduction and write one matrix per variable, curves by
   factors.

The number of factors to retain is set by the caller. The number of principal
components is accepted but does not change the result, because only the leading
component of each point feeds the second stage. Since that stage sums
autocovariances over lags one to five, it needs at least six curves, and fewer
are refused rather than silently fitted.

Use within DTC-IS
-----------------

The module is deployed in the DTC-IS environment and runs as a job through the
`Query API <https://query.dtc-ice-sheets.org>`_, usually on the output of the
Gaussian Random Field module. In Use Case 1 it compresses each melt season of a
supraglacial lake record into coefficients that serve as a predictor of river
discharge in the Watson catchment.

The module reduces; it does not predict. Relating the coefficients to climate
drivers, or to something measured elsewhere, belongs to whatever consumes the
output.
