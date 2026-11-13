.. DTC Python Template documentation master file, created by
   sphinx-quickstart on Tue Jan 21 10:28:47 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Functional Time Series Module
=============================

The ``functional_time_series`` module is a sub-component of the
`Digital Twin Component for Ice Sheets (DTC-IS) <https://dtc-ice-sheets.org/>`_, an ESA-funded, EO-driven digital twin
project.

The ``functional_time_series`` module reduces a panel of space-time ice-sheet data to a few coefficients per
cycle. Each point's series is cut into curves - melt seasons, typically - and the whole panel is compressed, so
that a season at every location becomes a handful of numbers usable as a predictor alongside other variables.
It usually runs on the output of the
`gaussian_random_field <https://github.com/DTC-Ice-Sheets/gaussian_random_field>`_ module.

It was developed to provide hydropower insights for the Watson River catchment on the Greenland ice sheet,
where it compresses each melt season of a supraglacial lake record into a predictor of river discharge.

For a more interactive way of exploring the system, please visit the digital twin's
`dashboard <https://dashboards.dtc-ice-sheets.org/>`_ and/or explore the `notebooks <https://colab.research.google.com/github/DTC-Ice-Sheets/dtc_is_notebooks/>`_.

The remainder of this documentation page provides an overview of the ``functional_time_series`` module and details
on how to get it running locally.

.. toctree::
   :maxdepth: 1

   Overview <overview>
   Quick Start <quick-start>
   API Reference <api>
