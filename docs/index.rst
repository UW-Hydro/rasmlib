Welcome to the rasmlib documentation
====================================

RASM
----------------------------------------

The `Regional Arctic System Model (RASM)`__ is a high resolution, regional, coupled atmosphere - land - sea ice - ocean model that uses the `Community Earth System Model (CESM)`__ coupling infrastructure (CPL7) over a Pan-Arctic domain. RASM is composed of five model components:
  - `Weather Research and Forecasting (WRF) atmospheric model`__
  - `Variable Infiltration Capacity (VIC) hydrology model`__
  - `RVIC streamflow routing model`__
  - `Parallel Ocean Program (POP) model`__
  - `Los Alamos Sea Ice model (CICE)`__

.. _RASM: http://www.oc.nps.edu/NAME/RASM.htm
.. _CESM: http://www2.cesm.ucar.edu/
.. _WRF: http://www.wrf-model.org/index.php
.. _VIC: http://www.hydro.washington.edu/Lettenmaier/Models/VIC/index.shtml
.. _RVIC: https://github.com/jhamman/RVIC
.. _POP: http://www.cesm.ucar.edu/models/cesm1.0/pop2/
.. _CICE: http://oceans11.lanl.gov/trac/CICE
__ RASM_
__ CESM_
__ WRF_
__ VIC_
__ RVIC_
__ POP_
__ CICE_

Purpose of rasmlib
----------------------------------------
`rasmlib` is a collection of post-processing, analysis, and plotting tools used to support the development of RASM.  Many of the scripts and functions in this package were developed specifically for RASM and are likely not applicable to other projects.

Documentation
=================

.. toctree::
   :maxdepth: 1

   installation
   post_processing
   analysis
   api
   api_post
   api_analysis
   scripts

Modules
----------------------------------------
`rasmlib` is organized into two modules:

- `analysis`
- `post_processing`

Command Line Utilities:
----------------------------------------
`rasmlib` includes two main command line utilities:

- `rasm_analysis`
- `rasm_post_process`

License
----------------------------------------
`rasmlib` is available under the GNU GPL v3.0 license.

