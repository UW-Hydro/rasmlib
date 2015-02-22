rasmlib:  RASM analysis and post-processing tools
=========

[![Build Status](https://travis-ci.org/jhamman/rasmlib.svg)](https://travis-ci.org/jhamman/rasmlib) [![Documentation Status](https://readthedocs.org/projects/rasmlib/badge/?version=develop)](https://readthedocs.org/projects/rasmlib/?badge=develop)

`rasmlib` documentation: http://rasmlib.readthedocs.org/en/develop/

## RASM
The [Regional Arctic System Model (RASM)](http://www.oc.nps.edu/NAME/RASM.htm) is a high resolution, regional, coupled atmosphere - land - sea ice - ocean model that uses the [Community Earth System Model (CESM)](http://www2.cesm.ucar.edu/) coupling infrastructure (CPL7) over a Pan-Arctic domain. RASM is composed of five model components:

  - [Weather Research and Forecasting (WRF) atmospheric model](http://www.wrf-model.org/index.php)
  - [Variable Infiltration Capacity (VIC) hydrology model](http://www.hydro.washington.edu/Lettenmaier/Models/VIC/index.shtml)
  - [RVIC streamflow routing model](https://github.com/jhamman/RVIC)
  - [Parallel Ocean Program (POP) model](http://www.cesm.ucar.edu/models/cesm1.0/pop2/)
  - [Los Alamos Sea Ice model (CICE)](http://oceans11.lanl.gov/trac/CICE)

## Purpose of rasmlib
`rasmlib` is a collection of post-processing, analysis, and plotting tools used to support the development of RASM.  Many of the scripts and functions in this package were developed specifically for RASM and are likely not applicable to other projects.

## History
`rasmlib` development began in 2012 to support the development of RASM.  Joe Hamman has done most of the development with support and contributions from of [Bart Nijssen](@bartnijssen) and [Mu Xiao](@muxiao).

## License
rasmlib analysis and post processing toolbox. Copyright (C) 2015, Joe Hamman

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

rasmlib includes portions of pandas, matplotlib, xray, netCDF, and other projects. See those projects for licenses specific use policies related to their licesnses.
