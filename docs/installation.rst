.. _installing:

##########################
Installing
##########################

Dependencies:
----------------------------------------
- Python 2.6 or later
- netCDF4_ : python/numpy interface for NetCDF4

analysis module dependencies:
----------------------------------------
  - numpy_ : Fundamental package needed for scientific computing with Python.
  - scipy_: Open-source software for mathematics, science, and engineering.
  - pandas_ : Flexible and powerful data analysis / manipulation library for Python.
  - xarray_ :  N-D labeled arrays and datasets in Python.
  - matplotlib basemap_ : The matplotlib basemap toolkit is a library for plotting 2D data on maps in Python.
  - seaborn_ : Statistical data visualization using matplotlib.

post_process module dependencies:
----------------------------------------
  - nco_ :  Command line utility that manipulates and analyzes data stored in netCDF-accessible formats, including DAP, HDF4, and HDF5.
  - `nco-bindings`__ : Python bindings to nco comamnd line utility.

Before you install `rasmlib`, be sure you have the required dependencies (python and netCDF4) installed. The easiest way to do so is to use the Anaconda_ python distribution.

To install `rasmlib`, use `setuptools`:

    python setup.py install


.. _netCDF4: https://github.com/Unidata/netcdf4-python
.. _numpy: http://www.numpy.org/
.. _scipy: http://docs.scipy.org/doc/
.. _pandas: http://pandas.pydata.org/
.. _xarray: http://xray.readthedocs.org/en/stable/index.html
.. _basemap: http://matplotlib.org/basemap/index.html
.. _seaborn: http://web.stanford.edu/~mwaskom/software/seaborn/index.html
.. _nco: http://nco.sourceforge.net/
.. _ncopy: https://github.com/jhamman/nco-bindings
.. _anaconda: https://store.continuum.io/cshop/anaconda/
__ ncopy_
