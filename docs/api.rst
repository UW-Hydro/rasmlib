
.. _api:

.. automodule:: rasmlib

#####################
API
#####################

This page provides an auto-generated summary of rasmlib's API.

Modules
-----------------------

.. autosummary::
   :toctree: generated/

      analysis
      post_processing

calendar
-----------------------
.. currentmodule:: rasmlib.calendar

Functions:

.. autosummary::
   :toctree: generated/

      next_month
      prev_month
      next_day
      prev_day
      leap_year
      get_dpm

Attributes:

.. autosummary::
   :toctree: generated/

      HOURSPERDAY
      SECSPERHOUR
      MINSPERHOUR
      SECSPERMINUTE
      MINSPERDAY
      SECSPERDAY
      MONTHSPERYEAR
      dpm
      seasons

io
-----------------------
.. currentmodule:: rasmlib.io

Functions:

.. autosummary::
   :toctree: generated/

      read_config
      config_type
      make_tarfile
      get_data_files_namelist
      get_variables_namelist
      get_datasets
      read_domain
      get_time_units

utils
-----------------------
.. currentmodule:: rasmlib.utils

Functions:

.. autosummary::
   :toctree: generated/

      argsort
      clean_dir
      clean_file
      make_directories
      multiple_replace
      custom_strptime
      custom_strftime
      partition
      chunks
