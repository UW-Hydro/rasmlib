.. _post_processing:

###########
Post Processing
###########

After each CESM / RASM simulation, it is necessary to do a certain amount of post processing on the output files before comprehensive analysis can be performed.  This page provides a brief explination of how a few of the RASM model components are post processed and the tools used from the `rasmlib` package.

Tools
----------------------------------------
`rasm_post_process`:  a command line script that performs batch post processing on model output.  It can be run using multiple cores using Python `multiprocessing` module.  See the `Scripts` page for information on usage of this script.

Models
----------------------------------------

**VIC**

Current VIC output requires two steps of post processing.
Step 1:  The time variable and filename timestamp must be corrected.
Step 2:  Model output is concatenated to 3 levels of model output:
    - monthly mean diurnal cycles: requires hourly model output timestep.
    - daily timeseries:  requires daily or subdaily model output timestep.
    - monthly timeseries: performed for any model output timestep

**RVIC, CPL, WRF**

Output from these models is post processed in the same way as for the **VIC** model except there is no need to adjust the time variables or timesteps.
