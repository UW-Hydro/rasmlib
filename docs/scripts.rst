.. _scripts:

###########
Scripts
###########

`rasmlib` includes two command line scripts:

1.  **rasm_post_process**:  this script runs a set of utilities to concatenate and adjust RASM output files.  Its implementation is fairly generaric and has been used to post process output from the VIC, RVIC, ATM, and CPL components.
2.  **rasm_analysis**:  this script is currently underdevelopment and is intended to provide a very flexible interface to doing batch analysis on RASM output.

-------------------------------------------------------------------------------

rasm_post_process
-----------------------
usage: rasm_post_process [-h] [--output_preset OUTPUT_PRESET]
                         [--processed_dir PROCESSED_DIR] [-np NUMOFPROC]
                         config_file short_term_archive

Generic RASM history file post processing script

positional arguments:
  config_file           Input configuration file
  short_term_archive    Case short term archive directory

optional arguments:
  -h, --help            show this help message and exit
  --output_preset OUTPUT_PRESET
                        Input configuration preset
  --processed_dir PROCESSED_DIR
                        Input configuration file
                        (default=$WORKDIR/processed/$RUN/$COMPONENT)
  -np NUMOFPROC, --numofproc NUMOFPROC
                        Number of processors used to run job

rasm_analysis
-----------------------
usage: rasm_analysis [-h] config

positional arguments:
  config      Input Configuration File

optional arguments:
  -h, --help  show this help message and exit
