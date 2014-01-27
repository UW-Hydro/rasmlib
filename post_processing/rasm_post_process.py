#!/usr/bin/env python
"""RASM Postprocessing script for VIC/RVIC/Coupler fields

The script does the following to a series of netCDF files:
 * Reads command line options
 * Reads the configuration file
 * Can adjust the timestamp
 * Calculates any or all of the following averages
    - timeseries of mean monthly diurnal cycles from hourly data files
    - timeseries of daily means
    - timeseries of daily means

Note that any existing file will be overwritten.

The script relies on the nco and nco.py utilities for some of its steps:
 * http://nco.sourceforge.net
"""

from __future__ import print_function
import argparse
import glob
import os
import means
import share
from adjust_timestamp import adjust_timestamp


def main():
    """ """
    # ---------------------------------------------------------------- #
    # Read command Line
    config_file, short_term_archive, output_preset, \
        processed_dir, numofproc = process_command_line()
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Initialize
    files, config_dict = process_init(config_file, short_term_archive,
                                      output_preset, processed_dir)

    # config_dict['options']['numofproc'] = numofproc
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Run
    results = process_run(files, config_dict)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Final
    process_final(results, config_dict)
    # ---------------------------------------------------------------- #
    return


def process_init(config_file, short_term_archive,
                 output_preset, processed_dir):
    # ---------------------------------------------------------------- #
    # Read Configuration files
    config_dict = share.read_config(config_file)

    if output_preset not in config_dict:
        raise ValueError('Output_preset is not in the configuration file, \
                         check your spelling')
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Setup Directory Structures
    casename = os.path.basename(os.path.normpath(short_term_archive))

    if processed_dir:
        outdir = processed_dir
    else:
        if "WORKDIR" in os.environ:
            workdir = os.environ.get('WORKDIR')
            procdir = os.path.join(workdir, 'processed')
            outdir = os.path.join(procdir, casename)
        else:
            raise ValueError('ERROR: Either provide an path to \
                             --processed_dir or set your $WORKDIR \
                             environment variable')

    dirs = ['temp', 'monthly_mean_timeseries']
    if output_preset not in 'monthly':
        dirs.append('daily_mean_timeseries')
    if output_preset == 'hourly':
        dirs.append('monthly_mean_diurnal_cycle')

    processed_case_dir = os.path.join(outdir, casename)
    processed_comp_dir = os.path.join(processed_case_dir,
                                      config_dict['options']['component'])
    directories = share.make_directories(processed_comp_dir, dirs)

    config_dict['options']['directories'] = directories
    config_dict['options']['directories']['processed_case_dir'] = processed_case_dir
    config_dict['options']['directories']['processed_comp_dir'] = processed_comp_dir
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Get list of files
    format = config_dict['format'][output_preset]
    format = format.replace('{CASE}', casename)
    globformat = format.replace('%Y', '????')
    globformat = globformat.replace('%m', '??')
    globformat = globformat.replace('%d', '??')
    globformat = globformat.replace('%s', '?????')

    hist_dir = os.path.join(short_term_archive,
                            config_dict['options']['subpath'])
    print('Glob format: {0}. Histdir: {1}'.format(globformat, hist_dir))
    files = glob.glob(os.path.join(hist_dir, globformat))
    print('Found {0} files from glob'.format(len(files)))

    filelist = [share.Histfile(filename=f, fname_format=format) for f in files]
    filelist.sort(key=lambda r: r.filedate)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Pack some options into the config dict
    config_dict['options']['casename'] = casename
    config_dict['options']['fname_format'] = format
    config_dict['options']['output_preset'] = output_preset
    config_dict['options']['timestep'] = output_preset
    # ---------------------------------------------------------------- #

    return filelist, config_dict


def process_run(filelist, config_dict):
    """ """
    results = {}
    # ---------------------------------------------------------------- #
    # Unpack some useful stuff
    options = config_dict['options']
    directories = config_dict['options']['directories']
    fname_format = options['fname_format']
    output_preset = config_dict['options']['output_preset']
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Adjust timestep if needed
    if abs(options['timestamp_offset']) > 0:
        filelist = adjust_timestamp(filelist,
                                    timestep=output_preset,
                                    nsteps=options['timestamp_offset'],
                                    fname_format=fname_format,
                                    destdir=directories['temp'],
                                    calendar=options['calendar'])
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Calculate monthly_mean_diurnal_cycle
    if output_preset == 'hourly':
        function = 'monthly_mean_diurnal_cycle'
        varlist = config_dict[output_preset]['monthly_mean_diurnal_cycle']
        results[function] = means.monthly_mean_diurnal_cycle(filelist, options,
                                                             variables=varlist)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Calculate daily mean timeseries
    if output_preset in ['daily', 'hourly']:
        function = 'daily_mean_timeseries'
        varlist = config_dict[output_preset]['daily_mean_timeseries']
        results[function] = means.daily_mean_timeseries(filelist, options,
                                                        variables=varlist)
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Calculate monthly mean timeseries
    function = 'monthly_mean_timeseries'
    varlist = config_dict[output_preset]['monthly_mean_timeseries']
    results[function] = means.monthly_mean_timeseries(filelist, options,
                                                      variables=varlist)
    # ---------------------------------------------------------------- #

    return results


def process_final(results, config_dict):
    """ """

    print('Finalizing rasm_post_process.py\n')
    output_preset = config_dict['options']['output_preset']
    component = config_dict['options']['component']
    processed_case_dir = config_dict['options']['directories']['processed_case_dir']
    processed_comp_dir = config_dict['options']['directories']['processed_comp_dir']

    # ---------------------------------------------------------------- #
    # Make a compressed tar.gz file from the archive
    tarfile_name = os.path.join(processed_case_dir,
                                "{0}.tar.gz".format(component))
    share.make_tarfile(tarfile_name, processed_comp_dir)
    # ---------------------------------------------------------------- #
    # ---------------------------------------------------------------- #

    # Clean up
    if config_dict['options']['clean']:
        tempdir = config_dict['options']['directories']['temp']
        share.clean_dir(tempdir)
        print("Cleaned up temporary directory: {}".format(tempdir))
    # ---------------------------------------------------------------- #

    # ---------------------------------------------------------------- #
    # Report on what we did
    print('Results from rasm_post_process.py:')
    for key, val in results.iteritems():
        print("{}: {}".format(key, val))
    print("\n")

    print('Variables included in each level of analysis were:')
    for key, val in config_dict[output_preset].iteritems():
        print("{}: {}".format(key, val))
    print("\n")

    print('Completed rasm_post_process.py.  Options were:')
    for key, val in config_dict['options'].iteritems():
        print("{}: {}".format(key, val))
    print("\n")

    print('tarfile is here: {0}'.format(tarfile_name))
    print("\n")

    print('done')
    # ---------------------------------------------------------------- #

    return


def process_command_line():
    """
    Read Command Line Arguments
    """
    # Parse arguments
    parser = argparse.ArgumentParser(description='Generic RASM history file \
                                     post processing script')
    parser.add_argument("config_file", type=str,
                        help="Input configuration file")
    parser.add_argument("short_term_archive", type=str,
                        help="Case short term archive directory")
    parser.add_argument("--output_preset", type=str,
                        help="Input configuration preset",
                        default='daily')
    parser.add_argument("--processed_dir", type=str,
                        help="Input configuration file \
                        (default=$WORKDIR/processed/$RUN/$COMPONENT)",
                        default=None)
    parser.add_argument("-np", "--numofproc", type=int,
                        help="Number of processors used to run job", default=1)

    args = parser.parse_args()

    return args.config_file, args.short_term_archive, args.output_preset,\
        args.processed_dir, args.numofproc


# -------------------------------------------------------------------- #
if __name__ == "__main__":
    main()
# -------------------------------------------------------------------- #
