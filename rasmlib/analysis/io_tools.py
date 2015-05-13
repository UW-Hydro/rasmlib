"""
rasmlib analysis input/output utilities
"""
import numpy as np
import pandas as pd
from datetime import datetime


def read_r_arctic_net_data(data_file, attrs_file):
    '''read R-ArcticNET dataset into Pandas.DataFrame

    return data, drainage_area, attributes
    '''
    def make_dates(years):
        dates = []
        for year in years:
            for month in range(1, 13):
                dates.append(datetime(year, month, 1))
        return dates

    # open attributes file
    atts = pd.read_csv(attrs_file, engine='python')

    # open observations file(s)
    data = pd.read_csv(data_file, low_memory=False)

    # Flatten monthly structure
    _r_months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    data_series = {}
    basin_area = {}
    for i, point in enumerate(atts['PointID']):
        name = atts['Name'][i].replace("'", '').replace(" ", "_").encode()

        inds = np.nonzero(data['Point_ID'] == point)[0]
        flow_table = np.empty((len(inds), len(_r_months)))

        for j, month in enumerate(_r_months):
            flow_table[:, j] = data[month][inds]

        times = make_dates(data['Year'][inds])

        da = atts['DArea_effective'][i] * 1000. * 1000.  # m2

        flows = flow_table.flatten() / da * 1000 * 86400  # mm/d
        basin_area[name] = da
        data_series[name] = pd.Series(flows, index=times)

    observations = pd.DataFrame(data_series)

    observations[observations <= 0] = np.nan

    return observations, basin_area, atts
