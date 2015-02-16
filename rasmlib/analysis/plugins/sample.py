"""
2x2 plotting analysis for 2 datasets pluggin

|  -------------  |  -------------  |
|  contours from  |  pcolor from    |
|  both datasets  |  dataset 1      |
|  -------------  |  -------------  |
|  pcolor diff    |  pcolor from    |
|  both datasets  |  dataset 2      |
|  -------------  |  -------------  |

colorbar location = bottom
"""


class SampleException(Exception):
    pass

_NDATASETS = 2
_NPANNELS = 4


def run(cases, compares, domain, **kwargs):
    """plugin run function"""

    case_names = cases.keys()
    compare_names = compares.keys()

    dsets = cases.values()+compares.values()

    if len(dsets) != _NDATASETS:
        raise SampleException('Incorrect number of datasets provided')

    # get_monthly_means(*dsets)
    # get_seasonal_means()
    # get_annual_means()
    # get_full_means()

    return


def __plot():

    return
