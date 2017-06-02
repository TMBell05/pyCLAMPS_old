# Standard Libraries
import csv

# 3rd Party Libs
import numpy as np
import xarray


def concat_files(files, concat_dim='time'):
    """
    concatenates files based on a certain dimension and returns the data in a dictionary
    :param files: Files to concatenate
    :param concat_dim: Dimension to concatenate alone
    :return: Dictionary of  var_name:data_array
    """
    arr = xarray.open_mfdataset(files, autoclose=True, concat_dim=concat_dim, decode_times=False)
    data = {}

    for key in arr.keys():
        data[key] = arr[key].values
        data[key] = np.ma.masked_where(np.isnan(data[key]), data[key])

    arr.close()

    return data


def read_wind_profile(profile):
    """
    Reads in a wind profile processed by the Halo lidar
    :param profile: File to process
    :return:
    """

    with open(profile, 'rb') as f:
        reader = csv.reader(f, delimiter=' ', skipinitialspace=True)

        for i, line in enumerate(reader):
            if i == 0:
                num_hgt = int(line[0])
                z = np.zeros(num_hgt)
                dir = np.zeros(num_hgt)
                spd = np.zeros(num_hgt)
            else:
                z[i-1] = line[0]
                dir[i-1] = line[1]
                spd[i-1] = line[2]

    return z, dir, spd
