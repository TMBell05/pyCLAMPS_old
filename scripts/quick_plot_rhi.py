import matplotlib; matplotlib.use('agg')

from argparse import ArgumentParser
from datetime import datetime
from os.path import join, isfile

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from pyclamps.plotting import rhi_plot

terrain_nc = '/Users/tbupper90/Data/clamps/newa_perdigao_map_topo.nc'

parser = ArgumentParser()
parser.add_argument('-i', dest='in_files', nargs='*')
parser.add_argument('-o', dest='out_dir')
parser.add_argument('-t', dest='terrain', default=None)
args = parser.parse_args()
for f in args.in_files:
    nc = netCDF4.Dataset(f)

    scans = []
    counter = 0

    # Figure out how many scans there are
    last = None
    for az in nc['azimuth'][:]:
        if last is None:
            last = az
        elif az != last:
            last = az
            scans.append(counter)
            counter += 1
        else:
            scans.append(counter)

    scans = np.asarray(scans)

    for scan in range(0, counter):
        ind = np.where(scans == scan)

        # Get the date
        time = datetime.utcfromtimestamp(nc['base_time'][0] + nc['time_offset'][ind][0])

        # Get the grid figured out
        elev = np.deg2rad(nc['elevation'][ind])
        rng_m = nc['height'][:] * 1e3
        az = np.round(np.mean(nc['azimuth'][ind]))

        # Get the data
        vel = nc['velocity'][ind]
        intensity = nc['intensity'][ind]
        vel = np.ma.masked_where(intensity < 1.01, vel)

        # Figure out the image name
        image_name = 'RHI_horiz_{}_{}.png'.format(az, time.strftime("%Y%m%d_%H%M"))
        # image_name = 'RHI_{}_{}.png'.format(az, time.strftime("%Y%m%d_%H%M"))
        image_name = join(args.out_dir, image_name)

        # Testing horiz component view
        tmp_elev = np.meshgrid(elev, rng_m)[0].transpose()
        h_vel = vel * np.cos(tmp_elev)
        del tmp_elev

        if not isfile(image_name):
            print(image_name)
            rhi_plot(elev, rng_m, h_vel, az, time, vmax=10, vmin=-10, xlim=(-2800, 2800), ylim=(0, 2500),
                     terrain_file=args.terrain)
            plt.savefig(image_name)
            plt.close()

    nc.close()
