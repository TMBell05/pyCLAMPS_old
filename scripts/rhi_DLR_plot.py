import matplotlib; matplotlib.use('agg')

from argparse import ArgumentParser
from datetime import datetime
from os.path import join, isfile, basename

import matplotlib.pyplot as plt
import netCDF4
import numpy as np
from pyclamps.plotting import rhi_plot

terrain_nc = '/Users/tbell/Data/perdigao/terrain/newa_perdigao_map_topo.nc'

dlr_lat_lons = {'DLR85': (39.714533, -7.731142),
                'DLR86': (39.712467, -7.735108),
                'DLR89': (39.709867, -7.746794)}

parser = ArgumentParser()
parser.add_argument('-i', dest='in_files', nargs='*')
parser.add_argument('-o', dest='out_dir')
parser.add_argument('-t', dest='terrain', default=terrain_nc)
args = parser.parse_args()
for f in args.in_files:
    nc = netCDF4.Dataset(f)

    # Get the coords of the appropriate lidar
    coord_key = basename(f).split('_')[0]
    lon, lat = dlr_lat_lons[coord_key]

    for key in nc.groups.keys():
        sweep = nc.groups[key]

        # Get the datetime of the first ray
        time = datetime.utcfromtimestamp(sweep['time'][:][0])

        # Get the grid figured out
        elev = np.deg2rad(sweep['elevation'])
        rng_m = sweep['range'][:][0]
        az = np.round(np.mean(sweep['azimuth']))

        # Get the data
        vel = sweep['rws']
        # intensity = nc['intensity'][
        # vel = np.ma.masked_where(intensity < 1.01, vel)

        # Figure out the image name
        image_name = '{}_RHI_horiz_{}_{}.png'.format(coord_key, az, time.strftime("%Y%m%d_%H%M"))
        # image_name = 'RHI_{}_{}.png'.format(az, time.strftime("%Y%m%d_%H%M"))
        image_name = join(args.out_dir, image_name)

        # Testing horiz component view
        tmp_elev = np.meshgrid(elev, rng_m)[0].transpose()
        h_vel = vel * np.cos(tmp_elev)
        if az > 180.: h_vel = -h_vel
        del tmp_elev

        # if not isfile(image_name):
        print(image_name)

        rhi_plot(elev, rng_m, -h_vel, az, time, vmax=10, vmin=-10, xlim=(-500, 2000), ylim=(0, 1500),
                 terrain_file=args.terrain, lat_0=lon, lon_0=lat)
        plt.savefig(image_name)
        plt.close()

    nc.close()
