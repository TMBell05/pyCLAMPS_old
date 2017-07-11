import matplotlib; matplotlib.use('agg')

# Python standard libs
import logging
import os
from argparse import ArgumentParser
from datetime import datetime
from os.path import join, isfile

# Downloaded libs
import matplotlib.pyplot as plt
import netCDF4
import numpy as np

# Custom libs
from pyclamps.plotting import rhi_plot
from pyclamps.io import concat_files
from pyclamps.utils import rotate


parser = ArgumentParser()
parser.add_argument('-i', dest='in_files', nargs='*')
parser.add_argument('-o', dest='out_dir')
parser.add_argument('-t', dest='terrain', default=None)
parser.add_argument('-v', dest='vads', nargs='*', default=None)
args = parser.parse_args()

if args.vads is not None:
    # Get the vads
    vads = concat_files(args.vads, concat_dim='time')

    # Get the times now to save time later
    vad_times = []
    for t in vads['time']: vad_times.append(datetime.utcfromtimestamp(t))
    vad_times = np.array(vad_times)

else:
    vads = None

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
        # if az > 180.: h_vel = -h_vel

        del tmp_elev

        # If the image has already been created, don't need to do anything
        if isfile(image_name):
            continue

        # Make the RHI plot
        ax, z_0 = rhi_plot(elev, rng_m, h_vel, az, time, vmax=10, vmin=-10, xlim=(-2300, 2300), ylim=(0, 1800),
                           terrain_file=args.terrain)

        if vads is not None:
            # Get the closest time
            ind = (np.abs(vad_times-time)).argmin()
            # Rotate the coordinate system so +x is in the direction of the beam

            if az > 180.:
                result = rotate(vads['u'][ind], vads['v'][ind], vads['w'][ind], np.deg2rad(-270. - (az - 360.)), 0, 0)
            else:
                result = rotate(vads['u'][ind], vads['v'][ind], vads['w'][ind], np.deg2rad(-270. - az), 0, 0)
            u, w = result[:, :, 0], result[:, :, 2]

            # Plot the vector
            q = ax.quiver(0, vads['hgt'][ind] + z_0, u, 0, scale_units='inches', scale=10)
            plt.quiverkey(q, .9, .8, 5, '5 $ms^{-1}$')

            # plot barbs on the RHS
            ax.barbs(np.repeat(2000, vads['hgt'][ind].size), vads['hgt'][ind] + z_0, vads['u'][ind], vads['v'][ind])

        plt.savefig(image_name)
        plt.close()

    nc.close()
