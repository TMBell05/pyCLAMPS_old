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
from pyproj import Proj
from scipy.interpolate import RectBivariateSpline

# Custom libs
from pyclamps.__init__ import tower_ids, perdigao_clamps_lat, perdigao_clamps_lon
from pyclamps.plotting import rhi_plot
from pyclamps.utils import rotate, get_tower_from_ucar_nc


parser = ArgumentParser()
parser.add_argument('-i', dest='in_file', help="RHI files")
parser.add_argument('-o', dest='out_dir', help='Directory to store images')
parser.add_argument('-t', dest='terrain', help='Terrain file for cross section and tower height')
parser.add_argument('--towers-file', dest='towers_file', help='Directory containing the tower data')
args = parser.parse_args()


# Make an interpolation and projection so we can get the height of the terrain using lat lon for later
with netCDF4.Dataset(args.terrain) as terrain_nc:
    proj = Proj(terrain_nc['x'].proj4)
    terrain_height = RectBivariateSpline(terrain_nc['x'][0, :], terrain_nc['y'][:, 0], terrain_nc['z'][:].transpose())

# Now we need to pull in the tower data
cs = 'valley'
towers = []
for id in tower_ids[cs]:
    print id
    try:
        tower = get_tower_from_ucar_nc(args.towers_file, id)
        towers.append(tower)
    except Exception as e:
        print e

# Go through the files and process them
nc = netCDF4.Dataset(args.in_file)

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

# Loop through the scans in the file and plot them
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
    image_name = join(args.out_dir, image_name)

    # Testing horiz component view
    tmp_elev = np.meshgrid(elev, rng_m)[0].transpose()
    h_vel = vel * np.cos(tmp_elev)
    # if az > 180.: az *= -1.

    del tmp_elev

    # If the image has already been created, don't need to do anything
    if isfile(image_name):
        continue

    # Make the RHI plot
    ax, z_0 = rhi_plot(elev, rng_m, h_vel, az, time, vmax=10, vmin=-10, xlim=(-1000, 1000), ylim=(300, 700),
                       terrain_file=args.terrain)

    '''Now we need to plot the tower data'''

    for tower in towers:

        # Get the base height of the tower and change the coord system to be relative to the radar
        x, y = proj(tower.lon, tower.lat)
        x_0, y_0 = proj(perdigao_clamps_lon, perdigao_clamps_lat)

        x_0 -= x
        y_0 -= y
        z_0 = terrain_height(x, y)

        # Get the closest time
        ind = (np.abs(tower.time - time)).argmin()

        for height in tower['u'].heights:
            # Rotate the coordinate system so +x is in the direction of the beam
            if az > 180.:
                points = rotate(x_0, y_0, z_0, np.deg2rad(-270. - (az - 360.)), 0, 0)
                result = rotate(tower['u'][height][ind], tower['v'][height][ind], tower['w'][height][ind],
                                np.deg2rad(-270. - (az - 360.)), 0, 0)
            else:
                points = rotate(x_0, y_0, z_0, np.deg2rad(-270. - az), 0, 0)
                result = rotate(tower['u'][height][ind], tower['v'][height][ind], tower['w'][height][ind],
                                np.deg2rad(-270. - az), 0, 0)
            u, w = result[:, 0], result[:, 2]
            x, z_0 = points[:, 0], points[:, 2]

            # Plot the vector
            q = ax.quiver(x, height + z_0, u, w)#, scale_units='inches', scale=10)

    '''Save the figure'''
    plt.quiverkey(q, .9, .8, 5, '5 $ms^{-1}$')
    plt.savefig(image_name)
    plt.close()

nc.close()
