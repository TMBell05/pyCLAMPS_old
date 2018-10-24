import matplotlib
matplotlib.use('agg')

import netCDF4
import argparse
import os
import matplotlib.pyplot as plt
import numpy as np

from datetime import datetime
from metpy.plots import SkewT
from metpy.units import units
from metpy.calc.thermo import vapor_pressure, saturation_vapor_pressure, dewpoint
from scipy.interpolate import interp1d

from pyclamps.io import concat_files

# Parse the arguments
parser = argparse.ArgumentParser()
parser.add_argument('-i', action='store', dest='aeri_file', help="Aeri file")
parser.add_argument('-o', action='store', dest='out_dir', help="Output directory")
parser.add_argument('-v', action='store', dest='vads', default=None, help="VAD files", nargs='*')
args = parser.parse_args()

# Create the out dir if it doesn't exist
if not os.path.isdir(args.out_dir):
    os.makedirs(args.out_dir)

# Open the netcdf
aeri_data = netCDF4.Dataset(args.aeri_file)

# Get the times from the vads so we can match them to the time in the aeri file
if args.vads is not None:
    vad_data = concat_files(args.vads)
    vad_times = np.asarray([datetime.utcfromtimestamp(d) for d in (vad_data['time'])])

# Create a sonde for each time
for i in range(aeri_data.dimensions['time'].size):

    # Get the date
    date = datetime.utcfromtimestamp(aeri_data['base_time'][:]
                                     + aeri_data['time_offset'][i])

    # Extract the pressure and temperature data
    T = aeri_data['temperature'][i] * units.C
    p = aeri_data['pressure'][i] * units.mbar
    z = aeri_data['height'][:] * 1.e3 + 296 - 484

    # Calculate the dewpoint
    e = vapor_pressure(p, aeri_data['waterVapor'][i]*1e-3)
    Td = dewpoint(e)

    # make the plot
    fig, ax = plt.subplots(1, figsize=(9, 9))
    # skew = SkewT(fig, rotation=0)
    # plt = skew

    # Plot temperature and dewpoint
    ax.plot(T, z, 'r')
    ax.plot(Td, z, 'g')

    if args.vads is not None:
        # Find the closest vad to the aeri time and get the data
        ind = np.abs((vad_times - date)).argmin()
        u = vad_data['u'][ind].filled(np.nan)
        v = vad_data['v'][ind].filled(np.nan)
        h = vad_data['hgt'][ind].filled(np.nan)

        # Eliminate the vad heights above the aeri data
        max_h = max(aeri_data['height'][:]*1e3)
        ind = np.where(h < max_h)
        u = u[ind]
        v = v[ind]
        h = h[ind] + 296 - 484

        # Convert the heights to pressure by interpolating from the aeri derived pressures
        # h2p = interp1d(aeri_data['height'][:]*1e3, aeri_data['pressure'][i])
        # vad_p = h2p(h)

        # Plot the barbs
        ax.barbs(np.full_like(h, 45), h, u, v)

    # Add the relevant special lines
    ax.set_ylim(-200, 1500)
    ax.set_xlim(0, 50)
    ax.grid()
    # skew.plot_mixing_lines()
    ax.set_title("OU AERI " + date.isoformat())

    plt.savefig(os.path.join(args.out_dir, date.strftime('aeri_sonde_bl_%Y%m%d_%H%M%S.png')))
    plt.close()
