# System Libraries
import os

# MPL gets special treatment for serverside processing
import matplotlib; matplotlib.use('agg')
import matplotlib.pyplot as plt

# 3rd party Libraries
import cmocean
import numpy as np
from colormaps import create_colormap

# From this library
from .utils import jet_max

# Make some colormaps for later
cm_ws = create_colormap(1000, base='ncl_helix', name='helix', reverse=False, white=False)
cm_bias = create_colormap(1000, base='ncl_temp_diff_18lev', name='tempdiff', reverse=False, white=False)


def rhi_plot(elev, rng_m, vel, az, time, vmin=-5, vmax=5,
             xlim=(-7500, 7500), ylim=(0, 7500),path=None):
    # ind = np.where(scans == 90)
    #
    # # Get the grid figured out
    # elev = np.deg2rad(nc['elevation'][ind])
    # rng_m = nc['height'][:] * 1e3

    elev, rng_m = np.meshgrid(elev, rng_m)

    # Sort the elevations so it plots right....
    sort = np.argsort(elev, axis=1)[0, :]
    elev = elev[:, sort]
    rng_m = rng_m[:, sort]

    x_m = rng_m * np.cos(elev)
    y_m = rng_m * np.sin(elev)

    vel = vel[sort, :].transpose()

    # Make the plot
    plt.figure(figsize=(10, 5))
    plt.pcolor(x_m, y_m, vel, vmin=vmin, vmax=vmax)
    plt.colorbar()
    plt.xlim(xlim)
    plt.ylim(ylim)
    plt.title('RHI {} {}'.format(az, time.isoformat()))

    if path is not None:
        plt.savefig(path)
    else:
        plt.show()

    plt.close()


def time_height(t, z, field, fieldtype, path, jet_overlay=0, jet_range=.75, wsfield=0, fieldmin=0, fieldmax=25, zmin=0,
                zmax=1501, zstep=500, ):
    # function requires a time axis, height axis, field to plot, field type, & path to save
    # fieldtypes(w,ws,wd,pt,q,bias)

    if fieldtype == 'w':
        cm = cmocean.cm.delta
        cbar_label = 'vertical velocity [m/s]'
    elif fieldtype == 'ws':
        # cm = 'inferno_r'
        cm = cm_ws
        cbar_label = 'windspeed [m/s]'
    elif fieldtype == 'wd':
        cm = cmocean.cm.phase
        cbar_label = 'wind direction [deg]'
    elif fieldtype == 'pt':
        cm = cmocean.cm.thermal
        # cm = 'jet'
        cbar_label = 'potential temperature [K]'
    elif fieldtype == 'q':
        cm = cmocean.cm.haline
        cbar_label = 'q [g/kg]'
    elif fieldtype == 'bias':
        cm = cm_bias
        cbar_label = 'bias [WRF-obs]'

    Field = np.ma.array(field, mask=np.isnan(field))
    c = plt.pcolormesh(t, z, Field.transpose(), vmin=fieldmin, vmax=fieldmax, cmap=cm)
    c.cmap.set_bad('grey', 1.0)
    cb = plt.colorbar(c, use_gridspec=True)
    cb.set_label(cbar_label)
    plt.ylim((zmin, zmax))
    plt.yticks(np.arange(zmin, zmax, zstep))
    plt.ylabel('Height [m]')
    plt.xlabel('Time [UTC]')
    if jet_overlay > 0:
        jet_height = []
        jet_top = []  # jet_range% of max speed above max
        jet_bot = []  # jet_range% of max speed below
        if isinstance(wsfield, np.ndarray) == False:
            print("No wind field provided... set wsfield")
        else:
            z_maxes, ws_maxes = jet_max(t, z, wsfield)
            # print z_maxes
            for i in range(wsfield.shape[0]):
                # print i
                if z_maxes[i] < 40 or np.all(np.isnan(wsfield[i, :100])):
                    jet_top.append(np.nan)
                    jet_bot.append(np.nan)
                    jet_height.append(np.nan)
                else:
                    wsmax = ws_maxes[i]
                    jet_height.append(z_maxes[i])
                    if np.isnan(wsmax):
                        wsmax = 0.
                    wsrange = wsmax * jet_range
                    # print wsrange
                    top_index = np.where(wsfield[i, :100] > wsrange)[0][-1]
                    bot_index = np.where(wsfield[i, :100] > wsrange)[0][0]
                    jet_top.append(z[top_index])
                    jet_bot.append(z[bot_index])

                    # wsmax = np.max(wsfield[i,:100])
                    # if np.isnan(wsmax):
                    #    jet_top.append(np.nan)
                    #    jet_bot.append(np.nan)
                    #    jet_height.append(np.nan)
                    # else:
                    #    wsmax_index = np.where(wsfield[i,:100]==wsmax)[0][0]
                    #    jet_height.append(z[wsmax_index])
                    #    wsmax75 = wsmax*.75
                    #    top_index = np.where(wsfield[i,:100]>wsmax75)[0][-1]
                    #    bot_index = np.where(wsfield[i,:100]>wsmax75)[0][0]
                    #    jet_top.append(z[top_index])
                    #    jet_bot.append(z[bot_index])
            # print(z_maxes)
            # print(jet_top)
            # print(jet_bot)
            plt.plot(t, z_maxes, color='k', ls='-', lw=4)
            plt.plot(t, jet_top, color='k', ls='--', lw=4)
            plt.plot(t, jet_bot, color='k', ls='--', lw=4)
            print("-- plotting jet overlay")
    plt.tight_layout()
    plt.savefig(path)
    print("Figure saved")
    plt.close()


def vad_component_plot(u, v, w, hgt, rmse, r_sq, date, prefix, out_dir,
                       vel_lim=(-30,30), hgt_lim=(0,1000)):
    fig = plt.figure(1, figsize=(25, 7))
    fig.suptitle(date.isoformat())
    plt.subplot(1, 5, 1)
    plt.plot(u, hgt, label=date.strftime("%H:%M"))
    plt.title("U velocity vs Height")
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Height (m)')
    plt.xlim(vel_lim)
    plt.ylim(hgt_lim)
    plt.legend()

    plt.subplot(1, 5, 2)
    plt.plot(v, hgt)
    plt.title("V velocity vs Height")
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Height (m)')
    plt.xlim(vel_lim)
    plt.ylim(hgt_lim)

    plt.subplot(1, 5, 3)
    plt.plot(w, hgt)
    plt.title("W velocity vs Height")
    plt.xlabel('Velocity (m/s)')
    plt.ylabel('Height (m)')
    plt.xlim((-5, 5))
    plt.ylim(hgt_lim)

    plt.subplot(1, 5, 4)
    plt.plot(rmse, hgt)
    plt.title("RMS vs Height")
    plt.xlabel('RMS')
    plt.ylabel('Height (m)')
    plt.xlim((0, 10))
    plt.ylim(hgt_lim)

    plt.subplot(1, 5, 5)
    plt.plot(r_sq, hgt)
    plt.title("$R^2$ vs Height")
    plt.xlabel('$R^2$')
    plt.ylabel('Height (m)')
    plt.xlim((.8, 1))
    plt.ylim(hgt_lim)

    # Create profile base on start time of first profile
    image_name = "{prefix}_{date}.png"
    image_name = image_name.format(prefix=prefix, date=date.strftime("%Y%m%d_%H%M%S"))
    image_name = os.path.join(out_dir, image_name)

    # Save the image
    plt.savefig(image_name)
    plt.clf()
    plt.close()