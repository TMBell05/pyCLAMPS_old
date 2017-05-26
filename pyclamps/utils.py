"""
Collection of utility functions for processing CLAMPS data
"""
import numpy as np
import xarray
import smtplib
import os

from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email.mime.text import MIMEText
from email.mime.image import MIMEImage
from email.utils import COMMASPACE, formatdate
from email import Encoders
import ConfigParser

# Constants
Re = 6371000
R43 = Re * 4.0 / 3.0


def concat_files(files, concat_dim='time'):
    arr = xarray.open_mfdataset(files, autoclose=True, concat_dim=concat_dim, decode_times=False)
    data = {}

    for key in arr.keys():
        data[key] = arr[key].values

    arr.close()

    return data


def jet_max(t, z, ws, buf=12):
    # need time, height  (assumes constant dz), and wind speed, buffer of hieght for max to move in time
    print('Finding max...')
    z_llj = np.full_like(t, 0)  # store jet max heights in time
    ws_llj = np.full_like(t, 0)  # store jet max ws in time
    for i in range(0, len(t)):
        wsmax = np.nanmax(ws[i, :22])
        if np.isnan(wsmax):
            z_llj[i] = np.nan
            ws_llj[i] = np.nan
        else:
            wsmax_index = np.where(ws[i, :22] == wsmax)[0][0]
            z_llj[i] = z[wsmax_index]
            ws_llj[i] = wsmax
            ##TIM METHOD######
            #   dz  = z[1]-z[0]
            #   for i in range(0,len(t)):
            #       shear = np.full_like(ws[i,:],np.nan)
            #       for j in range(0,len(z)-1):
            #           shear[j] = (ws[i,j+1]-ws[i,j])/dz
            #       if np.all(np.isnan(shear)):
            #           max_index = 0
            #       else:
            #           max_index = np.where(shear<(-.001))[0][0]-1 #after Tim's threshold, bonin 2015 in BLM
            #       #print max_index
            # z_llj[i] = z[max_index]
            #       ws_llj[i] = ws[i,max_index]

    return z_llj, ws_llj


def list_to_masked_array(in_list, mask_value):
    a = np.array(in_list)
    return np.ma.masked_where(a == mask_value, a)


def ray_height(rng, elev, H0=0, R1=R43):
    """
    Center of radar beam height calculation.
    Rinehart (1997), Eqn 3.12, Bech et al. (2003) Eqn 3
    INPUT::
    -----
    r : float
        Range from radar to point of interest [m]
    elev : float
        Elevation angle of radar beam [deg]
    H0 : float
        Height of radar antenna [m]
    R1 : float
        Effective radius
    OUTPUT::
    -----
    H : float
        Radar beam height [m]
    USAGE::
    -----
    H = ray_height(r,elev,H0,[R1=6374000.*4/3])
    NOTES::
    -----
    If no Effective radius is given a "standard atmosphere" is assumed,
       the 4/3 approximation.
    Bech et al. (2003) use a factor ke that is the ratio of earth's radius
       to the effective radius (see r_effective function) and Eqn 4 in B03
    """

    # Convert earth's radius to km for common dN/dH values and then
    # multiply by 1000 to return radius in meters
    hgt = np.sqrt(rng ** 2 + R1 ** 2 + 2 * rng * R1 * np.sin(np.deg2rad(elev)))
    hgt = hgt - R1 + H0

    return hgt


def send_mail(send_from, send_to, subject, text, files=None,
              data_attachments=None, server="smtp.mail.me.com", port=587,
              tls=True, html=False, images=None,
              username=None, password=None,
              config_file=None, config=None):

    if files is None:
        files = []

    if images is None:
        images = []

    if data_attachments is None:
        data_attachments = []

    if config_file is not None:
        config = ConfigParser.ConfigParser()
        config.read(config_file)

    if config is not None:
        server = config.get('smtp', 'server')
        port = config.get('smtp', 'port')
        tls = config.get('smtp', 'tls').lower() in ('true', 'yes', 'y')
        username = config.get('smtp', 'username')
        password = config.get('smtp', 'password')

    msg = MIMEMultipart('related')
    msg['From'] = send_from
    msg['To'] = send_to if isinstance(send_to, basestring) else COMMASPACE.join(send_to)
    msg['Date'] = formatdate(localtime=True)
    msg['Subject'] = subject

    msg.attach( MIMEText(text, 'html' if html else 'plain') )

    for f in files:
        part = MIMEBase('application', "octet-stream")
        part.set_payload( open(f,"rb").read() )
        Encoders.encode_base64(part)
        part.add_header('Content-Disposition', 'attachment; filename="%s"' % os.path.basename(f))
        msg.attach(part)

    for f in data_attachments:
        part = MIMEBase('application', "octet-stream")
        part.set_payload( f['data'] )
        Encoders.encode_base64(part)
        part.add_header('Content-Disposition', 'attachment; filename="%s"' % f['filename'])
        msg.attach(part)

    for (n, i) in enumerate(images):
        fp = open(i, 'rb')
        msgImage = MIMEImage(fp.read())
        fp.close()
        msgImage.add_header('Content-ID', '<image{0}>'.format(str(n+1)))
        msg.attach(msgImage)

    smtp = smtplib.SMTP(server, int(port))
    if tls:
        smtp.starttls()

    if username is not None:
        smtp.login(username, password)
    smtp.sendmail(send_from, send_to, msg.as_string())
    smtp.close()

