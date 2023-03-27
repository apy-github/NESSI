"""
    NESSI –Numerical Empirical Sun as a Star Integrator– is a numerical code
    that calculates the integrated solar spectrum provided with a 
    center-to-limb-variation and a disk-center spectrum.
    Copyright (C) 2023  Adur Pastor Yabar

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful, 
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import numpy as _np
from scipy.interpolate import RectBivariateSpline as _RectBivariateSpline
from scipy.interpolate import interp1d as _interp1d

def get_grid(nn):

  rarr = _np.arange(nn) / (nn-1.)
  rad_cen = (rarr[0:-1] + rarr[1:]) / 2.
  rad_wid = (rarr[1:] - rarr[0:-1])

  dr = _np.nanmean(rad_wid)

  for itrn, itrv in enumerate(rad_cen):

    inn = _np.int32(_np.round(2.*_np.pi*itrv/dr))
    ita = _np.arange(inn) / (inn-1.) * 2. * _np.pi
    if (itrn==0):
      azi_cen = (ita[0:-1] + ita[1:]) / 2.
      azi_wid = (ita[1:] - ita[0:-1])
      rad_cen = _np.ones((ita.size-1,)) * itrv
      drad = _np.ones((ita.size-1,)) * ((itrv + dr/2.)**2-(itrv - dr/2.)**2)
    else:
      azi_cen = _np.hstack([azi_cen, (ita[0:-1] + ita[1:]) / 2.])
      azi_wid = _np.hstack([azi_wid, (ita[1:] - ita[0:-1])])
      rad_cen = _np.hstack([rad_cen, _np.ones((ita.size-1,)) * itrv])
      drad = _np.hstack([drad, _np.ones((ita.size-1,)) * ((itrv + dr/2.)**2-(itrv - dr/2.)**2)])

  area = azi_wid / 2. * drad#(erad**2-irad**2)

  nn = rad_cen.size
  pts = _np.zeros((3,nn), dtype='f8')
  pts[0,:] = rad_cen[:] * _np.cos(azi_cen[:])
  pts[1,:] = rad_cen[:] * _np.sin(azi_cen[:])
  pts[2,:] = _np.sqrt(1.-_np.nansum(pts[0:2,:]**2, axis=(0,)))

  return rad_cen, area, pts
#
#
#Taken from ISP library:
#  https://github.com/ISP-SST/ISPy/tree/master/ISPy
#
from astropy.table import Table
from astropy.io import fits
#
def _limbdarkening(wave, mu=1.0, nm=False):
    """
    Return limb-darkening factor given wavelength and viewing angle
    mu=cos(theta)
    Parameters
    ----------
    wave : float or array_like
        scalar or 1D array with wavelength(s).
    mu : float, optional
        cosine of heliocentric viewing angle (default 1.0 -> disc centre)
    nm : bool, optional
        input wavelength units are nanometers (default False)
    Returns
    -------
    factor : float or array_like
        scaling factor(s) to be applied for given input wavelengths. Has as many
        elements as `wave`.
    Example
    -------
    >>> factor = limbdarkening(630.25, mu=0.7, nm=True)
    :Author:
        Gregal Vissers (ISP/SU 2020)
    """

    #this_dir, this_filename = os.path.split(__file__)
    #DATA_PATH = os.path.join(this_dir, "../data/limbdarkening_Neckel_Labs_1994.fits")

    wave = np.atleast_1d(wave)  # Ensure input is iterable

    #table = Table(fits.getdata(DATA_PATH))
    table = Table(fits.getdata("limbdarkening_Neckel_Labs_1994.fits"))
    wavetable = np.array(table['wavelength'])
    if nm is False:
        wavetable *= 10.

    # Get table into 2D numpy array
    Atable = np.array([ table['A0'], table['A1'], table['A2'],
        table['A3'], table['A4'], table['A5'] ])

    factor = np.zeros((wave.size), dtype='float64')
    for ii in range(6):
      Aint = np.interp(wave, wavetable, Atable[ii,:])
      factor += Aint * mu**ii

    return factor
#

def _get_rot_coeffs(ref="s90"):

  if (ref.lower()=="n51"):
    coeff = [14.368, -2.69, 0.00]
  elif (ref.lower()=="h84"):
    coeff = [14.522, -2.84, 0.00]
  elif (ref.lower()=="h86"):
    coeff = [14.551, -2.87, 0.00]
  elif (ref.lower()=="s84"):
    coeff = [14.050, -1.492, -2.606]
  elif (ref.lower()=="s90"):
    coeff = [14.71, -2.39, -1.78]
  elif (ref.lower()=="k93"):
    coeff = [14.42, -2.00, -2.09]
  elif (ref.lower()=="b01"):
    coeff = [14.60, -3.00, 0.00]
  elif (ref.lower()=="t75"):
    coeff = [14.23, -0.40, 0.00]
  else:
    raise Exception("Unknown solar rotation profile: %s" % (ref))

  return coeff

def _rotate(x, y, z, ang, axis='x'):

  ca = _np.cos(ang)
  sa = _np.sin(ang)

  if (axis=='x'):
    tmp = y * ca + sa * z
    nz = ca * z - sa * y
    ny = tmp * 1.
    nx = x * 1.
  elif (axis=='y'):
    tmp = x * ca + sa * z
    nz = ca * z - sa * x
    nx = tmp * 1.
    ny = y * 1.
  elif (axis=='z'):
    tmp = x * ca + sa * y
    ny = ca * y - sa * x
    nx = tmp * 1.
    nz = z * 1.
  else:
    raise Exception('axis MUST be one of: "x", "y", or "z"')
  return nx, ny, nz

def _get_sun_vrot(p0, b0, n=101, pts=None, rotref="s90"):

  if (type(pts)==type(None)):
    _, _, pts = get_grid(n)

  x = pts[0].copy()
  y = pts[1].copy()
  z = pts[2].copy()

  odims = x.shape
  x = x.flatten()
  y = y.flatten()
  z = z.flatten()

  ang = p0 / 180. * _np.pi

  x, y, z = _rotate(x, y, z, ang, axis='z')

  lat = _np.arctan2(y, _np.sqrt(x**2+z**2))
  lon = _np.arctan2(x, z)

  a, b, c = _get_rot_coeffs(ref=rotref)

  sa2 = _np.sin(lat) ** 2

  omega = a + b * sa2 + c * sa2**2

  rsun = 695700. #km
  rsunred = _np.cos(lat) * rsun

  circ = 2. * _np.pi * rsunred
  res = circ * (omega / 360. / (24.*3600.))   # km/s

  vz = res * _np.cos(lat) * _np.sin(lon)
  vx = res * _np.cos(lat) * _np.cos(lon)
  vy = res * 0.

  ang = b0 / 180. * _np.pi

  vx, vy, vz = _rotate(vx, vy, vz, ang, axis='x')

  return vz.reshape(*odims)

#
#
#
class sun_as_a_star(object):

  def __init__(self, muclv, wavclv, clv, wavdc, dc \
        , nr=51, rotref="s90" \
        , include_red_cclv=False, include_full_cclv=False):

    assert(wavdc.size==dc.size)
    assert(wavclv.size==clv.shape[1])
    assert(muclv.size==clv.shape[0])

    self.nw = dc.size
    self.wave = wavdc.copy()
    self.dc = dc.copy()
    wclv = wavclv.copy()
    rclv = _np.sin(_np.arccos(muclv))
    clv = clv.copy()

    ww = _np.argsort(rclv)
    rclv = rclv[ww] * 1.
    clv = clv[ww,:] * 1.
    ww = _np.argsort(wclv)
    wclv = wclv[ww] * 1.
    clv = clv[:,ww] * 1.

    self.fclv = _RectBivariateSpline(rclv, wclv, clv)#test)
    self.fdc = _interp1d(wavdc, dc, bounds_error=False, fill_value="extrapolate")

    self.rad_cen, self.area, self.pts = get_grid(nr)
    self.cclv_desc = "native"
    if (include_red_cclv==True):
      self.cclv_desc = "reduced"
      print("\t[] Including red cclv representation.")
      mu = _np.cos(_np.arcsin(self.rad_cen))
      self.cclv = _np.zeros((self.wave.size, mu.size), dtype="f8")
      for itw in range(self.wave.size):
        self.cclv[itw,:] = _limbdarkening(self.wave[itw], mu)

    elif (include_full_cclv==True):
      self.cclv_desc = "full"
      print("\t[] Including full cclv representation.")
      self.mu = _np.cos(_np.arcsin(self.rad_cen))

    self.rotref = rotref

    return

  def update_vrot(self, p0, b0, return_vrot=False):

    vrot = _get_sun_vrot(p0, b0, pts=self.pts, rotref=self.rotref)

    c = 2.99792458e5
    self.dwave = self.wave[None,:] * vrot[:,None] / c

    if (return_vrot!=False):
      return vrot

    return

  def get_integration(self):

    res = self.dc * 0.
    tarea = _np.nansum(self.area)
    cclv = 1.
    for itw in range(self.nw):

      itwav = self.wave[itw] - self.dwave[:,itw]
      tmp = self.fclv.ev(self.rad_cen, itwav)
      tmpdc = self.fdc(itwav)
      if (self.cclv_desc=="reduced"):
        cclv = self.cclv
      elif (self.cclv_desc=="full"):
        cclv = _limbdarkening(itwav, self.mu)
      else:
        cclv = 1.

      res[itw] = _np.sum(cclv * tmp * tmpdc * self.area) / tarea

    return res
