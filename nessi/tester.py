#!/usr/bin/env python3
# coding: utf-8

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
def test():

  import numpy as np
  from nessi import integrator as nss

  path = "%s/%s" % ("/".join(__file__.split("/")[0:-1]), "data", )
  data = np.load("%s/syn_clv.npz" % (path,))

  wav = data['wave']
  dc = data["spectra"]
  mu = data['mu']
  clv = data['clv']

  sol = data['integration']

  saas = nss.sun_as_a_star()
  saas.update_clv(mu,wav,clv,wav,dc)
  saas.update_vrot(0.,0.)
  test_si = saas.get_integration()

  assert np.nanstd(test_si-sol) < 1.e-16, "\n\tTest Failed!"
  print("\n\tCongratulations! Test passed!")

  return 0
#
#
#
def load_data():

  import numpy as np

  path = "%s/%s" % ("/".join(__file__.split("/")[0:-1]), "data", )
  data = np.load("%s/syn_clv.npz" % (path,))

  wav = data['wave']
  dc = data["spectra"]
  mu = data['mu']
  clv = data['clv']

  return wav, dc, mu, clv

def test_fov_accuracy():

  from nessi import integrator as nss
  from nessi.tester import load_data

  import numpy as np

  from matplotlib import pyplot as pl
  pl.ioff()

  ntest = 1

  wav, dc, mu, clv = load_data()
  saas = nss.sun_as_a_star(nr=101)
  saas.update_clv(mu,wav,clv,wav,dc)
  saas.update_vrot(0.,0.)
  test_si = saas.get_integration()
  
  fgs = []
  axs = []

  for itn in range(ntest):
    fg, ax = pl.subplots(clear=True, num=itn+1)
    fgs.append(fg)
    axs.append(ax)

  for itn in range(ntest):
    axs[itn].plot(test_si, label="Sun as a star reference spectra")
  
  # Example 0: Block the whole Sun (minimum accuracy test):
  

  ax = axs[0]
  
  # Example 0b: Block the whole Sun (cartesian limited accuracy test):
  
  for i in [10, 100, 1000,]:
  
    x1d = np.linspace(-1,1,i+11)                     
    y1d = np.linspace(-1,1,i+21)                     
    xx, yy = np.meshgrid(x1d,y1d)                      
  
    rr = xx**2+yy**2
    xx[rr>1] = np.nan
    yy[rr>1] = np.nan
  
    fov_spectra = (xx*0.)[None,1:,1:]+np.ones((test_si.size,))[:,None,None]
    fov_spectra *= 0
  
    tmp = saas.get_diff_spectra_fov(xx,yy,fov_spectra)
  
    ax.plot(test_si + tmp, label="Total Eclipse (%is)" % (i,))
    resid = test_si + tmp
    print("\tn=%i ; rel. error: %.4e" % (i, np.nansum(resid*resid)/np.nansum(test_si*test_si),))
  
  ax.legend()
  pl.show()

  return


def test_fov_flat_emission():

  from nessi import integrator as nss
  from nessi.tester import load_data

  import numpy as np

  from matplotlib import pyplot as pl
  pl.ioff()

  ntest = 1

  wav, dc, mu, clv = load_data()
  saas = nss.sun_as_a_star(nr=101)
  saas.update_clv(mu,wav,clv,wav,dc)
  saas.update_vrot(0.,0.)
  test_si = saas.get_integration()
  
  fgs = []
  axs = []

  for itn in range(ntest):
    fg, ax = pl.subplots(clear=True, num=itn+1)
    fgs.append(fg)
    axs.append(ax)

  for itn in range(ntest):
    axs[itn].plot(test_si, label="Sun as a star reference spectra")
  
  # Example 1: Constant white emission spectra from a small chunk:
  
  x1d = np.linspace(0.4,0.5,111)                     
  y1d = np.linspace(0.4,0.5,121)                     
  xx, yy = np.meshgrid(x1d,y1d)                      
  
  fov_spectra = (xx*0.)[None,1:,1:]+np.ones((test_si.size,))[:,None,None]
  
  tmp = saas.get_diff_spectra_fov(xx,yy,fov_spectra)
  
  ax = axs[0]
  ax.plot(test_si + tmp, label="FOV (1st quadrant)")
  ax.legend()

  pl.show()
  return
  

def test_east_west_hemispheres():

  from nessi import integrator as nss
  from nessi.tester import load_data

  import numpy as np

  from matplotlib import pyplot as pl
  pl.ioff()

  ntest = 2

  wav, dc, mu, clv = load_data()
  saas = nss.sun_as_a_star(nr=101)
  saas.update_clv(mu,wav,clv,wav,dc)
  saas.update_vrot(0.,0.)
  test_si = saas.get_integration()
  
  fgs = []
  axs = []

  for itn in range(ntest):
    fg, ax = pl.subplots(clear=True, num=itn+1)
    fgs.append(fg)
    axs.append(ax)

  for itn in range(ntest):
    axs[itn].plot(test_si, label="Sun as a star reference spectra")
  
  
  # Example 2: Block the East side:
  
  x1d = np.linspace(-1,0.,111)                     
  y1d = np.linspace(-1,1.,121)                     
  xx, yy = np.meshgrid(x1d,y1d)                      
  
  rr = xx**2+yy**2
  xx[rr>1] = np.nan
  yy[rr>1] = np.nan
  
  fov_spectra = (xx*0.)[None,1:,1:]+np.ones((test_si.size,))[:,None,None]
  fov_spectra *= 0
  
  tmp = saas.get_diff_spectra_fov(xx,yy,fov_spectra)
  
  ax = axs[0]
  ax.plot(test_si + tmp, label="East Eclipse")
  mtest = saas._get_mspectra(mask=np.int16(saas.pts[0,:]>0))
  ax.plot(mtest, label="East Eclipse (brute-force check)")
  ax.legend()
  
  # Example 3: Block the West side:
  
  x1d = np.linspace(0.,1.,111)                     
  y1d = np.linspace(-1,1.,121)                     
  xx, yy = np.meshgrid(x1d,y1d)                      
  
  rr = xx**2+yy**2
  xx[rr>1] = np.nan
  yy[rr>1] = np.nan
  
  fov_spectra = (xx*0.)[None,1:,1:]+np.ones((test_si.size,))[:,None,None]
  fov_spectra *= 0
  
  tmp = saas.get_diff_spectra_fov(xx,yy,fov_spectra)
  
  ax = axs[1]
  #uts.nplot(test_si + tmp,noerase=1, plkwargs={"label":"West Eclipse"})
  ax.plot(test_si + tmp, label="West Eclipse")
  mtest = saas._get_mspectra(mask=np.int16(saas.pts[0,:]<0))
  ax.plot(mtest, label="West Eclipse (brute-force check)")
  ax.legend()
 
  pl.show()
  return

def load_data_fov(nw):

  import numpy as np

  nx = 110
  ny = 120
  x1d = np.linspace(0.4,0.5,nx+1) # Notice! +1
  y1d = np.linspace(0.4,0.5,ny+1) # Notice! +1
  xx, yy = np.meshgrid(x1d,y1d)
  
  fov_spectra = np.ones((nw,ny,nx),dtype="f8")

  return xx, yy, fov_spectra

