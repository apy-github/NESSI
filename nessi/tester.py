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
