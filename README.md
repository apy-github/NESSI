# NESSI -Numerical-Empirical-Sun-as-a-Star-Integrator-

This set of routines perform the numerical integration of the solar disc for a user provided spectral dependant center-to-limb-variation (CLV) and a disk center spectrum.

The CLV may or maynot include the continuum CLV. In the latter, it is possible to use Neckel & Labs 1994 continuum CLV.

Requirements:

  numpy
  scipy

Installation:

  `python3 -m pip install . --user`

Test:

  `from nessi import tester`
  `tester.test()`

Minimal example:
  `from nessi.tester import load_data`
  `from nessi import integrator as saas`
  `wav, dc, mu, clv = load_data()`
  `saas = nss.sun_as_a_star()`
  `saas.update_clv(mu,wav,clv,wav,dc)`
  `saas.update_vrot(0.,0.)`
  `test_si = saas.get_integration()`

