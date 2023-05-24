
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
