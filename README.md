# NESSI -Numerical-Empirical-Sun-as-a-Star-Integrator-

This set of routines perform the numerical integration of the solar disc for a user provided spectral dependant center-to-limb-variation (CLV) and a disk center spectrum.

The CLV may or maynot include the continuum CLV. In the latter, it is possible to use Neckel & Labs 1994 continuum CLV.

## Requirements:

  ```
  numpy
  scipy
  ```

## Installation:

  ```
  python3 -m pip install . --user
  ```

## Test:

  ```
  from nessi import tester
  tester.test()
  ```

## Minimal example:

  Load some default data to get use to the data types required:

  ```
  from nessi.tester import load_data
  
  wav, dc, mu, clv = load_data()
  ```

  Import the python module that performs the integration:
  
  ```
  from nessi import integrator as saas
  ```
  
  Initializate the main class required to perform the resolved solar disk integration:
  
  ```
  saas = nss.sun_as_a_star(nr=101)
  ```
  
  "nr" specifies the number of radial segments and azimuthal segments to use in the computation.
  
  Load the observed/synthetic center to limb variation (clv) and the disk center reference spectrum (dc):
  
  ```
  saas.update_clv(mu,wav,clv,wav,dc)
  ```
  
  Update the solar differential rotation velocity profile:
  
  ```
  saas.update_vrot(0.,0.)
  ```
  
  the inputs for the function are the usual P0 and B0 angles.
  
  Finally, calculate the integrated spectrum by:
  
  ```
  test_si = saas.get_integration()
  ```

