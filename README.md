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

## Uninstall NESSI:

  ```
  python3 -m pip uninstall nessi
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
  from nessi import integrator as nss
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
  
  mu and wav are 1d arrays with the mu, wav values for which the clv is evaluated.
  Let set nmu to the size of mu and nwav the size of wav
  clv is a 2d array of shape: (nmu, nwav) with the CLV to be used. If the class initialization includes any of include\_red\_cclv or include\_full\_cclv then clv must only include the CLV as respect to the disk center profile and without any continuum CLV (as it is included from Neckel & Labs 1992).
  Update the solar differential rotation velocity profile:
  
  ```
  saas.update_vrot(0.,0.)
  ```
  
  the inputs for the function are the usual P0 and B0 angles.
  
  Finally, calculate the integrated spectrum by:
  
  ```
  test_si = saas.get_integration()
  ```

  If you want to obtain the spectrum for each grid point employed in the integration you can call:

  ```
  full_spectrum = saas.get_spectra()
  ```

  that returns a 2d array with dimension (nwav, npxl), where npxl are the number of grid points employed.


## Minimal example for FOV substitution:

  Since we need to calculate the Sun as a Star spectra first, we can proceed as in the "Minimal example":

  ```
  from nessi.tester import load_data
  
  wav, dc, mu, clv = load_data()
  from nessi import integrator as nss
  saas = nss.sun_as_a_star(nr=101)
  saas.update_clv(mu,wav,clv,wav,dc)
  saas.update_vrot(0.,0.)
  test_si = saas.get_integration()
  ```

  If you want to evaluate the impact a certain FOV has over the Sun as a Star spectra, there is a function called: "get_diff_spectra_fov", which actually calculates the differential impact that a given spectra for a FOV has on the Sun as a Star spectra.

  To do so, you need to provide not only the spectra (nw,ny,nx; with nw the number of wavelengths in test_si -to be updated at some point-, nx and ny the number of pixels in each spatial dimension) but also the coordinates of the vertices of the FOV pixels (ny+1,nx+1):

  ```
  from nessi.tester import load_data_fov

  xx, yy, fov_spectra = load_data_fov(test_si.size)
  print(xx.shape)
  print(yy.shape)
  print(fov_spectra.shape)
  ```
  Once with these three arrays, the function returns the difference between the provided spectra and the one that the code would use in the whole Sun integration. This differential spectrum is weighted for the whole visible solar disk so if you want to have the Sun as a Star spectrum taking into account the spectra for the FOV you are considering ("fov_si" below) you just need to add the reference Sun as a Star spectrum ("test_si") and the weighted difference ("diff").
  ```
  diff = saas.get_diff_spectra_fov(xx,yy,fov_spectra)

  fov_si = test_si + diff
  ```
  

