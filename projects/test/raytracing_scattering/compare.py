from astropy.io import fits
from astropy import units as u
import numpy as np


"""
----------------------------------------------
Raytracing and scattering 
----------------------------------------------

Test whether scattering in a dust simulation
(from radiation field) looks the same as
scattering in a dust_mc simulation.

In order to pass this test, all detectors
(with different orientations) have to yield
the same stellar + scattered + emitted flux 
for each wavelength.
(Should get better with more photons)
"""


def read_data(sed_fits_file):
    fits_header = fits.getheader(sed_fits_file)
    fits_data = fits.getdata(sed_fits_file)

    nr_wave = np.shape(fits_data)[2]
    sed_wavelengths = np.zeros(nr_wave) * u.m
    for i_wave in range(nr_wave):
        sed_wavelengths[i_wave] = float(fits_header[f'HIERARCH WAVELENGTH{i_wave+1}']) * u.m

    _stokes = ['I', 'Q', 'U', 'V']
    sed_data = {}
    for i_s, i_stokes in enumerate(_stokes):
        sed_data[i_stokes] = fits_data[i_s,0,:] * u.Jy

    return sed_data


def compare():
    dust_sed_data = read_data('projects/test/raytracing_scattering/dust/data/polaris_detector_nr0001_sed.fits.gz')
    dust_rt_sed_data = read_data('projects/test/raytracing_scattering/dust_rt/data/polaris_detector_nr0001_sed.fits.gz')
    dust_mc_sed_data = read_data('projects/test/raytracing_scattering/dust_mc/data/polaris_detector_nr0001_sed.fits.gz')

    max_rel_diff = np.max(np.abs( (dust_sed_data['I'] + dust_mc_sed_data['I']) / dust_rt_sed_data['I'] - 1.0 ))
    if max_rel_diff > 1e-3:
        raise Exception(f'Test failed: Stokes I does not match (max. relative difference = {max_rel_diff})')

    mc_polarization = np.sqrt(dust_mc_sed_data['Q']**2 + dust_mc_sed_data['U']**2 + dust_mc_sed_data['V']**2) / (dust_mc_sed_data['I'] + dust_sed_data['I'])
    rt_polarization = np.sqrt(dust_rt_sed_data['Q']**2 + dust_rt_sed_data['U']**2 + dust_rt_sed_data['V']**2) / dust_rt_sed_data['I']
    max_abs_diff = np.max(np.abs( mc_polarization - rt_polarization ))
    if max_abs_diff > 1e-3:
        raise Exception(f'Test failed: Polarization does not match (max. absolute difference = {max_abs_diff})')

    max_polarization = np.max(mc_polarization)
    if max_polarization > 1e-3:
        raise Exception(f'Test failed: Polarization is too large (max. value = {max_polarization})')

    return True


if __name__ == '__main__':
    res = compare()
    if res:
        print('Test passed')
