#include "SourceLaser.hpp"


bool CSourceLaser::initSource(uint id, uint max, bool use_energy_density)
{
    // Initial output
    cout << CLR_LINE << flush;
    cout << "-> Initiating source laser         \r" << flush;

    if(use_energy_density)
    {
        // For using energy density, only the photon number is required
        cout << "- Source (" << id + 1 << " of " << max << ") LASER: " << float(L)
             << " [W] (wavelength = " << float(wl) << " [m], FWHM = " << float(fwhm) << " [m])" << endl;
        cout << "    photons FOR EACH wavelength: " << float(nr_of_photons) << "      " << endl;
        // cout << "\nERROR: Laser source cannot be used for dust temperature calculation!\n" << endl;
        // return false;
    }
    else
    {
        // Init variables
        dlist star_emi;
        double tmp_luminosity;

        for(uint w = 0; w < getNrOfWavelength(); w++)
        {
            double line_shape =
                1 / sqrt(2 * PI * sigma_sq) * exp(-pow(wavelength_list[w] - wl, 2) / (2 * sigma_sq));
            star_emi.push_back(L * line_shape);
        }

        tmp_luminosity = CMathFunctions::integ(wavelength_list, star_emi, 0, getNrOfWavelength() - 1);
        // L = tmp_luminosity;

        // For using energy density, only the photon number is required
        cout << "- Source (" << id + 1 << " of " << max << ") LASER: " << float(L)
             << " [W] (wavelength = " << float(wl) << " [m], FWHM = " << float(fwhm) << " [m])" << endl;
        cout << "    photons: " << float(nr_of_photons) << "      " << endl;

        double fr;
        lam_pf.resize(getNrOfWavelength());
        for(uint l = 0; l < getNrOfWavelength(); l++)
        {
            fr = CMathFunctions::integ(wavelength_list, star_emi, 0, l) / tmp_luminosity;
            lam_pf.setValue(l, fr, double(l));
        }
    }

    return true;
}

void CSourceLaser::createNextRay(photon_package * pp, CRandomGenerator * rand_gen)
{
    // Init variables
    StokesVector tmp_stokes_vector;
    uint wID = 0;

    pp->setDirection(dir);
    pp->setPosition(pos);

    if(pp->getDustWavelengthID() != MAX_UINT)
        wID = pp->getDustWavelengthID();
    else
    {
        wID = lam_pf.getXIndex(rand_gen->getRND());
        // Mol3D uses the upper value of the wavelength interval,
        // used for the selection of the emitting wavelengths from source!
        pp->setWavelength(wavelength_list[wID + 1], wID + 1);
    }

    double line_shape =
        1 / sqrt(2 * PI * sigma_sq) * exp(-pow(wavelength_list[wID] - wl, 2) / (2 * sigma_sq));
    tmp_stokes_vector = L / double(nr_of_photons) * line_shape * StokesVector(1.0, q, u, 0);

    pp->setStokesVector(tmp_stokes_vector);
    pp->initCoordSystem();
}

void CSourceLaser::createDirectRay(photon_package * pp, CRandomGenerator * rand_gen, Vector3D dir_obs)
{
    // Init variables
    StokesVector tmp_stokes_vector;

    pp->setDirection(dir);
    pp->setPosition(pos);

    double dir_err = 1e-10;
    if(dir_obs * dir > 1. - dir_err && pp->getDustWavelengthID() != MAX_UINT)
    {
        uint wID = pp->getDustWavelengthID();
        double line_shape =
            1 / sqrt(2 * PI * sigma_sq) * exp(-pow(wavelength_list[wID] - wl, 2) / (2 * sigma_sq));
        tmp_stokes_vector = L * line_shape * StokesVector(1.0, q, u, 0);
    }

    pp->setStokesVector(tmp_stokes_vector);

    // Set direction of the photon package to the observer
    if(dir_obs.length() > 0)
    {
        pp->initCoordSystem();
    }
}

void CSourceLaser::setParameter(parameters & param, uint p)
{
    dlist values = param.getLaserSources();

    pos = Vector3D(values[p], values[p + 1], values[p + 2]);
    dir = Vector3D(values[p + 3], values[p + 4], values[p + 5]);
    dir.normalize();

    L = values[p + 6];
    wl = values[p + 7];
    fwhm = values[p + 8];
    q = values[p + 9];
    u = values[p + 10];

    nr_of_photons = (ullong)values[p + NR_OF_LASER_SOURCES - 1];

    // FWHM to sigma^2
    sigma_sq = fwhm * fwhm / (8.0 * log(2.0));
}
