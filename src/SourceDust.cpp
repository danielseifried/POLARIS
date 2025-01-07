#include "SourceDust.hpp"
#include "SourceBasic.hpp"


bool CSourceDust::initSource(uint id, uint max, bool use_energy_density)
{
    if(!use_energy_density)
    {
        cout << "\nERROR: The dust source for radiation field calculation "
             << "can only be used with energy density!" << endl;
        return false;
    }

    // Init variables
    ulong nr_of_cells = grid->getMaxDataCells();
    ulong nr_of_wavelengths = getNrOfWavelength();
    photon_package pp = photon_package();
    cell_prob = new prob_list[nr_of_wavelengths];
    total_energy = new double[nr_of_wavelengths];
    uint max_counter = nr_of_wavelengths * nr_of_cells;
    ullong per_counter = 0;
    float last_percentage = 0;

    // Show Initial message
    cout << CLR_LINE;
    cout << "-> Initiating dust grain emission          \r" << flush;

    for(uint w = 0; w < nr_of_wavelengths; w++)
    {
        // Init variables
        cell_prob[w].resize(nr_of_cells + 1);

        // Set wavelength of photon package
        pp.setWavelength(wavelength_list[w], w);

        // Set total energy to zero and starting value of prob_list
        total_energy[w] = 0;
        cell_prob[w].setValue(0, total_energy[w]);

        #if (DUST_EMI_PROB)
            // Causes problems. Find better solution!
            //#pragma omp parallel for schedule(dynamic)
            for(long i_cell = 0; i_cell < long(nr_of_cells); i_cell++)
            {
                // Increase counter used to show progress
                per_counter++;

                // Calculate percentage of total progress
                float percentage = 100.0 * float(per_counter) / float(max_counter);

                // Show only new percentage number if it changed
                if((percentage - last_percentage) > PERCENTAGE_STEP)
                {
                    //#pragma omp critical
                    {
                        cout << "-> Calculating prob. distribution for dust source: " << percentage << " [%]    \r"
                            << flush;
                        last_percentage = percentage;
                    }
                }

                // Put photon package into current cell
                pp.setPositionCell(grid->getCellFromIndex(i_cell));

                // Get total energy of thermal emission
                total_energy[w] += dust->getTotalCellEmission(grid, pp);

                // Add energy to probability distribution
                cell_prob[w].setValue(i_cell + 1, total_energy[w]);
            }

            // Normalize probability distribution
            cell_prob[w].normalize(total_energy[w]);
        #endif
    }

    #if (!DUST_EMI_PROB)
        // reduce the total number of photons, so that each cell launches the same amount of photons,
        // i.e. (n_photon % n_cell) should be zero
        nr_of_photons -= nr_of_photons % grid->getMaxDataCells();
        nr_of_photons_per_cell = ullong(nr_of_photons / double(nr_of_cells));
    #endif

    // Show information
    cout << CLR_LINE;
    cout << "- Source (" << id + 1 << " of " << max
         << ") DUST: photons FOR EACH wavelength: " << float(nr_of_photons) << "      " << endl;

    return true;
}

bool CSourceDust::initSource(uint w)
{
    // Init variables
    ulong nr_of_cells = grid->getMaxDataCells();
    photon_package pp = photon_package();
    uint max_counter = nr_of_cells;
    ullong per_counter = 0;
    float last_percentage = 0;

    // Init variables
    total_energy = new double[getNrOfWavelength()];
    cell_prob = new prob_list[getNrOfWavelength()];
    cell_prob[w].resize(nr_of_cells + 1);

    // Set wavelength of photon package
    pp.setWavelength(wavelength_list[w], w);

    // Set total energy to zero and starting value of prob_list
    total_energy[w] = 0;
    cell_prob[w].setValue(0, total_energy[w]);

    // Show Initial message
    cout << "-> Initiating dust grain emission          \r" << flush;

    #if (DUST_EMI_PROB)
        // Causes problems. Find better solution!
        //#pragma omp parallel for schedule(dynamic)
        for(long i_cell = 0; i_cell < long(nr_of_cells); i_cell++)
        {
            // Increase counter used to show progress
            //#pragma omp atomic update
            per_counter++;

            // Calculate percentage of total progress
            float percentage = 100.0 * float(per_counter) / float(max_counter);

            // Show only new percentage number if it changed
            if((percentage - last_percentage) > PERCENTAGE_STEP)
            {
                //#pragma omp critical
                {
                    cout << "-> Calculate prob. distribution for dust source: " << percentage << " [%]    \r"
                        << flush;
                    last_percentage = percentage;
                }
            }

            // Put photon package into current cell
            pp.setPositionCell(grid->getCellFromIndex(i_cell));

            // Get total energy of thermal emission
            total_energy[w] += dust->getTotalCellEmission(grid, pp);

            // Add energy to probability distribution
            cell_prob[w].setValue(i_cell + 1, total_energy[w]);
        }

        // Normalize probability distribution
        cell_prob[w].normalize(total_energy[w]);
    #else
        // reduce the total number of photons, so that each cell launches the same amount of photons,
        // i.e. (n_photon % n_cell) should be zero
        nr_of_photons -= nr_of_photons % grid->getMaxDataCells();
        nr_of_photons_per_cell = ullong(nr_of_photons / double(nr_of_cells));
    #endif

    return true;
}

void CSourceDust::createNextRay(photon_package * pp, CRandomGenerator * rand_gen)
{
    pp->setRandomDirection(rand_gen->getRND(), rand_gen->getRND());

    // Set wavelength of photon package
    uint w = pp->getDustWavelengthID();

    // Get index of current cell
    ulong i_cell = 0;
    #if (DUST_EMI_PROB)
        // Get random number
        double rnd = rand_gen->getRND();

        i_cell = cell_prob[w].getIndex(rnd);
    #else
        i_cell = ulong(pp->getPhotonID() % grid->getMaxDataCells());
    #endif

    // Put photon package into current cell
    pp->setPositionCell(grid->getCellFromIndex(i_cell));

    // Set random position in cell
    grid->setRndPositionInCell(pp, rand_gen);

    // Set Stokes vector of photon package
    double energy = 0;
    #if (DUST_EMI_PROB)
        energy = total_energy[w] / double(nr_of_photons);
    #else
        energy = dust->getTotalCellEmission(grid, *pp) / double(nr_of_photons_per_cell);
    #endif

    // Set Stokes Vector
    pp->setStokesVector(StokesVector(energy, 0, 0, 0));

    // Init coordinate System for polarization
    pp->initCoordSystem();
}

void CSourceDust::createDirectRay(photon_package * pp, CRandomGenerator * rand_gen, Vector3D dir_obs)
{
    // Set wavelength of photon package
    uint w = pp->getDustWavelengthID();

    // Get random number
    double rnd = rand_gen->getRND();

    // Get index of current cell
    ulong i_cell = cell_prob[w].getIndex(rnd);

    // Put photon package into current cell
    pp->setPositionCell(grid->getCellFromIndex(i_cell));

    // Set random position in cell
    grid->setRndPositionInCell(pp, rand_gen);

    // Set Stokes vector of photon package
    double energy = total_energy[w] / PIx4;

    // Set Stokes Vector
    pp->setStokesVector(StokesVector(energy, 0, 0, 0));

    // Set direction of the photon package to the observer
    if(dir_obs.length() > 0)
    {
        pp->setDirection(dir_obs);
        pp->initCoordSystem();
    }
}

void CSourceDust::setParameter(parameters & param, uint p)
{
    nr_of_photons = param.getNrOfDustPhotons();
}

ullong CSourceDust::getNrOfPhotons()
{
    return nr_of_photons;
}
