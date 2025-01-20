/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "SourceGas.hpp"

bool CSourceGas::initSource(uint id, uint max, bool use_energy_density)
{
    // Init variables
    ulong nr_of_cells = grid->getMaxDataCells();

    // Show information
    cout << "- Source GAS (lvl population):" << endl;
    cout << "    photons per cell: " << float(nr_of_photons) << " (" << nr_of_cells << " cells)" << endl;

    return true;
}

void CSourceGas::createNextRayToCell(photon_package * pp, CRandomGenerator * rand_gen, ulong i_cell, bool cell_as_border)
{
    pp->setRandomDirection(rand_gen->getRND(), rand_gen->getRND());

    // Put photon package into current cell
    pp->setPositionCell(grid->getCellFromIndex(i_cell));

    // Set random position in cell
    grid->setRndPositionInCell(pp, rand_gen);

    // Save final position as position of last interaction
    pp->setBackupPosition(pp->getPosition());

    if(cell_as_border)
    {
        // Go one cell outward
        grid->next(pp);

        // Invert direction
        pp->multDirection(-1);

        // Go into cell again
        grid->next(pp);
    }
    else
    {
        // Move photon along the path outwards the grid
        pp->addPosition(pp->getDirection() * grid->getMaxLength());

        // Invert direction
        pp->multDirection(-1);
    }

    // // Init coordinate System for polarization
    pp->initCoordSystem();
}

void CSourceGas::setParameter(parameters & param, uint p)
{
    nr_of_photons = param.getMCLvlPopNrOfPhotons();
}

ullong CSourceGas::getNrOfPhotons()
{
    return nr_of_photons;
}
