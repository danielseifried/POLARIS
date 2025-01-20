/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CSOURCE_DUST_H
#define CSOURCE_DUST_H

#include "DustMixture.hpp"
#include "Vector3D.hpp"
#include "MathFunctions.hpp"
#include "Matrix2D.hpp"
#include "Parameters.hpp"
#include "SourceBasic.hpp"
#include "Stokes.hpp"
#include "Typedefs.hpp"

class CSourceDust : public CSourceBasic
{
public:
    CSourceDust(void)
    {
        total_energy = 0;
        cell_prob = 0;

        nr_of_photons_per_cell = 0;

        source_id = SRC_DUST;
    }

    ~CSourceDust(void)
    {
        if(cell_prob != 0)
            delete[] cell_prob;
        if(total_energy != 0)
            delete[] total_energy;
    }

    bool initSource(uint w);

    bool initSource(uint id, uint max, bool use_energy_density);

    void createNextRay(photon_package * pp, CRandomGenerator * rand_gen);
    void createDirectRay(photon_package * pp, CRandomGenerator * rand_gen, Vector3D dir_obs);

    void setParameter(parameters & param, uint p);

    ullong getNrOfPhotons();

private:
    double * total_energy;
    prob_list * cell_prob;
    ullong nr_of_photons_per_cell;
};

#endif /* CSOURCE_DUST_H */
