/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CSOURCE_GAS_H
#define CSOURCE_GAS_H

#include "DustMixture.hpp"
#include "Vector3D.hpp"
#include "MathFunctions.hpp"
#include "Matrix2D.hpp"
#include "Parameters.hpp"
#include "SourceBasic.hpp"
#include "Stokes.hpp"
#include "Typedefs.hpp"

class CSourceGas : public CSourceBasic
{
public:
    CSourceGas(void)
    {
        source_id = SRC_GAS_LVL;
    }

    bool initSource(uint id, uint max, bool use_energy_density);
    void createNextRayToCell(photon_package * pp, CRandomGenerator * rand_gen, ulong i_cell, bool cell_as_border);

    void setParameter(parameters & param, uint p);

    ullong getNrOfPhotons();
};

#endif /* CSOURCE_GAS_H */
