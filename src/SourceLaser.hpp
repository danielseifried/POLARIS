#pragma once

#include "DustMixture.hpp"
#include "Vector3D.hpp"
#include "MathFunctions.hpp"
#include "Matrix2D.hpp"
#include "Parameters.hpp"
#include "SourceBasic.hpp"
#include "Stokes.hpp"
#include "Typedefs.hpp"


class CGridBasic;
class photon_package;

#ifndef CSOURCE_LASER_H
#define CSOURCE_LASER_H


class CSourceLaser : public CSourceBasic
{
  public:
    CSourceLaser()
    {
        source_id = SRC_LASER;
    }

    ~CSourceLaser()
    {}

    bool initSource(uint id, uint max, bool use_energy_density);

    void createNextRay(photon_package * pp, CRandomGenerator * rand_gen);
    void createDirectRay(photon_package * pp, CRandomGenerator * rand_gen, Vector3D dir_obs);

    void setParameter(parameters & param, uint p);

  protected:
    Vector3D dir;
    double wl, sigma_sq, fwhm;
};

#endif