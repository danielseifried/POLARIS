#pragma once

#include "Matrix2D.hpp"
#include "Stokes.hpp"
#include "Typedefs.hpp"


#ifndef RANDOM_GENERATOR_H
#define RANDOM_GENERATOR_H

class CRandomGenerator
{
  public:
    CRandomGenerator()
    {
        // The standard seed values as proposed by George Marsaglia
        // These will be used is init is not called
        // x
        KISS_state[0] = 1234567890987654321ULL;
        // y
        KISS_state[1] = 362436362436362436ULL;
        // z
        KISS_state[2] = 1066149217761810ULL;
        // c
        KISS_state[3] = 123456123456123456ULL;
    }

    void init(ullong seed);

    ullong CONG(ullong current_state);

    double getRND();

    double getRNDnormal(double mu, double sigma);

  private:
    ullong KISS_state[4];
};

#endif
