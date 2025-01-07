#pragma once

#include "Matrix2D.hpp"
#include "Stokes.hpp"
#include "Typedefs.hpp"


#ifndef CMATH_INTERP_H
#define CMATH_INTERP_H

class interp
{
  public:
    interp()
    {
        N = 0;
    }

    interp(uint size)
    {
        N = size - 1;
        x.resize(size);
        y.resize(size);
    }

    uint size() const;

    void resize(uint size);

    void setValue(uint pos, double _x, double _y);

    void addValue(double _x, double _y);

    double getLinear(uint i, double v) const;

    double getValue(double v, uint interpolation = LINEAR) const;

  private:
    uint N;
    dlist x;
    dlist y;
};

#endif
