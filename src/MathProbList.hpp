/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CMATH_PROBLIST_H
#define CMATH_PROBLIST_H

#include "Matrix2D.hpp"
#include "Stokes.hpp"
#include "Typedefs.hpp"

class prob_list
{
public:
    prob_list()
    {
        N = 0;
        x = 0;
    }

    prob_list(uint size)
    {
        N = size - 1;
        x = new double[size];

        for(uint i = 0; i < size; i++)
            x[i] = 0;
    }

    ~prob_list()
    {
        if(x != 0)
            delete[] x;
    }

    uint size();

    double X(uint i);

    void resize(uint size);

    void setValue(uint pos, double _x);

    uint getIndex(double v);

    void normalize(double integValue);

    friend prob_list operator-(prob_list & list1, prob_list & list2);

private:
    uint N;
    double * x;
};

#endif /* CMATH_PROBLIST_H */
