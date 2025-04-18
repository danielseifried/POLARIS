/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CROSS_SECTIONS_H
#define CROSS_SECTIONS_H

#include "Typedefs.hpp"

class cross_sections
{
public:
    cross_sections()
    {
        Cext = 0;
        Cpol = 0;
        Cabs = 0;
        Cpabs = 0;
        Ccirc = 0;
        Csca = 0;
    }

    void operator/=(double w);

    void operator*=(double w);

    cross_sections operator*(double w);

    void operator+=(const cross_sections & cs);

    void operator=(const double val);

    double Cext;
    double Cpol;
    double Cabs;
    double Cpabs;
    double Csca;
    double Ccirc;
};

#endif /* CROSS_SECTIONS_H */
