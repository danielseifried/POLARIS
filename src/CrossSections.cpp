#include "CrossSections.hpp"


void cross_sections::operator/=(double w)
{
    Cext /= w;
    Cpol /= w;
    Cabs /= w;
    Cpabs /= w;
    Ccirc /= w;
    Csca /= w;
}

void cross_sections::operator*=(double w)
{
    Cext *= w;
    Cpol *= w;
    Cabs *= w;
    Cpabs *= w;
    Ccirc *= w;
    Csca *= w;
}

cross_sections cross_sections::operator*(double w)
{
    Cext *= w;
    Cpol *= w;
    Cabs *= w;
    Cpabs *= w;
    Ccirc *= w;
    Csca *= w;

    return *this;
}

void cross_sections::operator+=(const cross_sections & cs)
{
    Cext += cs.Cext;
    Cpol += cs.Cpol;
    Cabs += cs.Cabs;
    Cpabs += cs.Cpabs;
    Ccirc += cs.Ccirc;
    Csca += cs.Csca;
}

void cross_sections::operator=(const double val)
{
    Cext = val;
    Cpol = val;
    Cabs = val;
    Cpabs = val;
    Ccirc = val;
    Csca = val;
}
