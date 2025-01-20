/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "MathInterp.hpp"

uint interp::size() const
{
    return N + 1;
}

void interp::resize(uint size)
{
    N = size - 1;
    x.resize(size);
    y.resize(size);
}

void interp::setValue(uint pos, double _x, double _y)
{
#ifdef DEBUG
    if(x == 0)
    {
        cout << ERROR_LINE << "Linear interpolation was not initiated!" << endl;
        return;
    }
#endif
    x[pos] = _x;
    y[pos] = _y;
}

void interp::addValue(double _x, double _y)
{
    x.push_back(_x);
    y.push_back(_y);
}

double interp::getLinear(uint i, double v) const
{
    double t = v - x[i];
    return y[i] + t * (y[i + 1] - y[i]) / (x[i + 1] - x[i]);
}

double interp::getValue(double v, uint interpolation) const
{
    if(N == 0)
        return y[0];

    if(v < x[0])
        switch(interpolation)
        {
            case CONST:
                return y[0];
                break;

            case LINEAR:
                return getLinear(0, v);
                break;
        }
    else if(v > x[N])
        switch(interpolation)
        {
            case CONST:
                return y[N];
                break;

            case LINEAR:
                return getLinear(N - 1, v);
                break;
        }
    else
    {
        uint min = 0;

        if(v != x[0])
            min = lower_bound(x.begin(), x.end(), v) - x.begin() - 1;

        switch(interpolation)
        {
            case CONST:
                return y[min+1];
                break;

            case LINEAR:
                return getLinear(min, v);
                break;
        }
    }

    return 0;
}
