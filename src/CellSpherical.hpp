/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CELL_SPHERICAL_H
#define CELL_SPHERICAL_H

#include "CellBasic.hpp"

class cell_sp : public cell_basic
{
public:
    cell_sp()
    {
        rID = 0;
        phID = 0;
        thID = 0;
        data = 0;
        id = 0;
    }

    cell_sp(uint _rID, uint _phID, uint _thID)
    {
        rID = _rID;
        phID = _phID;
        thID = _thID;
        data = 0;
        id = 0;
    }

    ~cell_sp()
    {
        if(data != 0)
            delete[] data;

        data = 0;
    }

    void setRID(uint id);

    void setPhID(uint id);

    void setThID(uint id);

    uint getRID() const;

    uint getPhID() const;

    uint getThID() const;

private:
    uint rID, phID, thID;
};

#endif /* CELL_SPHERICAL_H */
