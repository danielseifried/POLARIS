/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef CELL_CYLINDRICAL_H
#define CELL_CYLINDRICAL_H

#include "CellBasic.hpp"

class cell_cyl : public cell_basic
{
public:
    cell_cyl()
    {
        rID = 0;
        phID = 0;
        zID = 0;
        data = 0;
        id = 0;
    }

    ~cell_cyl()
    {
        if(data != 0)
            delete[] data;

        data = 0;
    }

    void setRID(uint id);

    void setPhID(uint id);

    void setZID(uint id);

    uint getRID() const;

    uint getPhID() const;

    uint getZID() const;

private:
    uint rID, phID, zID;
};

#endif /* CELL_CYLINDRICAL_H */
