#pragma once
#include "Typedefs.hpp"
#include "Vector3D.hpp"


#ifndef CELL_BASIC_H
#define CELL_BASIC_H

class cell_basic
{
  public:
    cell_basic()
    {
        data = 0;
    }

    virtual ~cell_basic()
    {
        if(data != 0)
            delete[] data;

        data = 0;
    }

    bool isValid();

    double getData(uint i) const;

    void setData(uint i, double d);

    void updateData(uint i, double d);

    void convertData(uint i, double c);

    void resize(uint size);

    void setID(uint _id);

    virtual ulong getUniqueID() const;

    void updateID(uint _id);

  protected:
    double * data;
    uint id;
};

#endif
