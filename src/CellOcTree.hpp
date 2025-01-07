#pragma once
#include "CellBasic.hpp"


#ifndef CELL_OCTREE_H
#define CELL_OCTREE_H

class cell_oc : public cell_basic
{
  public:
    cell_oc()
    {
        parent = 0;
        children = 0;
        level = 0;
        id = 0;
        unique_id = 0;
        x_min = 0;
        y_min = 0;
        z_min = 0;
        length = 0;
        data = 0;
    }

    ~cell_oc()
    {
        if(data != 0)
            delete[] data;

        data = 0;
    }

    void setXmin(double _x_min);

    void setYmin(double _y_min);

    void setZmin(double _z_min);

    void setLength(double _length);

    void setLevel(uchar _level);

    void setParent(cell_oc * _parent);

    void setChildren(cell_oc * _children);

    void setUniqueID(ulong _unique_id);

    double getXmin() const;

    double getYmin() const;

    double getZmin() const;

    double getXmax() const;

    double getYmax() const;

    double getZmax() const;

    double getLength() const;

    uchar getLevel();

    /*
    get the position ID within the parent cube (0-7)
    */
    uint getID() const;

    /*
    get the actual ID of the cell (globally unique)
    */
    ulong getUniqueID() const;

    cell_oc * getParent();

    cell_oc * getChildren();

    cell_oc * getChild(uint i);

  private:
    double x_min, y_min, z_min, length;
    cell_oc * children;
    cell_oc * parent;
    uchar level;
    ulong unique_id;
};

#endif
