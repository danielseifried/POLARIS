/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "CellOcTree.hpp"

void cell_oc::setXmin(double _x_min)
{
    x_min = _x_min;
}

void cell_oc::setYmin(double _y_min)
{
    y_min = _y_min;
}

void cell_oc::setZmin(double _z_min)
{
    z_min = _z_min;
}

void cell_oc::setLength(double _length)
{
    length = _length;
}

void cell_oc::setLevel(uchar _level)
{
    level = _level;
}

void cell_oc::setParent(cell_oc * _parent)
{
    parent = _parent;
}

void cell_oc::setChildren(cell_oc * _children)
{
    children = _children;
}

void cell_oc::setUniqueID(ulong _unique_id)
{
    unique_id = _unique_id;
}

double cell_oc::getXmin() const
{
    return x_min;
}

double cell_oc::getYmin() const
{
    return y_min;
}

double cell_oc::getZmin() const
{
    return z_min;
}

double cell_oc::getXmax() const
{
    return x_min + length;
}

double cell_oc::getYmax() const
{
    return y_min + length;
}

double cell_oc::getZmax() const
{
    return z_min + length;
}

double cell_oc::getLength() const
{
    return length;
}

uchar cell_oc::getLevel()
{
    return level;
}

uint cell_oc::getID() const
{
    return id;
}

ulong cell_oc::getUniqueID() const
{
    return unique_id;
}

cell_oc * cell_oc::getParent()
{
    return parent;
}

cell_oc * cell_oc::getChildren()
{
    return children;
}

cell_oc * cell_oc::getChild(uint i)
{
    return &children[i];
}
