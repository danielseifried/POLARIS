/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "CellCylindrical.hpp"

void cell_cyl::setRID(uint id)
{
    rID = id;
}

void cell_cyl::setPhID(uint id)
{
    phID = id;
}

void cell_cyl::setZID(uint id)
{
    zID = id;
}

uint cell_cyl::getRID() const
{
    return rID;
}

uint cell_cyl::getPhID() const
{
    return phID;
}

uint cell_cyl::getZID() const
{
    return zID;
}
