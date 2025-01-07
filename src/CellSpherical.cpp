#include "CellSpherical.hpp"


void cell_sp::setRID(uint id)
{
    rID = id;
}

void cell_sp::setPhID(uint id)
{
    phID = id;
}

void cell_sp::setThID(uint id)
{
    thID = id;
}

uint cell_sp::getRID() const
{
    return rID;
}

uint cell_sp::getPhID() const
{
    return phID;
}

uint cell_sp::getThID() const
{
    return thID;
}
