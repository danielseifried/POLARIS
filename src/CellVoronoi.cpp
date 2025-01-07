#include "CellVoronoi.hpp"


void cell_vo::initNeighbors(short nr)
{
    nr_neighbors = nr;
    neighbors = new int[nr];

    for(ushort i = 0; i < nr; i++)
        neighbors[i] = 0;
}

void cell_vo::setNeighbor(uint pos, int id)
{
    if(pos > 500)
        return;

    neighbors[pos] = id;
}

void cell_vo::setCenter(double cx, double cy, double cz)
{
    center.set(cx, cy, cz);
}

void cell_vo::setVolume(double v)
{
    volume = v;
}

Vector3D cell_vo::getCenter() const
{
    return center;
}

double cell_vo::getX()
{
    return center.X();
}

double cell_vo::getY()
{
    return center.Y();
}

double cell_vo::getZ()
{
    return center.Z();
}

double cell_vo::getVolume() const
{
    return volume;
}

int cell_vo::getNeighborID(uint pos)
{
    return neighbors[pos];
}

ushort cell_vo::getNrOfNeighbors()
{
    return nr_neighbors;
}
