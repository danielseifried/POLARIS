/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#include "Photon.hpp"
#include "GridBasic.hpp"

void photon_package::setStokesVector(StokesVector st, uint i)
{
    if(i == MAX_UINT)
        multi_stokes[i_spectral] = st;
    else
        multi_stokes[i] = st;
}

const StokesVector & photon_package::getStokesVector(uint i) const
{
    if(i == MAX_UINT)
        return multi_stokes[i_spectral];
    else
        return multi_stokes[i];
}

StokesVector * photon_package::getStokesVector(uint i)
{
    if(i == MAX_UINT)
        return &multi_stokes[i_spectral];
    else
        return &multi_stokes[i];
}

StokesVector * photon_package::getMultiStokes()
{
    return multi_stokes;
}

void photon_package::setSpectralID(uint i)
{
    // Current index of the photon package (default wavelength or frequency)
    if(i < nr_of_spectral_bins)
        i_spectral = i;
}

uint photon_package::getSpectralID() const
{
    return i_spectral;
}

uint photon_package::getNrOfSpectralBins() const
{
    return nr_of_spectral_bins;
}

uint photon_package::getDustWavelengthID(uint i) const
{
    if(trans_frequency > 0)
        return wID[0];
    else if(i == MAX_UINT)
        return wID[i_spectral];
    else
        return wID[i];
}

void photon_package::setWavelength(double wave, uint _wID)
{
    wavelength[i_spectral] = wave;
    wID[i_spectral] = _wID;
}

double photon_package::getWavelength(uint i) const
{
    if(trans_frequency > 0)
        return wavelength[0];
    else if(i == MAX_UINT)
        return wavelength[i_spectral];
    else
        return wavelength[i];
}

double photon_package::getFrequency() const
{
    if(trans_frequency > 0)
        return getTransFrequency() + getVelocity() / con_c * getTransFrequency();
    return con_c / getWavelength();
}

void photon_package::setVelocity(double vel)
{
    velocity[i_spectral] = vel;
}

double photon_package::getVelocity(uint i) const
{
    if(i == MAX_UINT)
        return velocity[i_spectral];
    else
        return velocity[i];
}

double * photon_package::getVelocities() const
{
    return velocity;
}

double photon_package::getTransFrequency() const
{
    return trans_frequency;
}

void photon_package::setD(Matrix2D _mD)
{
    mD = _mD;
}

const Matrix2D & photon_package::getD() const
{
    return mD;
}

const Vector3D & photon_package::getPosition() const
{
    return pos;
}

const Vector3D & photon_package::getDirection() const
{
    return ez;
}

double photon_package::getTmpPathLength() const
{
    return tmp_path;
}

cell_basic * photon_package::getPositionCell()
{
    return cell_pos;
}

const cell_basic * photon_package::getPositionCell() const
{
    return cell_pos;
}

void photon_package::setRandomDirection(double r1, double r2, uint exponentThetaBias)
{
    ez.rndDir(r1, r2, exponentThetaBias);
}

void photon_package::setRandomDirectionTRUST(double r1, double r2)
{
    ez.rndDirTRUST(r1, r2);
}

void photon_package::initCoordSystem()
{
    double phi = Vector3D::atan3(ez.Y(), -ez.X());
    double theta = acos(ez.Z());
    double cos_phi = cos(phi);
    double sin_phi = sin(phi);
    double cos_theta = ez.Z();
    double sin_theta = sin(theta);

    mD = CMathFunctions::getRotationMatrix(cos_phi, sin_phi, cos_theta, sin_theta);

    ex = mD * Vector3D(1, 0, 0);
    ey = mD * Vector3D(0, 1, 0);
}

void photon_package::updateCoordSystem()
{
    double phi = Vector3D::atan3(ez.Y(), -ez.X());
    double theta = acos(ez.Z());
    double cos_phi = cos(phi);
    double sin_phi = sin(phi);
    double cos_theta = ez.Z();
    double sin_theta = sin(theta);

    Matrix2D D_help = CMathFunctions::getRotationMatrix(cos_phi, sin_phi, cos_theta, sin_theta);
    mD = mD * D_help;

    ex = mD * Vector3D(1, 0, 0);
    ey = mD * Vector3D(0, 1, 0);
    ez = mD * Vector3D(0, 0, 1);
}

void photon_package::updateCoordSystem(double phi, double theta)
{
    double cos_phi = cos(phi);
    double sin_phi = sin(phi);
    double cos_theta = cos(theta);
    double sin_theta = sin(theta);

    Matrix2D D_help = CMathFunctions::getRotationMatrix(cos_phi, sin_phi, cos_theta, sin_theta);
    mD = mD * D_help;

    ex = mD * Vector3D(1, 0, 0);
    ey = mD * Vector3D(0, 1, 0);
    ez = mD * Vector3D(0, 0, 1);
}

void photon_package::adjustPosition(Vector3D _pos, double _len)
{
    tmp_path = _len;
    pos = _pos + tmp_path * ez;

    // Calculate cell by position
    if(tmp_path > 0)
        dirID = MAX_UINT;
}

void photon_package::setPosition(Vector3D val)
{
    pos = val;
}

void photon_package::addPosition(Vector3D val)
{
    pos += val;
}

void photon_package::setBackupPosition(Vector3D val)
{
    backup_pos = val;
}

void photon_package::updateBackupPosition()
{
    backup_pos = pos;
}

void photon_package::resetPositionToBackupPos()
{
    pos = backup_pos;
}

bool photon_package::reachedBackupPosition()
{
    if(ez * (backup_pos - pos) <= 0)
        return true;
    return false;
}

bool photon_package::reachedBackupPosition(Vector3D tmp_pos)
{
    if(ez * (backup_pos - tmp_pos) <= 0)
        return true;
    return false;
}

const Vector3D & photon_package::getBackupPosition() const
{
    return backup_pos;
}

void photon_package::setDetectorProjection()
{
    double tmp_vec = ey * pos;
    pos.setX(ex * pos);
    pos.setY(tmp_vec);
}

void photon_package::setRelativePosition(double tx, double ty, double tz)
{
    pos = tx * ex + ty * ey + tz * ez;
}

void photon_package::setDirection(Vector3D val)
{
    ez = val;
}

void photon_package::addDirection(Vector3D val)
{
    ez += val;
}

void photon_package::multDirection(double val)
{
    ez *= val;
}

const Vector3D & photon_package::getEX() const
{
    return ex;
}

const Vector3D & photon_package::getEY() const
{
    return ey;
}

const Vector3D & photon_package::getEZ() const
{
    return ez;
}

void photon_package::setEX(Vector3D _e)
{
    ex = _e;
}

void photon_package::setEY(Vector3D _e)
{
    ey = _e;
}

void photon_package::setEZ(Vector3D _e)
{
    ez = _e;
}

void photon_package::setCoordinateSystem(const Vector3D & _ex, const Vector3D & _ey, const Vector3D & _ez)
{
    ex = _ex;
    ey = _ey;
    ez = _ez;
}

void photon_package::getCoordinateSystem(Vector3D * _ex, Vector3D * _ey, Vector3D * _ez) const
{
    *_ex = ex;
    *_ey = ey;
    *_ez = ez;
}

void photon_package::setTmpPathLength(double val)
{
    tmp_path = val;
}

void photon_package::setPositionCell(cell_basic * val)
{
    cell_pos = val;
}

void photon_package::setDirectionID(uint val)
{
    dirID = val;
}

uint photon_package::getDirectionID()
{
    return dirID;
}

void photon_package::setPhotonID(ullong _photonID)
{
    photonID = _photonID;
}

ullong photon_package::getPhotonID()
{
    return photonID;
}
