#pragma once

#include "CellBasic.hpp"
#include "MathFunctions.hpp"
#include "Stokes.hpp"
#include "Typedefs.hpp"
#include "Vector3D.hpp"


#ifndef PHOTON_H
#define PHOTON_H

class photon_package
{
  public:
    photon_package(uint nr_bins = 1)
    {
        // Photon package for dust simulations

        // Init variables
        cell_pos = 0;
        tmp_path = 0;
        i_spectral = 0;
        velocity = 0;
        trans_frequency = 0;
        dirID = MAX_UINT;

        // Set total number of wavelength saved in photon package (usually 1)
        nr_of_spectral_bins = nr_bins;

        // Init pointer arrays
        wavelength = new double[nr_of_spectral_bins];
        wID = new uint[nr_of_spectral_bins];
        multi_stokes = new StokesVector[nr_of_spectral_bins];

        // Set initial values
        for(uint i_spectral = 0; i_spectral < nr_of_spectral_bins; i_spectral++)
        {
            wavelength[i_spectral] = 0;

            // Dust wavelength index has to be MAX_UINT if not set in routine
            wID[i_spectral] = MAX_UINT;
        }
    }

    photon_package(double trans_freq, uint _wID, uint nr_bins = 1)
    {
        // Photon package for spectral line simulations

        // Init variables
        cell_pos = 0;
        tmp_path = 0;
        i_spectral = 0;
        dirID = MAX_UINT;

        // Set number of velocity channels saved in photon package
        nr_of_spectral_bins = nr_bins;

        // Set rest frequency of transition (used for dust properties)
        trans_frequency = trans_freq;

        // Init pointer arrays
        velocity = new double[nr_of_spectral_bins];
        multi_stokes = new StokesVector[nr_of_spectral_bins];

        // Set initial values
        for(uint i_spectral = 0; i_spectral < nr_of_spectral_bins; i_spectral++)
        {
            velocity[i_spectral] = 0;
        }

        // Set wavelength for dust contribution
        wavelength = new double[1];
        wavelength[0] = con_c / getTransFrequency();

        // Set wavelength index to get the right dust properties
        wID = new uint[1];
        wID[0] = _wID;
    }

    ~photon_package()
    {
        if(multi_stokes != 0)
            delete[] multi_stokes;

        if(wID != 0)
            delete[] wID;

        if(wavelength != 0)
            delete[] wavelength;

        if(velocity != 0)
            delete[] velocity;
    }

    void setStokesVector(StokesVector st, uint i = MAX_UINT);
    const StokesVector & getStokesVector(uint i = MAX_UINT) const;

    StokesVector * getStokesVector(uint i = MAX_UINT);

    StokesVector * getMultiStokes();

    void setSpectralID(uint i);

    uint getSpectralID() const;

    uint getNrOfSpectralBins() const;

    uint getDustWavelengthID(uint i = MAX_UINT) const;

    void setWavelength(double wave, uint _wID);

    double getWavelength(uint i = MAX_UINT) const;

    double getFrequency() const;

    void setVelocity(double vel);

    double getVelocity(uint i = MAX_UINT) const;

    double * getVelocities() const;

    double getTransFrequency() const;

    void setD(Matrix2D _mD);

    const Matrix2D & getD() const;

    const Vector3D & getPosition() const;

    const Vector3D & getDirection() const;

    double getTmpPathLength() const;

    cell_basic * getPositionCell();

    const cell_basic * getPositionCell() const;

    void setRandomDirection(double r1, double r2, uint exponentThetaBias = 1);

    void setRandomDirectionTRUST(double r1, double r2);

    void initCoordSystem();

    void updateCoordSystem();

    void updateCoordSystem(double phi, double theta);

    void adjustPosition(Vector3D _pos, double _len);
    void setPosition(Vector3D val);

    void addPosition(Vector3D val);

    void setBackupPosition(Vector3D val);

    void updateBackupPosition();

    void resetPositionToBackupPos();

    bool reachedBackupPosition();

    bool reachedBackupPosition(Vector3D tmp_pos);

    const Vector3D & getBackupPosition() const;

    void setDetectorProjection();

    void setRelativePosition(double tx, double ty, double tz);

    void setDirection(Vector3D val);

    void addDirection(Vector3D val);

    void multDirection(double val);

    /*
    get r-axis, based on O. Fischer (1993)
    */
    const Vector3D & getEX() const;

    /*
    get l-axis, based on O. Fischer (1993)
    */
    const Vector3D & getEY() const;

    /*
    get p-axis, based on O. Fischer (1993)
    */
    const Vector3D & getEZ() const;

    void setEX(Vector3D _e);

    void setEY(Vector3D _e);

    void setEZ(Vector3D _e);

    void setCoordinateSystem(const Vector3D & _ex, const Vector3D & _ey, const Vector3D & _ez);

    void getCoordinateSystem(Vector3D * _ex, Vector3D * _ey, Vector3D * _ez) const;

    void setTmpPathLength(double val);

    void setPositionCell(cell_basic * val);

    void setDirectionID(uint val);

    uint getDirectionID();

    void setPhotonID(ullong _photonID);

    ullong getPhotonID();

  private:
    Vector3D pos;
    Vector3D backup_pos;
    Vector3D ex; // r-axis, based on O. Fischer (1993)
    Vector3D ey; // l-axis
    Vector3D ez; // p-axis
    Matrix2D mD;
    StokesVector * multi_stokes;

    ullong photonID;

    uint dirID;
    uint i_spectral;
    uint nr_of_spectral_bins;

    double tmp_path;
    double trans_frequency;

    double * velocity;
    double * wavelength;
    uint * wID;

    cell_basic * cell_pos;
};

#endif
