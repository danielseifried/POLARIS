/************************************************************************************
*                      POLARIS: POLArized RadIation Simulator                       *
*                         Copyright (C) 2018 Stefan Reissl                          *
************************************************************************************/

#ifndef STOKESVECTOR_H
#define STOKESVECTOR_H

#include "Matrix2D.hpp"
#include "Typedefs.hpp"

class StokesVector
{
public:
    StokesVector()
    {
        sI = 0;
        sQ = 0;
        sU = 0;
        sV = 0;
        sT = 0;
        sSp = 0;
    }

    StokesVector(double val)
    {
        sI = val;
        sQ = val;
        sU = val;
        sV = val;
        sT = val;
        sSp = val;
    }

    StokesVector(double I, double Q, double U, double V)
    {
        sI = I;
        sQ = Q;
        sU = U;
        sV = V;
        sT = 0;
        sSp = 0;
    }

    StokesVector(double I, double Q, double U, double V, double T)
    {
        sI = I;
        sQ = Q;
        sU = U;
        sV = V;
        sT = T;
        sSp = 0;
    }

    StokesVector(double I, double Q, double U, double V, double T, double Sp)
    {
        sI = I;
        sQ = Q;
        sU = U;
        sV = V;
        sT = T;
        sSp = Sp;
    }

    StokesVector(const StokesVector & st)
    {
        sI = st.I();
        sQ = st.Q();
        sU = st.U();
        sV = st.V();
        sT = st.T();
        sSp = st.Sp();
    }

    ~StokesVector(void)
    {}

    // linearly polarized intensity
    double iPol();

    // totaly polarized intensity
    double tPol();

    // degree of linear polarization
    double linPol();

    // degree of circular polarization
    double circPol();

    // polarization angle
    double getAngle();

    void setI(double _I);

    void setQ(double _Q);

    void setU(double _U);

    void setV(double _V);

    void setT(double _T);

    void setSp(double _Sp);

    void set(double _I, double _Q, double _U, double _V, double _T);

    void set(double _I, double _Q, double _U, double _V);

    void set(StokesVector _S);

    void addI(double _I);

    void addQ(double _Q);

    void addU(double _U);

    void addV(double _V);

    void addT(double _T);

    void addSp(double _Sp);

    void addS(StokesVector _S);

    void multI(double _I);

    void multQ(double _Q);

    void multU(double _U);

    void multV(double _V);

    void multT(double _T);

    void multSp(double _Sp);

    void multS(double _S);

    double I() const;

    double Q() const;

    double U() const;

    double V() const;

    double T() const;

    double Sp() const;

    void rot(double phi);

    void rot(double sin_phi, double cos_phi);

    bool isConsistent();

    void normalize();

    void clear();

    void resetIntensity();

    void depolarize();

    StokesVector & operator=(const StokesVector & ex);

    StokesVector & operator=(double v);

    StokesVector operator+(const StokesVector & ex) const;

    StokesVector operator-(const StokesVector & ex) const;

    StokesVector & operator+=(const StokesVector & ex);

    StokesVector & operator-=(const StokesVector & ex);

    StokesVector & operator*=(double val);

    StokesVector & operator*=(const Matrix2D & dM);

    StokesVector & operator/=(double val);

    friend ostream & operator<<(ostream & out, const StokesVector & ex);

    friend StokesVector operator*(const Matrix2D & dM, const StokesVector & v);

    friend StokesVector operator*(const StokesVector & v, const StokesVector & u);

    friend StokesVector operator*(const StokesVector & v, double val);

    friend StokesVector operator/(const StokesVector & v, double val);

private:
    double sI, sQ, sU, sV, sT, sSp;
};

#endif /* STOKESVECTOR_H */
