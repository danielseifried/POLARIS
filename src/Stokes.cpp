#include "Stokes.hpp"
#include "Vector3D.hpp"


double StokesVector::iPol()
{
    return sqrt(sU * sU + sQ * sQ);
}

double StokesVector::tPol()
{
    return sqrt(sU * sU + sQ * sQ + sV * sV);
}

double StokesVector::linPol()
{
    if(sI != 0)
        return sqrt(sU * sU + sQ * sQ) / sI;

    return 0;
}

double StokesVector::circPol()
{
    if(sI != 0)
        return sV / sI;

    return 0;
}

double StokesVector::getAngle()
{
    return 0.5 * Vector3D::angle(sQ, sU);
}

void StokesVector::setI(double _I)
{
    sI = _I;
}

void StokesVector::setQ(double _Q)
{
    sQ = _Q;
}

void StokesVector::setU(double _U)
{
    sU = _U;
}

void StokesVector::setV(double _V)
{
    sV = _V;
}

void StokesVector::setT(double _T)
{
    sT = _T;
}

void StokesVector::setSp(double _Sp)
{
    sSp = _Sp;
}

void StokesVector::set(double _I, double _Q, double _U, double _V, double _T)
{
    sI = _I;
    sQ = _Q;
    sU = _U;
    sV = _V;
    sT = _T;
    sSp = 0;
}

void StokesVector::set(double _I, double _Q, double _U, double _V)
{
    sI = _I;
    sQ = _Q;
    sU = _U;
    sV = _V;
    sT = 0;
    sSp = 0;
}

void StokesVector::set(StokesVector _S)
{
    sI = _S.I();
    sQ = _S.Q();
    sU = _S.U();
    sV = _S.V();
    sT = _S.T();
    sSp = _S.Sp();
}

void StokesVector::addI(double _I)
{
    sI += _I;
}

void StokesVector::addQ(double _Q)
{
    sQ += _Q;
}

void StokesVector::addU(double _U)
{
    sU += _U;
}

void StokesVector::addV(double _V)
{
    sV += _V;
}

void StokesVector::addT(double _T)
{
    sT += _T;
}

void StokesVector::addSp(double _Sp)
{
    sSp += _Sp;
}

void StokesVector::addS(StokesVector _S)
{
    sI += _S.I();
    sQ += _S.Q();
    sU += _S.U();
    sV += _S.V();
    sT += _S.T();
    sSp += _S.Sp();
}

void StokesVector::multI(double _I)
{
    sI *= _I;
}

void StokesVector::multQ(double _Q)
{
    sQ *= _Q;
}

void StokesVector::multU(double _U)
{
    sU *= _U;
}

void StokesVector::multV(double _V)
{
    sV *= _V;
}

void StokesVector::multT(double _T)
{
    sT *= _T;
}

void StokesVector::multSp(double _Sp)
{
    sSp *= _Sp;
}

void StokesVector::multS(double _S)
{
    sI *= _S;
    sQ *= _S;
    sU *= _S;
    sV *= _S;
}

double StokesVector::I() const
{
    return sI;
}

double StokesVector::Q() const
{
    return sQ;
}

double StokesVector::U() const
{
    return sU;
}

double StokesVector::V() const
{
    return sV;
}

double StokesVector::T() const
{
    return sT;
}

double StokesVector::Sp() const
{
    return sSp;
}

void StokesVector::rot(double phi)
{
    double tQ = sQ, tU = sU;
    double s = sin(2.0 * phi), c = cos(2.0 * phi);
    // The scattering part is based on O. Fischer (1993)
    //
    // The equations
    // Q' = Q*cos(2*a) + U*sin(2*a)
    // U' = -Q*sin(2*a) + U*cos(2*a)
    // transform the stokes parameters if the rotation is a counterclockwise
    // rotation by the angle a > 0 when locking in the direction of propagation.
    //
    // However, when rotating the r- and l-axis by the scattering angle phi > 0,
    // the rotation is defined clockwise when locking in the direction of propagation.
    // Therefore, it is
    // Q' = Q*cos(2*phi) - U*sin(2*phi)
    // U' = Q*sin(2*phi) + U*cos(2*phi)
    sQ = tQ * c - tU * s;
    sU = tQ * s + tU * c;
}

void StokesVector::rot(double sin_phi, double cos_phi)
{
    double tQ = sQ, tU = sU;
    sQ = tQ * cos_phi - tU * sin_phi;
    sU = tQ * sin_phi + tU * cos_phi;
}

bool StokesVector::isConsistent()
{
    if(sI * sI < sQ * sQ + sU * sU + sV * sV)
        return false;

    if(sI < 0)
        return false;

    if(sT < 0)
        return false;

    if(sSp < 0)
        return false;

    return true;
}

void StokesVector::normalize()
{
    sQ /= sI;
    sU /= sI;
    sV /= sI;
    sI = 1;
}

void StokesVector::clear()
{
    sI = 0;
    sQ = 0;
    sU = 0;
    sV = 0;
    sT = 0;
    sSp = 0;
}

void StokesVector::resetIntensity()
{
    sI = 0;
    sQ = 0;
    sU = 0;
    sV = 0;
}

void StokesVector::depolarize()
{
    sQ = 0;
    sU = 0;
    sV = 0;
}

StokesVector & StokesVector::operator=(const StokesVector & ex)
{
    sI = ex.I();
    sQ = ex.Q();
    sU = ex.U();
    sV = ex.V();
    sT = ex.T();
    sSp = ex.Sp();

    return *this;
}

StokesVector & StokesVector::operator=(double v)
{
    sI = v;
    sQ = v;
    sU = v;
    sV = v;
    sT = v;
    sSp = v;
    return *this;
}

StokesVector StokesVector::operator+(const StokesVector & ex) const
{
    return StokesVector(sI + ex.I(), sQ + ex.Q(), sU + ex.U(), sV + ex.V(), sT + ex.T(), sSp + ex.Sp());
}

StokesVector StokesVector::operator-(const StokesVector & ex) const
{
    return StokesVector(sI - ex.I(), sQ - ex.Q(), sU - ex.U(), sV - ex.V(), sT - ex.T(), sSp - ex.Sp());
}

StokesVector & StokesVector::operator+=(const StokesVector & ex)
{
    sI += ex.I();
    sQ += ex.Q();
    sU += ex.U();
    sV += ex.V();
    sT += ex.T();
    sSp += ex.Sp();
    return *this;
}

StokesVector & StokesVector::operator-=(const StokesVector & ex)
{
    sI -= ex.I();
    sQ -= ex.Q();
    sU -= ex.U();
    sV -= ex.V();
    sT -= ex.T();
    sSp -= ex.Sp();
    return *this;
}

StokesVector & StokesVector::operator*=(double val)
{
    sI *= val;
    sQ *= val;
    sU *= val;
    sV *= val;
    return *this;
}

StokesVector & StokesVector::operator*=(const Matrix2D & dM)
{
    double tmp_sI = sI * dM(0, 0) + sQ * dM(0, 1) + sU * dM(0, 2) + sV * dM(0, 3);
    double tmp_sQ = sI * dM(1, 0) + sQ * dM(1, 1) + sU * dM(1, 2) + sV * dM(1, 3);
    double tmp_sU = sI * dM(2, 0) + sQ * dM(2, 1) + sU * dM(2, 2) + sV * dM(2, 3);
    double tmp_sV = sI * dM(3, 0) + sQ * dM(3, 1) + sU * dM(3, 2) + sV * dM(3, 3);
    sI = tmp_sI;
    sQ = tmp_sQ;
    sU = tmp_sU;
    sV = tmp_sV;
    return *this;
}

StokesVector & StokesVector::operator/=(double val)
{
    sI /= val;
    sQ /= val;
    sU /= val;
    sV /= val;
    return *this;
}

ostream & operator<<(ostream & out, const StokesVector & ex)
{
    out << " I: " << ex.I() << " Q: " << ex.Q();
    out << " U: " << ex.U() << " V: " << ex.V() << " T: " << ex.T() << " ";
    return out;
}

StokesVector operator*(const Matrix2D & dM, const StokesVector & v)
{
    return StokesVector(v.I() * dM(0, 0) + v.Q() * dM(0, 1) + v.U() * dM(0, 2) + v.V() * dM(0, 3),
                        v.I() * dM(1, 0) + v.Q() * dM(1, 1) + v.U() * dM(1, 2) + v.V() * dM(1, 3),
                        v.I() * dM(2, 0) + v.Q() * dM(2, 1) + v.U() * dM(2, 2) + v.V() * dM(2, 3),
                        v.I() * dM(3, 0) + v.Q() * dM(3, 1) + v.U() * dM(3, 2) + v.V() * dM(3, 3));
}

StokesVector operator*(const StokesVector & v, const StokesVector & u)
{
    return StokesVector(
        v.I() * u.I(), v.Q() * u.Q(), v.U() * u.U(), v.V() * u.V(), v.T() * u.T(), v.Sp() * u.Sp());
}

StokesVector operator*(const StokesVector & v, double val)
{
    return StokesVector(v.I() * val, v.Q() * val, v.U() * val, v.V() * val, v.T(), v.Sp());
}

StokesVector operator/(const StokesVector & v, double val)
{
    return StokesVector(v.I() / val, v.Q() / val, v.U() / val, v.V() / val, v.T(), v.Sp());
}
