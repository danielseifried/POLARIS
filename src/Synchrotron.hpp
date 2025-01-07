#pragma once

#include "MathFunctions.hpp"
#include "MathSpline.hpp"
#include "SyncParameters.hpp"
#include "Matrix2D.hpp"
#include "Typedefs.hpp"


#ifndef SYNCHROTRON_H
#define SYNCHROTRON_H

// class for the physics of sync. RT
class CSynchrotron
{
  public:
    CSynchrotron()
    {
        initGamma();
        initBesselK();
    };

    ~CSynchrotron()
    {
        if(p != 0)
        {
            delete[] p;
            p = 0;
        }

        if(q != 0)
        {
            delete[] q;
            q = 0;
        }

        if(c != 0)
        {
            delete[] c;
            c = 0;
        }
    };

    /*
    calculation of sync. coefficients on the basis of approximation functions (taken
    from Pandya 2016 and Dexter 2016)
    */
    syn_param get_Thermal_Parameter(double n_e, double T_e, double l, double B, double theta);
    syn_param get_Power_Law_Parameter(double n_e,
                                      double l,
                                      double B,
                                      double theta,
                                      double g_min,
                                      double g_max,
                                      double p);

  private:
    // Gamma numerator coefficients for approximation over the interval (1,2)
    double * p;

    // Gamma denominator coefficients for approximation over the interval (1,2)
    double * q;

    // Asymptotic series LogGamma, Abramowitz and Stegun 6.1.41
    double * c;

    double Euler_gamma;
    double halfLogTwoPi;

    spline BesselK_0;
    spline BesselK_1;
    spline BesselK_2;

    double BesselK(uint n, double x);

    // initiate Gamma coefficients
    void initGamma();

    // initiate Bessel points for spline interpolation
    void initBesselK();

    double LogGamma(double x);
    double Gamma(double x);

    // additional correction function for alpha_V (see Reissl et al. 2018)
    double corr(double theta);

    double Gamma_I_p(double g_min, double g_max, double p);

    double Gamma_A_p(double g_min, double g_max, double p);

    double getI_Q_p(double p);

    double getI_V_p(double p, double tan_theta);

    double getA_Q_p(double p);

    double getA_V_p(double p, double sin_theta, double cos_theta);

    double getK_Q_th(double x);

    double getK_V_th(double x);

    double getI_I_th(double x);

    double getI_Q_th(double x, double Theta);

    double getI_V_th(double x, double theta, double Theta);

    double planck_hz(double nu, double Theta);
};

#endif
