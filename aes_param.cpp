#include "math.h"
#include "aes_param.h"

AES_Param::AES_Param() {}

AES_Param::AES_Param(long double * param)
{
    //всё приходит в радианах обычно
    omega = param[0];
    i = param[1];
    w = param[2];
    a = param[3];
    e = param[4];
    fi = param[5];

    r_or.resize(3); v_or.resize(3);
    R.resize(3); V.resize(3);
    X.resize(6);
    A.resize(3, 3);
    recalc_orbit_param();
}

void AES_Param::recalc_orbit_param()
{
    recalc_p();
    long double u = w + fi;
    r_or[0] = p/(1 + e*cos(fi));
    r_or[1] = 0;
    r_or[2] = 0;
    v_or[0] = sqrt(mu/p)*e*sin(fi);
    v_or[1] = sqrt(mu/p)*(1 + e*cos(fi));
    v_or[2] = 0;
    A(0, 0) = cos(u)*cos(omega) - sin(u)*sin(omega)*cos(i);
    A(0, 1) = -sin(u)*cos(omega) - cos(u)*sin(omega)*cos(i);
    A(0, 2) = sin(i)*sin(omega);
    A(1, 0) = cos(u)*sin(omega) + sin(u)*cos(omega)*cos(i);
    A(1, 1) = -sin(u)*sin(omega) + cos(u)*cos(omega)*cos(i);
    A(1, 2) = -sin(i)*cos(omega);
    A(2, 0) = sin(u)*sin(i);
    A(2, 1) = cos(u)*sin(i);
    A(2, 2) = cos(i);
    R = A*r_or;
    V = A*v_or;
    for (int j = 0; j < 3; j++) { X[j] = R[j]; X[j+3] = V[j]; }
}

void AES_Param::recalc_p() { p = a*(1 - e*e); }
