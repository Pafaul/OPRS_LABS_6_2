#ifndef CUSTOM_H
#define CUSTOM_H
#include "math.h"
#include <fstream>
#include "model.h"
#include "aes_param.h"

struct grav_point
{
    double x, y, z;
    double mass;
    grav_point(double x, double y, double z, double mass) { this->x = x; this->y = y; this->z = z; this->mass = mass; }
};

class GravModel : public TModel
{
protected:
    AES_Param * AES;
    long double mu = 398600.436E9L,
                omega = 7.292115E-5L;
public:
    GravModel();
    GravModel(long double * init_param);
    bool run ( const TVector& X, long double t ) { return true; }
    void do_thing( const TVector& X, long double t ) {}
    void finalize() {}
};

class Central_Grav_Model : public GravModel
{

public:
    Central_Grav_Model();
    Central_Grav_Model( long double * init_param);
    void getRight( const TVector& X, long double t, TVector& Y );
};

class Normal_Grav_Model : public GravModel
{
protected:
    long double ro, fi, lambda;
    TVector J, C, gSphere, temp;
    TMatrix A;
public:
    Normal_Grav_Model();
    Normal_Grav_Model(long double * init_param);
    long double Pnm(int n, int m, long double fi);
    long double Pder(int n, long double fi);
    long double sigma(int m);
    void calcSphere(const TVector& X, long double t);
    void calcTransition(const TVector& X, long double ro, long double fi, long double lambda);
    void getRight( const TVector& X, long double t, TVector& Y);

};

class Point_Grav_Model : public GravModel
{
protected:
    double grav_points[3*60];
public:
    Point_Grav_Model();
    Point_Grav_Model(long double * init_param);
    void getRight( const TVector& X, long double t, TVector& Y );
};

#endif // CUSTOM_H
