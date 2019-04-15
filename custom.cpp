#include "custom.h"

GravModel::GravModel()
{

}

GravModel::GravModel(long double * init_param)
{
    AES = new AES_Param (init_param);
    X0=AES->X;
}

Central_Grav_Model::Central_Grav_Model(long double * init_param) : GravModel (init_param)
{
};

void Central_Grav_Model::getRight(const TVector &X, long double t, TVector &Y)
{
    //test
    long double ro = sqrt(pow(X[0], 2) + pow(X[1], 2) + pow(X[2], 2));
    Y[0] = X[3];
    Y[1] = X[4];
    Y[2] = X[5];
    Y[3] = -mu*X[0]/pow(ro, 3);
    Y[4] = -mu*X[1]/pow(ro, 3);
    Y[5] = -mu*X[2]/pow(ro, 3);
}

Normal_Grav_Model::Normal_Grav_Model(long double * init_param) : GravModel (init_param)
{
    J.resize(9); C.resize(9); gSphere.resize(3); temp.resize(3);
    A.resize(3, 3);
    for (int i = 0; i < J.size(); i++ ) { J[i] = 0; C[i] = 0; }
    J[2] = 1082.62575E-6L;
    J[4] = -2.37089E-6L;
    J[6] = 6.08E-9L;
    J[8] = -1.4E-11L;
    for (int i = 0; i < J.size(); i++)
        C[i] = -J[i]/sqrt(2*static_cast<long double>(i) + 1.0L);
}

long double Normal_Grav_Model::sigma(int m)
{
    if (m == 0) return 1/2.0L;
    else return 1.0L;
}

long double Normal_Grav_Model::Pnm(int n, int m, long double fi)
{
    if ( (n == m) && (n != 0) && (m != 0))
        return Pnm(n-1, m-1, fi)*cos(fi)*sqrt((2*n+1)/(2*n*sigma(m-1)));
    if ( n > m )
        return Pnm(n-1, m, fi)*sin(fi)*sqrt((4*n*n - 1)/(n*n - m*m)) -
                Pnm(n - 2, m, fi)*sin(fi)*sqrt(((pow(n-1, 2)-m*m)*(2*n+1))/((n*n-m*m)*(2*n-3)));
    if ( n < m )
        return 0.0L;
    if ( n == m && n == 0 && m == 0 )
        return 1.0L;
}

long double Normal_Grav_Model::Pder(int n, long double fi)
{
    return sqrt(1/2.0L*n*(n+1))*Pnm(n, 1, fi);
}

void Normal_Grav_Model::calcSphere(const TVector &X, long double t)
{
    long double temp_L = 0;
    ro = sqrt(pow(X[0], 2) + pow(X[1], 2) + pow(X[2], 2));
    fi = atan2(X[2], sqrt(pow(X[0], 2) + pow(X[1], 2)));
    lambda = atan2(X[1], X[0]);

    for (int i = 2; i < 9; i+=2)
        temp_L += ( i + 1 )* pow((AES->a/ro), i)*C[i]*Pder(i, fi);
    gSphere[0] = -mu/pow(ro, 2) * (1 + temp_L);
    for (int i = 2; i < 9; i+=2)
        temp_L += pow(AES->a/ro, i)*C[i]*Pder(i, fi);
    gSphere[1] = mu/pow(ro, 2)*temp_L;
    gSphere[2] = 0.0L;
    calcTransition(X, ro, fi, lambda);
}

void Normal_Grav_Model::calcTransition(const TVector &X, long double ro, long double fi, long double lambda)
{
    long double R0 = sqrt(X[0]*X[0] + X[1]*X[1]);
    A(0, 0) = X[0]/ro;
    A(0, 1) = -X[0]*X[2]/(ro*R0);
    A(0, 2) = -X[1]/R0;
    A(1, 0) = X[1]/ro;
    A(1, 1) = -X[1]*X[2]/(ro*R0);
    A(1, 2) = X[0]/R0;
    A(2, 0) = X[2]/ro;
    A(2, 1) = R0/ro;
    A(2, 2) = 0.0L;
}

void Normal_Grav_Model::getRight(const TVector &X, long double t, TVector &Y)
{
    calcSphere(X, t);
    temp = A*gSphere;
    for (int i = 0; i < 3; i++)
    {
        Y[i] = X[i+3];
        Y[i+3] = temp[i];
    }
}

Point_Grav_Model::Point_Grav_Model()
{
}

Point_Grav_Model::Point_Grav_Model(long double * init_param) : GravModel (init_param)
{
    std::ifstream file("grav_points.txt");
    int i = 0;
    int j = 0;
    while(file >> j >> grav_points[i*4] >> grav_points[i*4+1] >> grav_points[i*4+2] >> grav_points[i*4+3]) i ++;
    for (i = 0; i < 60; i++) for ( j = 0; j < 3; j++) grav_points[i*4+j+1] *= 1000;
}

void Point_Grav_Model::getRight(const TVector &X, long double t, TVector &Y)
{
    double ro = 0;
    Y[0] = X[3];
    Y[1] = X[4];
    Y[2] = X[5];
    for (int i = 0; i < 60; i++) {
        ro = 0;
        for (int j = 0; j < 3; j++) ro += pow(X[j]-grav_points[i*4+j+1], 2);
        ro = sqrt(ro);
        Y[3] += -mu*(X[0]-grav_points[i*4+1])/pow(ro, 3);
        Y[4] += -mu*(X[1]-grav_points[i*4+2])/pow(ro, 3);
        Y[5] += -mu*(X[2]-grav_points[i*4+3])/pow(ro, 3);
    }
}
