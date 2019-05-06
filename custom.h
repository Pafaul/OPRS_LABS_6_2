#ifndef CUSTOM_H
#define CUSTOM_H
#include "math.h"
#include <fstream>
#include <sstream>
#include "model.h"
#include "aes_param.h"

long double string_to_double(std::string s);

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
    long double grav_points[5*60] {1.00, -1917861.343, -1597.53455, 3389.08854, -1206.07844, 2.00, -7649.811, 5243.88105, 2173.09105, 31.68, 3.00, -23204717.367, -694.74764, -2543.95209, 3010.82934, 4.00, -8601525.203, 597.26, 4124.29275, 3032.72188, 5.00, -39262968.108, 1472.01984, 2496.13550, 2305.68074, 6.00, -37656928.613, -1395.72710, 522.58, 1219.80230, 7.00, 55417.233, -4729.43278, 241.06, -1386.75430, 8.00, -1122388.574, 1613.74542, -1389.46971, -4742.37198, 9.00, 11553321.883, -737.33302, 4223.69658, -89.84823, 10.00, -3212166.768, -3498.93844, 2123.53095, -2780.28445, 11.00, -127570.128, 2443.27977, -2075.43783, -3115.11405, 12.00, 537718.408, 1399.38281, -2714.24422, -1999.26660, 13.00, 159431.526, -1807.95598, 4038.16007, -1984.43463, 14.00, -9621707.493, 500.99, 2417.51767, 2081.23547, 15.00, 14954047.019, 564.56, 4169.62120, 3107.28151, 16.00, -1454978.047, -3606.22095, 2202.63232, 2555.27655, 17.00, 12173.463, 3494.94521, -2225.41218, 3504.56634, 18.00, 6235318.028, 2819.10268, 1215.19399, 3031.84264, 19.00, -4022993.036, 830.79, 1187.36734, 2665.51499, 20.00, -13233537.348, -2407.50792, -987.75326, 1168.97751, 21.00, -837955.460, 2745.40696, 3312.38455, 561.61, 22.00, 66944.487, 474.28, 1205.35801, -4546.03756, 23.00, 9533469.996, -1675.28817, -797.14875, -3048.67319, 24.00, 13126.101, 1951.39223, -5027.13641, -1609.06073, 25.00, 133490487.677, -1131.17750, -625.93077, 614.82, 26.00, 11202.018, 4552.11518, 1485.42650, -1997.43144, 27.00, -8306318.912, 32.86, 1869.79667, -2877.64543, 28.00, 1563860.673, -3588.59466, 2192.00472, 2532.32304, 29.00, -4110249.075, -1688.88208, -2509.59158, 997.40, 30.00, -708432.828, 142.71, -2365.17339, 3559.90481, 31.00, 115758.347, -2583.62219, 3876.78111, -1530.49081, 32.00, 47117.768, -4410.11412, -1491.17613, 1951.16572, 33.00, 3060709.563, -3512.99151, 2136.54197, -2789.00902, 34.00, 125694.630, 2887.19498, -1520.35002, -2929.00588, 35.00, -19909835.138, -1221.99104, -1430.21334, 2886.63494, 36.00, -25600026.158, -795.89525, 4094.06397, -148.43103, 37.00, 8998492.144, -859.93100, -2675.58249, 3049.08749, 38.00, 6460340.213, 134.47, -3777.76920, -737.29613, 39.00, 22635646.931, 2030.13210, 2786.99103, 1094.12609, 40.00, 14963823.854, -894.96828, 3946.71389, -241.43049, 41.00, 104200607.824, 1424.84214, 2420.76271, 1954.22532, 42.00, -29025.051, 3640.95418, 2890.38428, -1299.38754, 43.00, -25923827.511, -1666.65137, -1236.38163, -177.16683, 44.00, 4710.121, -2926.62377, 1305.39588, 4683.40996, 45.00, -72296468.140, 1704.30580, 2570.01705, 1437.94150, 46.00, 11488188.566, -1329.71266, -1480.67585, 3122.43463, 47.00, 30871.376, 2993.13075, -3221.87772, 2310.62018, 48.00, -6821341.595, 539.12, 4204.18554, 3168.98432, 49.00, 4053215.384, 520.11, -3574.66707, -271.75382, 50.00, 2364380.549, -1696.20051, 1611.26783, 2166.50206, 51.00, -5107799.845, -2806.42940, 90.03, -437.20618, 52.00, 1152914.935, 1613.29970, -1398.25374, -4728.63714, 53.00, 14079228.613, 1331.61419, -1428.22398, 234.33, 54.00, 24148.390, 2240.63930, 2260.34015, -3797.81205, 55.00, 14981240.134, -523.57471, -2437.64320, 3040.87357, 56.00, -7746811.230, 2751.33446, 1231.49383, 2975.90618, 57.00, -52402113.809, 490.28, -1222.02040, 299.60, 58.00, -9702053.744, 214.47, -3717.91494, -621.72796, 59.00, 6998066.913, 43.27, 1953.11104, -2943.24995, 60.00, -11022424.436, -1649.80829, -755.81754, -2992.38699};
public:
    Point_Grav_Model();
    Point_Grav_Model(long double * init_param);
    void getRight( const TVector& X, long double t, TVector& Y );
};

#endif // CUSTOM_H
