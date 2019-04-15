#ifndef AES_PARAM_H
#define AES_PARAM_H
#include "LA.h"

class AES_Param
{
protected:
    TVector r_or, v_or; //в оскулирующих элементах
    TMatrix A; //переход в геоцентрицеские координаты
public:
    AES_Param();
    AES_Param(long double * param);
    TVector R, V, X;
    void toGeocentic();
    void recalc_orbit_param();
    void recalc_p();
    long double omega, //долгота восходящего узла
                i, //наклонение орбиты
                w, //широта перицентра
                a, //большая полуось орбиты
                e, //эксцентриситет
                fi, //истинная аномалия на начало эксперимента
                p, //фокальный параметр
                mu = 398600.436E9L; //гравитационная постоянная
};

#endif // AES_PARAM_H
