#include <iostream>
#include "integrator.h"
#include "custom.h"
#include "model.h"

using namespace std;

void writeRes(TModel* model, std::string filename);

int main()
{
    TDormandPrinceIntegrator * integrator = new TDormandPrinceIntegrator();
    integrator->setPrecision(1E-16);
    long double * X0 = new long double[6];
    X0[0] = 0;
    X0[1] = 0;
    X0[2] = 0;
    X0[3] = 7000E3L;
    X0[4] = 0.0L;
    X0[5] = 0.0L;
    std::string filename = "";
    Central_Grav_Model * c_grav_model = new Central_Grav_Model(X0);
    Normal_Grav_Model * n_grav_model = new Normal_Grav_Model(X0);
    Point_Grav_Model * p_grav_model = new Point_Grav_Model(X0);
    std::cout << " c_grav_model " << std::endl;
    integrator->Run(c_grav_model); filename = "c_grav.txt"; writeRes(c_grav_model, filename); delete c_grav_model;
    std::cout << " n_grav_model " << std::endl;
    integrator->Run(n_grav_model); filename = "n_grav.txt"; writeRes(n_grav_model, filename); delete n_grav_model;
    std::cout << " p_grav_model " << std::endl;
    integrator->Run(p_grav_model); filename = "p_grav.txt"; writeRes(p_grav_model, filename); delete p_grav_model;
    std::cout << " finish " << std::endl;
    delete[] X0;
    delete integrator;



}

void writeRes(TModel* model, std::string filename){

    std::ofstream file(filename);

        //TMatrix Result = model->getResult();

        for (int i=0; i<model->Result.rowCount(); i++)
        {
            for (int j=0; j<model->Result.colCount(); j++)
                file << model->Result(i, j) << " ";

            file << std::endl;
        }

        file.close();
}
