#ifndef SM_HESTON_H
#define SM_HESTON_H
#include <cmath>
#include <algorithm>
#include <vector>
#include <ctime>
#include <iostream>
#include "analytic.h"

namespace Heston {
    namespace MC {
        namespace Direct {
            std::vector<double>
            call_price(double &S_0, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho,
                       double &T, int &M);

            std::vector<double>
            call_delta(double &S_0, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho,
                       double &T, int &M, double &delta_S);

            std::vector<double>
            call_gamma(double &S_0, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho,
                       double &T, int &M, double &delta_S);

            std::vector<double>
            call_vega_V_t(double &S_0, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho,
                      double &T, int &M, double &delta_sigma);

            std::vector<double>
            call_vega_V_0(double &S_0, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho,
                      double &T, int &M, double &delta_sigma);

            std::vector<double>
            call_vega_theta(double &S_0, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho,
                      double &T, int &M, double &delta_theta);
        }
        namespace Antithetic {
            std::vector<double>
            call_price(double &S_0, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho,
                       double &T, int &M);

            std::vector<double>
            call_delta(double &S_0, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho,
                       double &T, int &M, double &delta_S);

            std::vector<double>
            call_gamma(double &S_0, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho,
                       double &T, int &M, double &delta_S);

            std::vector<double>
            call_vega_V_t(double &S_0, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho,
                      double &T, int &M, double &delta_sigma);
        }
        namespace ControlVariate {
            std::vector<double>
            call_price(double &S_0, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho,
                       double &T, int &M);
        }
    }
}
#endif //SM_HESTON_H
