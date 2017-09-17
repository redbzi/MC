#ifndef SM_BASIC_H
#define SM_BASIC_H
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <vector>
#include <ctime>

namespace GBM {
    namespace MC {
        namespace Direct {
            std::vector<double> call_price(double &S_0, double &K, double &r, double &sigma, double &T, int &M);

            std::vector<double>
            call_delta(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double &d_S);

            std::vector<double>
            call_gamma(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double &d_S);

            std::vector<double>
            call_vega(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double &d_vega);
        }
        namespace Antithetic {
            std::vector<double> call_price(double &S_0, double &K, double &r, double &sigma, double &T, int &M);

            std::vector<double>
            call_delta(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double &d_S);

            std::vector<double>
            call_gamma(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double &d_S);

            std::vector<double>
            call_vega(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double &d_vega);
        }
        namespace ControlVariate {
            std::vector<double> call_price(double &S_0, double &K, double &r, double &sigma, double &T, int &M);
        }
        namespace ImportantSampling {
            std::vector<double> call_price(double &S_0, double &K, double &r, double &sigma, double &T, int &M);

            std::vector<double>
            call_delta(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double d_S);

            std::vector<double>
            call_gamma(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double d_S);

            std::vector<double>
            call_vega(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double d_vega);
        }
    }
}
#endif //SM_BASIC_H
