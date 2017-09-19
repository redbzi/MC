#ifndef MC_BASIC_H
#define MC_BASIC_H

#include <cmath>


namespace GBM {
    namespace Analytic {
        inline
        double
        norm_pdf(const double &x) {
            return (1. / (pow(2 * M_PI, 0.5))) * exp(-0.5 * x * x);
        }

        inline
        double
        norm_cdf(const double &x) {
            double k = 1. / (1. + 0.2316419 * x);
            double k_sum =
                    k * (0.319381530 + k * (-0.356563782 + k * (1.781477937 + k * (-1.821255978 + 1.330274429 * k))));

            if (x >= 0.0) {
                return (1. - (1. / (pow(2 * M_PI, 0.5))) * exp(-0.5 * x * x) * k_sum);
            } else {
                return 1. - norm_cdf(-x);
            }
        }

        inline
        double
        d_j(const int &j, const double &S, const double &K, const double &r, const double &v, const double &T) {
            return (log(S / K) + (r + (pow(-1, j - 1)) * 0.5 * v * v) * T) / (v * (pow(T, 0.5)));
        }

        inline
        double
        call_price(const double &S, const double &K, const double &r, const double &v, const double &T) {
            return S * norm_cdf(d_j(1, S, K, r, v, T)) - K * exp(-r * T) * norm_cdf(d_j(2, S, K, r, v, T));
        }

        inline
        double
        call_delta(const double &S, const double &K, const double &r, const double &v, const double &T) {
            return norm_cdf(d_j(1, S, K, r, v, T));
        }

        inline
        double
        call_gamma(const double &S, const double &K, const double &r, const double &v, const double &T) {
            return norm_pdf(d_j(1, S, K, r, v, T)) / (S * v * sqrt(T));
        }

        inline
        double
        call_vega(const double &S, const double &K, const double &r, const double &v, const double &T) {
            return S * norm_pdf(d_j(1, S, K, r, v, T)) * sqrt(T);
        }
    }

    namespace MC {
        namespace Basic {
            std::vector<double>
            call_price(const double &S_0, const double &K, const double &r, const double &sigma, const double &T,
                       const int &M);

            std::vector<double>
            call_delta(const double &S_0, const double &K, const double &r, const double &sigma, const double &T,
                       const int &M, const double &d_S);

            std::vector<double>
            call_gamma(const double &S_0, const double &K, const double &r, const double &sigma, const double &T,
                       const int &M, const double &d_S);

            std::vector<double>
            call_vega(const double &S_0, const double &K, const double &r, const double &sigma, const double &T,
                      const int &M, const double &d_sigma);
        }
        namespace Antithetic {
            std::vector<double>
            call_price(const double &S_0, const double &K, const double &r, const double &sigma, const double &T,
                       const int &M);

            std::vector<double>
            call_delta(const double &S_0, const double &K, const double &r, const double &sigma, const double &T,
                       const int &M, const double &d_S);

            std::vector<double>
            call_gamma(const double &S_0, const double &K, const double &r, const double &sigma, const double &T,
                       const int &M, const double &d_S);

            std::vector<double>
            call_vega(const double &S_0, const double &K, const double &r, const double &sigma, const double &T,
                      const int &M, const double &d_sigma);
        }
        namespace ControlVariate {
            std::vector<double>
            call_price(const double &S_0, const double &K, const double &r, const double &sigma, const double &T,
                       const int &M);
        }
        namespace ImportantSampling{
            std::vector<double>
            call_price(const double &S_0, const double &K, const double &r, const double &sigma, const double &T,
                       const int &M);

            std::vector<double>
            call_delta(const double &S_0, const double &K, const double &r, const double &sigma, const double &T,
                       const int &M, const double &d_S);

            std::vector<double>
            call_gamma(const double &S_0, const double &K, const double &r, const double &sigma, const double &T,
                       const int &M, const double &d_S);

            std::vector<double>
            call_vega(const double &S_0, const double &K, const double &r, const double &sigma, const double &T,
                      const int &M, const double &d_sigma);
        }
    }
}

#endif //MC_BASIC_H
