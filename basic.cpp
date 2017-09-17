#include "basic.h"

namespace GBM {
    namespace MC {
        namespace Direct {
            std::vector<double> call_price(double &S_0, double &K, double &r, double &sigma, double &T, int &M) {
                std::clock_t t = clock();
                double sum_price = 0.;
                double sum_2_price = 0.;
                double D = exp(-r * T);
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    double S_T = S_0 * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    sum_price += D * std::max(S_T - K, 0.);
                    sum_2_price += std::pow(D * std::max(S_T - K, 0.), 2);
                }
                double m = sum_price / (double) M;
                double v = ((1 / (double) M) * sum_2_price - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_delta(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double &d_S) {
                std::clock_t t = clock();
                double sum_delta = 0.;
                double sum_2_delta = 0.;
                double D = exp(-r * T);
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    double S_T_m = (S_0 - d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_h = (S_0 + d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double price_m = D * std::max(S_T_m - K, 0.);
                    double price_h = D * std::max(S_T_h - K, 0.);
                    sum_delta += (price_h - price_m) / (2 * d_S);
                    sum_2_delta += std::pow((price_h - price_m) / (2 * d_S), 2);
                }
                double m = (sum_delta) / (double) M;
                double v = ((1 / (double) M) * sum_2_delta - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_gamma(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double &d_S) {
                std::clock_t t = clock();

                double D = exp(-r * T);
                double price_l = 0.;
                double price_m = 0.;
                double price_h = 0.;
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    double S_T_l = (S_0 - d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_m =  S_0        * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_h = (S_0 + d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    price_l += D * std::max(S_T_l - K, 0.);
                    price_m += D * std::max(S_T_m - K, 0.);
                    price_h += D * std::max(S_T_h - K, 0.);
                }
                double esp_pl = price_l / M;
                double esp_pm = price_m / M;
                double esp_ph = price_h / M;
                double sum_gamma = (esp_ph - 2 * esp_pm + esp_pl) / (d_S * d_S);
                double sum_2_gamma = std::pow((price_h - 2 * price_m + price_l) / (d_S * d_S), 2);
                double m = sum_gamma;
                double v = ((1 / (double) M) * sum_2_gamma - m * m);
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_vega(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double &d_sigma) {
                std::clock_t t = clock();
                double sum_vega = 0.;
                double sum_2_vega = 0.;
                double D = exp(-r * T);
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    double S_T_m = S_0 * exp(T * (r - 0.5 * (sigma - d_sigma) * (sigma - d_sigma)) + (sigma - d_sigma) * sqrt(T) * W);
                    double S_T_h = S_0 * exp(T * (r - 0.5 * (sigma + d_sigma) * (sigma + d_sigma)) + (sigma + d_sigma) * sqrt(T) * W);
                    double price_m = D * std::max(S_T_m - K, 0.);
                    double price_h = D * std::max(S_T_h - K, 0.);
                    sum_vega += (price_h - price_m) / (2*d_sigma);
                    sum_2_vega += std::pow((price_h - price_m) / d_sigma, 2);
                }
                double m = (sum_vega) / (double) M;
                double v = ((1 / (double) M) * sum_2_vega - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }
        }
        namespace Antithetic {
            std::vector<double> call_price(double &S_0, double &K, double &r, double &sigma, double &T, int &M) {
                std::clock_t t = clock();
                double sum_price = 0.;
                double sum_2_price = 0.;
                double D = exp(-r * T);
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    double S_T_1 = S_0 * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_2 = S_0 * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * (-W));
                    sum_price += D * (std::max(S_T_1 - K, 0.) + std::max(S_T_2 - K, 0.)) / 2;
                    sum_2_price += std::pow(D * (std::max(S_T_1 - K, 0.) + std::max(S_T_2 - K, 0.)) / 2, 2);
                }
                double m = sum_price / (double) M;
                double v = ((1 / (double) M) * sum_2_price - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_delta(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double &d_S) {
                std::clock_t t = clock();
                double sum_delta = 0.;
                double sum_2_delta = 0.;
                double D = exp(-r * T);
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    double S_T_m_1 = (S_0 - d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_h_1 = (S_0 + d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_m_2 = (S_0 - d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * (-W));
                    double S_T_h_2 = (S_0 + d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * (-W));
                    double price_m_1 = D * std::max(S_T_m_1 - K, 0.);
                    double price_h_1 = D * std::max(S_T_h_1 - K, 0.);
                    double price_m_2 = D * std::max(S_T_m_2 - K, 0.);
                    double price_h_2 = D * std::max(S_T_h_2 - K, 0.);
                    sum_delta += ((price_h_1 - price_m_1) / (2 * d_S) + (price_h_2 - price_m_2) / (2 * d_S)) / 2;
                    sum_2_delta += std::pow(((price_h_1 - price_m_1) / (2 * d_S) + (price_h_2 - price_m_2) / (2 * d_S)) / 2, 2);
                }
                double m = (sum_delta) / (double) M;
                double v = ((1 / (double) M) * sum_2_delta - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_gamma(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double &d_S) {
                std::clock_t t = clock();
                double sum_gamma = 0.;
                double sum_2_gamma = 0.;
                double D = exp(-r * T);
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    double S_T_l_1 = (S_0 - d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_m_1 = S_0         * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_h_1 = (S_0 + d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);

                    double S_T_l_2 = (S_0 - d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * (-W));
                    double S_T_m_2 = S_0         * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * (-W));
                    double S_T_h_2 = (S_0 + d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * (-W));

                    double price_l_1 = D * std::max(S_T_l_1 - K, 0.);
                    double price_m_1 = D * std::max(S_T_m_1 - K, 0.);
                    double price_h_1 = D * std::max(S_T_h_1 - K, 0.);

                    double price_l_2 = D * std::max(S_T_l_2 - K, 0.);
                    double price_m_2 = D * std::max(S_T_m_2 - K, 0.);
                    double price_h_2 = D * std::max(S_T_h_2 - K, 0.);
                    sum_gamma += ((price_h_1 - 2 * price_m_1 + price_l_1) / (d_S * d_S) + (price_h_2 - 2 * price_m_2 + price_l_2) / (d_S * d_S)) / 2;
                    sum_2_gamma += std::pow(((price_h_1 - 2 * price_m_1 + price_l_1) / (d_S * d_S) + (price_h_2 - 2 * price_m_2 + price_l_2) / (d_S * d_S)) / 2, 2);
                }
                double m = sum_gamma / (double) M;
                double v = ((1 / (double) M) * sum_2_gamma - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_vega(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double &d_sigma) {
                std::clock_t t = clock();
                double sum_vega = 0.;
                double sum_2_vega = 0.;
                double D = exp(-r * T);
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    double S_T_m_1 = S_0 * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_h_1 = S_0 * exp(T * (r - 0.5 * (sigma + d_sigma) * (sigma + d_sigma)) + (sigma + d_sigma) * sqrt(T) * W);
                    double S_T_m_2 = S_0 * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * (-W));
                    double S_T_h_2 = S_0 * exp(T * (r - 0.5 * (sigma + d_sigma) * (sigma + d_sigma)) + (sigma + d_sigma) * sqrt(T) * (-W));
                    double price_m_1 = D * std::max(S_T_m_1 - K, 0.);
                    double price_h_1 = D * std::max(S_T_h_1 - K, 0.);
                    double price_m_2 = D * std::max(S_T_m_2 - K, 0.);
                    double price_h_2 = D * std::max(S_T_h_2 - K, 0.);
                    sum_vega += ((price_h_1 - price_m_1) / d_sigma + (price_h_2 - price_m_2) / d_sigma) / 2;
                    sum_2_vega += std::pow(((price_h_1 - price_m_1) / d_sigma + (price_h_2 - price_m_2) / d_sigma) / 2, 2);
                }
                double m = (sum_vega) / (double) M;
                double v = ((1 / (double) M) * sum_2_vega - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }
        }
        namespace ControlVariate {
            std::vector<double> call_price(double &S_0, double &K, double &r, double &sigma, double &T, int &M) {
                std::clock_t t = clock();
                //pilot simulation
                double rho = 0.;//correlation between Y and Z, just for information
                double c = 0.;
                int p = 1000;
                double D = exp(-r * T);
                double E_Y = 0.;
                double E_Y_2 = 0.;
                std::vector<double> V_Y;
                std::vector<double> V_Z;
                double Var_Z = 0.;
                //generate Y and Z samples
                for (int i(0); i<p; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    double X = S_0 * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    V_Y.push_back(D * std::max(X - K, 0.));
                    V_Z.push_back(D * X);
                    E_Y += V_Y[i] / p;
                    E_Y_2 += std::pow(V_Y[i], 2) / p;
                    Var_Z += std::pow((V_Z[i] - S_0), 2) / p;
                }
                for (int i(0); i<p; ++i) {
                    c += - (V_Y[i] - E_Y) * (V_Z[i] - S_0) / (p*Var_Z);
                    rho += (V_Y[i] - E_Y) * (V_Z[i] - S_0) / (p*sqrt(Var_Z * (E_Y_2 - std::pow(E_Y, 2))));
                }

                double sum_price = 0.;
                double sum_2_price = 0.;
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    double X = S_0 * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double Y = D * std::max(X - K, 0.);
                    double Z = D * X;
                    double theta_c = Y + c*(Z - S_0);
                    sum_price += theta_c;
                    sum_2_price += std::pow(theta_c, 2);
                }
                double m = sum_price / (double) M;
                double v = ((1 / (double) M) * sum_2_price - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration, rho};
                return res;
            }
        }
        namespace ImportantSampling{
            std::vector<double> call_price(double &S_0, double &K, double &r, double &sigma, double &T, int &M) {
                double sig = 0.7;
                double mu = 1.3;
                std::clock_t t = clock();
                double sum_price = 0.;
                double sum_2_price = 0.;
                double D = exp(-r * T);
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    W = W * sig + mu; //Change of measure
                    double S_T = S_0 * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double OptionPrice = std::max(S_T - K, 0.)*(1/sqrt(2*M_PI))*exp(-pow(W,2)/2)/((1/(sig*sqrt(2*M_PI)))*exp(-0.5*(pow((W-mu)/sig,2))));
                    sum_price += D * OptionPrice;
                    sum_2_price += std::pow(D * OptionPrice, 2);
                }

                double m = sum_price / (double) M;
                double v = ((1 / (double) M) * sum_2_price - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_delta(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double d_S) {
                double sig = 0.7;
                double mu = 1.3;
                std::clock_t t = clock();
                double sum_delta = 0.;
                double sum_2_delta = 0.;
                double D = exp(-r * T);
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    W = W * sig + mu; //Change of measure
                    double S_T_m = (S_0 - d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_h = (S_0 + d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double OptionPrice_m = D*std::max(S_T_m - K, 0.)*(1/sqrt(2*M_PI))*exp(-pow(W,2)/2)/((1/(sig*sqrt(2*M_PI)))*exp(-0.5*(pow((W-mu)/sig,2))));
                    double OptionPrice_h = D*std::max(S_T_h - K, 0.)*(1/sqrt(2*M_PI))*exp(-pow(W,2)/2)/((1/(sig*sqrt(2*M_PI)))*exp(-0.5*(pow((W-mu)/sig,2))));
                    sum_delta += (OptionPrice_h - OptionPrice_m) / (2 * d_S);
                    sum_2_delta += std::pow((OptionPrice_h - OptionPrice_m) / (2 * d_S), 2);
                }
                double m = (sum_delta) / (double) M;
                double v = ((1 / (double) M) * sum_2_delta - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_gamma(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double d_S) {
                double sig = 0.7;
                double mu = 1.3;
                std::clock_t t = clock();
                double sum_gamma = 0.;
                double sum_2_gamma = 0.;
                double D = exp(-r * T);
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    W = W * sig + mu; //Change of measure
                    double S_T_l = (S_0 - d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_m = S_0 * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double S_T_h = (S_0 + d_S) * exp(T * (r - 0.5 * sigma * sigma) + sigma * sqrt(T) * W);
                    double OptionPrice_l = D*std::max(S_T_l - K, 0.)*(1/sqrt(2*M_PI))*exp(-pow(W,2)/2)/((1/(sig*sqrt(2*M_PI)))*exp(-0.5*(pow((W-mu)/sig,2))));
                    double OptionPrice_m = D*std::max(S_T_m - K, 0.)*(1/sqrt(2*M_PI))*exp(-pow(W,2)/2)/((1/(sig*sqrt(2*M_PI)))*exp(-0.5*(pow((W-mu)/sig,2))));
                    double OptionPrice_h = D*std::max(S_T_h - K, 0.)*(1/sqrt(2*M_PI))*exp(-pow(W,2)/2)/((1/(sig*sqrt(2*M_PI)))*exp(-0.5*(pow((W-mu)/sig,2))));
                    sum_gamma += (OptionPrice_h - 2 * OptionPrice_m + OptionPrice_l) / (d_S * d_S);
                    sum_2_gamma += std::pow((OptionPrice_h - 2 * OptionPrice_m + OptionPrice_l) / (d_S * d_S), 2);
                }
                double m = (sum_gamma) / (double) M;
                double v = ((1 / (double) M) * sum_2_gamma - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }

            std::vector<double>
            call_vega(double &S_0, double &K, double &r, double &sigma, double &T, int &M, double d_sigma) {
                double sig = 0.7;
                double mu = 1.3;
                std::clock_t t = clock();
                double sum_vega = 0.;
                double sum_2_vega = 0.;
                double D = exp(-r * T);
                for (int i = 0; i < M; ++i) {
                    double X_1 = (double) rand() / RAND_MAX;
                    double X_2 = (double) rand() / RAND_MAX;
                    double W = sqrt(-2 * log(X_1)) * sin(2 * M_PI * X_2);
                    W = W * sig + mu; //Change of measure
                    double S_T_m = S_0 * exp(T * (r - 0.5 * (sigma - d_sigma) * (sigma - d_sigma)) + (sigma - d_sigma) * sqrt(T) * W);
                    double S_T_h = S_0 * exp(T * (r - 0.5 * (sigma + d_sigma) * (sigma + d_sigma)) + (sigma + d_sigma) * sqrt(T) * W);
                    double OptionPrice_m = D*std::max(S_T_m - K, 0.)*(1/sqrt(2*M_PI))*exp(-pow(W,2)/2)/((1/(sig*sqrt(2*M_PI)))*exp(-0.5*(pow((W-mu)/sig,2))));
                    double OptionPrice_h = D*std::max(S_T_h - K, 0.)*(1/sqrt(2*M_PI))*exp(-pow(W,2)/2)/((1/(sig*sqrt(2*M_PI)))*exp(-0.5*(pow((W-mu)/sig,2))));
                    sum_vega += (OptionPrice_h - OptionPrice_m) / (2*d_sigma);
                    sum_2_vega += std::pow((OptionPrice_h - OptionPrice_m) / d_sigma, 2);
                }
                double m = (sum_vega) / (double) M;
                double v = ((1 / (double) M) * sum_2_vega - m * m) / (double) M;
                double duration = (std::clock() - t) / (double) CLOCKS_PER_SEC;
                std::vector<double> res = {m, 1.96*sqrt(v / M), duration};
                return res;
            }
        }

    }
}