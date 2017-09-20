#include "plot.h"
#include "mc.h"


void Plot_GBM::set_values(SolvingType RES, Method MTD, double &S, double &K, double &r, double &sigma, double &T, int &M, double &d_S){
    switch (RES){
        case MC:
        switch (MTD) {
            case Direct:
                switch(m_X_type){
                    case Spot: {
                        double X = m_X_min;
                        while (X < m_X_max) {
                            m_pts_price.emplace_back(
                                    std::make_pair(X, GBM::MC::Basic::call_price(X, K, r, sigma, T, M)[0]));
                            m_pts_delta.emplace_back(
                                    std::make_pair(X, 100 * GBM::MC::Basic::call_delta(X, K, r, sigma, T, M, d_S)[0]));
                            m_pts_gamma.emplace_back(
                                    std::make_pair(X, 5000 * GBM::MC::Basic::call_gamma(X, K, r, sigma, T, M, d_S)[0]));
                            m_pts_vega.emplace_back(
                                    std::make_pair(X, GBM::MC::Basic::call_vega(X, K, r, sigma, T, M, d_S)[0]));

                            X += (m_X_max - m_X_min) / (double) m_N;
                        }
                        break;
                    }
                    case Volatility: {
                        double X = m_X_min;
                        while (X < m_X_max) {
                            m_pts_price.emplace_back(
                                    std::make_pair(X, GBM::MC::Basic::call_price(S, K, r, X, T, M)[0]));
                            m_pts_delta.emplace_back(
                                    std::make_pair(X, 100 * GBM::MC::Basic::call_delta(S, K, r, X, T, M, d_S)[0]));
                            m_pts_gamma.emplace_back(
                                    std::make_pair(X, 5000 * GBM::MC::Basic::call_gamma(S, K, r, X, T, M, d_S)[0]));
                            m_pts_vega.emplace_back(
                                    std::make_pair(X, GBM::MC::Basic::call_vega(S, K, r, X, T, M, d_S)[0]));

                            X += (m_X_max - m_X_min) / (double) m_N;
                        }
                        break;
                    }
                    case Maturity: {
                        double X = m_X_min;
                        while (X < m_X_max) {
                            m_pts_price.emplace_back(
                                    std::make_pair(X, GBM::MC::Basic::call_price(S, K, r, sigma, X, M)[0]));
                            m_pts_delta.emplace_back(
                                    std::make_pair(X, 100 * GBM::MC::Basic::call_delta(S, K, r, sigma, X, M, d_S)[0]));
                            m_pts_gamma.emplace_back(
                                    std::make_pair(X, 5000 * GBM::MC::Basic::call_gamma(S, K, r, sigma, X, M, d_S)[0]));
                            m_pts_vega.emplace_back(
                                    std::make_pair(X, GBM::MC::Basic::call_vega(S, K, r, sigma, X, M, d_S)[0]));

                            X += (m_X_max - m_X_min) / (double) m_N;
                        }
                        break;
                    }
                    default:
                        std::cout<< "X_axis not implemented" << std::endl;
                }
                break;
            case Antithetic:
                switch(m_X_type){
                    case Spot: {
                        double X = m_X_min;
                        while (X < m_X_max) {
                            m_pts_price.emplace_back(
                                    std::make_pair(X, GBM::MC::Antithetic::call_price(X, K, r, sigma, T, M)[0]));
                            m_pts_delta.emplace_back(std::make_pair(X, 100 *
                                                                       GBM::MC::Antithetic::call_delta(X, K, r, sigma,
                                                                                                       T, M, d_S)[0]));
                            m_pts_gamma.emplace_back(std::make_pair(X, 5000 *
                                                                       GBM::MC::Antithetic::call_gamma(X, K, r, sigma,
                                                                                                       T, M, d_S)[0]));
                            m_pts_vega.emplace_back(
                                    std::make_pair(X, GBM::MC::Antithetic::call_vega(X, K, r, sigma, T, M, d_S)[0]));

                            X += (m_X_max - m_X_min) / (double) m_N;
                        }
                        break;
                    }
                    case Volatility: {
                        double X = m_X_min;
                        while (X < m_X_max) {
                            m_pts_price.emplace_back(
                                    std::make_pair(X, GBM::MC::Antithetic::call_price(S, K, r, X, T, M)[0]));
                            m_pts_delta.emplace_back(
                                    std::make_pair(X, 100 * GBM::MC::Antithetic::call_delta(S, K, r, X, T, M, d_S)[0]));
                            m_pts_gamma.emplace_back(std::make_pair(X, 5000 *
                                                                       GBM::MC::Antithetic::call_gamma(S, K, r, X, T, M,
                                                                                                       d_S)[0]));
                            m_pts_vega.emplace_back(
                                    std::make_pair(X, GBM::MC::Antithetic::call_vega(S, K, r, X, T, M, d_S)[0]));

                            X += (m_X_max - m_X_min) / (double) m_N;
                        }
                        break;
                    }

                    case Maturity: {
                        double X = m_X_min;
                        while (X < m_X_max) {
                            m_pts_price.emplace_back(
                                    std::make_pair(X, GBM::MC::Antithetic::call_price(S, K, r, sigma, X, M)[0]));
                            m_pts_delta.emplace_back(std::make_pair(X, 100 *
                                                                       GBM::MC::Antithetic::call_delta(S, K, r, sigma,
                                                                                                       X, M, d_S)[0]));
                            m_pts_gamma.emplace_back(std::make_pair(X, 5000 *
                                                                       GBM::MC::Antithetic::call_gamma(S, K, r, sigma,
                                                                                                       X, M, d_S)[0]));
                            m_pts_vega.emplace_back(
                                    std::make_pair(X, GBM::MC::Antithetic::call_vega(S, K, r, sigma, X, M, d_S)[0]));

                            X += (m_X_max - m_X_min) / (double) m_N;
                        }
                        break;
                    }

                    default:
                        std::cout<< "X_axis not implemented" << std::endl;
                }
                break;
            default:
                std::cout<< "Monte-Carlo Reduction Variance Method not defined for plotting" << std::endl;

            }
            break;
        case Analytic:
            switch(m_X_type){
                case Spot: {
                    double X = m_X_min;
                    while (X < m_X_max) {
                        m_pts_price.emplace_back(std::make_pair(X, GBM::Analytic::call_price(X, K, r, sigma, T)));
                        m_pts_delta.emplace_back(std::make_pair(X, 100 * GBM::Analytic::call_delta(X, K, r, sigma, T)));
                        m_pts_gamma.emplace_back(
                                std::make_pair(X, 5000 * GBM::Analytic::call_gamma(X, K, r, sigma, T)));
                        m_pts_vega.emplace_back(std::make_pair(X, GBM::Analytic::call_vega(X, K, r, sigma, T)));

                        X += (m_X_max - m_X_min) / (double) m_N;
                    }
                    break;
                }

                case Volatility: {
                    double X = m_X_min;
                    while (X < m_X_max) {
                        m_pts_price.emplace_back(std::make_pair(X, GBM::Analytic::call_price(S, K, r, X, T)));
                        m_pts_delta.emplace_back(std::make_pair(X, 100 * GBM::Analytic::call_delta(S, K, r, X, T)));
                        m_pts_gamma.emplace_back(std::make_pair(X, 5000 * GBM::Analytic::call_gamma(S, K, r, X, T)));
                        m_pts_vega.emplace_back(std::make_pair(X, GBM::Analytic::call_vega(S, K, r, X, T)));

                        X += (m_X_max - m_X_min) / (double) m_N;
                    }
                    break;
                }

                case Maturity: {
                    double X = m_X_min;
                    while (X < m_X_max) {
                        m_pts_price.emplace_back(std::make_pair(X, GBM::Analytic::call_price(S, K, r, sigma, X)));
                        m_pts_delta.emplace_back(std::make_pair(X, 100 * GBM::Analytic::call_delta(S, K, r, sigma, X)));
                        m_pts_gamma.emplace_back(
                                std::make_pair(X, 5000 * GBM::Analytic::call_gamma(S, K, r, sigma, X)));
                        m_pts_vega.emplace_back(std::make_pair(X, GBM::Analytic::call_vega(S, K, r, sigma, X)));

                        X += (m_X_max - m_X_min) / (double) m_N;
                    }
                    break;
                }
                default:
                    std::cout<< "X_axis not implemented" << std::endl;
            }
            break;
        default:
            std::cout<< "SolvingType not defined for plotting" << std::endl;
    }
}

void Plot_GBM::show(){
    Gnuplot gp;
    gp << boost::format("set xrange [%1%:%2%]\nset yrange [0:100]\n")  % m_X_min % m_X_max;
    gp << "plot" << gp.file1d(m_pts_price) << "with points title 'GBM-Price'," <<
                    gp.file1d(m_pts_delta) << "with points title 'GBM-Delta'," <<
                    gp.file1d(m_pts_gamma) << "with points title 'GBM-Gamma'," <<
                    gp.file1d(m_pts_vega) << "with points title 'GBM-Vega'," << std::endl;
}


void Plot_Heston::set_values(SolvingType RES, Method MTD, double &S, double &V, double &K, double &r, double &kappa,
                             double &theta, double &sigma, double &rho, double &T, int &M, double &delta_S, double &delta_sigma){
    switch (RES){
        case MC:
            switch (MTD) {
                case Direct:
                    switch(m_X_type){
                        case Spot: {
                            double X = m_X_min;
                            while (X < m_X_max) {
                                m_pts_price.emplace_back(std::make_pair(X, Heston::MC::Direct::call_price(X, V, K, r,
                                                                                                          kappa,
                                                                                                          theta, sigma,
                                                                                                          rho,
                                                                                                          T, M)[0]));
                                m_pts_delta.emplace_back(
                                        std::make_pair(X, 100 * Heston::MC::Direct::call_delta(X, V, K, r,
                                                                                               kappa,
                                                                                               theta,
                                                                                               sigma, rho,
                                                                                               T, M,
                                                                                               delta_S)[0]));
                                m_pts_gamma.emplace_back(std::make_pair(X, 5000 *
                                                                           Heston::MC::Direct::call_gamma(X, V, K, r,
                                                                                                          kappa, theta,
                                                                                                          sigma, rho, T,
                                                                                                          M,
                                                                                                          delta_S)[0]));
                                m_pts_vega.emplace_back(std::make_pair(X,
                                                                       Heston::MC::Direct::call_vega_V_t(X, V, K, r,
                                                                                                         kappa,
                                                                                                         theta, sigma,
                                                                                                         rho,
                                                                                                         T, M,
                                                                                                         delta_sigma)[0]));

                                X += (m_X_max - m_X_min) / (double) m_N;
                            }
                            break;

                        }

                        case Maturity: {
                            double X = m_X_min;
                            while (X < m_X_max) {
                                m_pts_price.emplace_back(std::make_pair(X,
                                                                        Heston::MC::Direct::call_price(S, V, K, r,
                                                                                                       kappa,
                                                                                                       theta, sigma,
                                                                                                       rho,
                                                                                                       X, M)[0]));
                                m_pts_delta.emplace_back(
                                        std::make_pair(X, 100 * Heston::MC::Direct::call_delta(S, V, K, r,
                                                                                               kappa,
                                                                                               theta,
                                                                                               sigma, rho,
                                                                                               X, M,
                                                                                               delta_S)[0]));
                                m_pts_gamma.emplace_back(std::make_pair(X, 5000 *
                                                                           Heston::MC::Direct::call_gamma(S, V, K, r,
                                                                                                          kappa, theta,
                                                                                                          sigma, rho, X,
                                                                                                          M,
                                                                                                          delta_S)[0]));
                                m_pts_vega.emplace_back(std::make_pair(X,
                                                                       Heston::MC::Direct::call_vega_V_t(S, V, K, r,
                                                                                                         kappa,
                                                                                                         theta, sigma,
                                                                                                         rho,
                                                                                                         X, M,
                                                                                                         delta_sigma)[0]));

                                X += (m_X_max - m_X_min) / (double) m_N;
                            }
                            break;
                        }
                        default:
                            std::cout<< "X_axis not implemented" << std::endl;
                    }
                    break;
                case Antithetic:
                    switch(m_X_type){
                        case Spot: {
                            double X = m_X_min;
                            while (X < m_X_max) {
                                m_pts_price.emplace_back(
                                        std::make_pair(X, Heston::MC::Antithetic::call_price(X, V, K, r,
                                                                                             kappa, theta,
                                                                                             sigma, rho,
                                                                                             T, M)[0]));
                                m_pts_delta.emplace_back(std::make_pair(X, 100 *
                                                                           Heston::MC::Antithetic::call_delta(X, V, K,
                                                                                                              r,
                                                                                                              kappa,
                                                                                                              theta,
                                                                                                              sigma,
                                                                                                              rho,
                                                                                                              T, M,
                                                                                                              delta_S)[0]));
                                m_pts_gamma.emplace_back(std::make_pair(X, 5000 *
                                                                           Heston::MC::Antithetic::call_gamma(X, V, K,
                                                                                                              r,
                                                                                                              kappa,
                                                                                                              theta,
                                                                                                              sigma,
                                                                                                              rho,
                                                                                                              T, M,
                                                                                                              delta_S)[0]));
                                m_pts_vega.emplace_back(
                                        std::make_pair(X, Heston::MC::Antithetic::call_vega_V_t(X, V, K, r,
                                                                                                kappa,
                                                                                                theta,
                                                                                                sigma, rho,
                                                                                                T, M,
                                                                                                delta_sigma)[0]));

                                X += (m_X_max - m_X_min) / (double) m_N;
                            }
                            break;
                        }

                        case Maturity: {
                            double X = m_X_min;
                            while (X < m_X_max) {
                                m_pts_price.emplace_back(std::make_pair(X,
                                                                        Heston::MC::Direct::call_price(S, V, K, r,
                                                                                                       kappa,
                                                                                                       theta, sigma,
                                                                                                       rho,
                                                                                                       X, M)[0]));
                                m_pts_delta.emplace_back(
                                        std::make_pair(X, 100 * Heston::MC::Direct::call_delta(S, V, K, r,
                                                                                               kappa,
                                                                                               theta,
                                                                                               sigma, rho,
                                                                                               X, M,
                                                                                               delta_S)[0]));
                                m_pts_gamma.emplace_back(std::make_pair(X, 5000 *
                                                                           Heston::MC::Direct::call_gamma(S, V, K, r,
                                                                                                          kappa, theta,
                                                                                                          sigma, rho, X,
                                                                                                          M,
                                                                                                          delta_S)[0]));
                                m_pts_vega.emplace_back(std::make_pair(X,
                                                                       Heston::MC::Direct::call_vega_V_t(S, V, K, r,
                                                                                                         kappa,
                                                                                                         theta, sigma,
                                                                                                         rho,
                                                                                                         X, M,
                                                                                                         delta_sigma)[0]));

                                X += (m_X_max - m_X_min) / (double) m_N;
                            }
                            break;
                        }

                        default:
                            std::cout<< "X_axis not implemented" << std::endl;
                    }
                    break;
                default:
                    std::cout<< "Monte-Carlo Reduction Variance Method not defined for plotting" << std::endl;
            }
            break;
        default:
            std::cout<< "SolvingType not defined for plotting" << std::endl;
    }
}

void Plot_Heston::show(){
    Gnuplot gp;
    gp << boost::format("set xrange [%1%:%2%]\nset yrange [0:100]\n")  % m_X_min % m_X_max;
    gp << "plot" << gp.file1d(m_pts_price) << "with points title 'Heston-Price'," <<
       gp.file1d(m_pts_delta) << "with points title 'Heston-Delta'," <<
       gp.file1d(m_pts_gamma) << "with points title 'Heston-Gamma'," <<
       gp.file1d(m_pts_vega) << "with points title 'Heston-Vega'," << std::endl;
}


void show_both_price(Plot &p_GBM, Plot &p_Heston, double rho) {
    Gnuplot gp;
    gp << boost::format("set xrange [%1%:%2%]\nset yrange [0:110]\n set xlabel 'Spot'\n") % p_GBM.m_X_min % p_GBM.m_X_max;
    gp << "plot" <<
       gp.file1d(p_GBM.m_pts_price) << "with lines lc rgb \"black\" title 'price (BS)'," <<
       gp.file1d(p_Heston.m_pts_price) << boost::format("with dots lw 5 lc rgb \"red\" title 'price (Heston), rho=%1%'") %rho << std::endl;
}

void show_both_delta(Plot &p_GBM, Plot &p_Heston, double rho){
    Gnuplot gp;
    gp << boost::format("set xrange [%1%:%2%]\nset yrange [0:110]\n set xlabel 'Spot'\n") % p_GBM.m_X_min % p_GBM.m_X_max;
    gp << "plot" <<
       gp.file1d(p_GBM.m_pts_delta) << "with lines lc rgb \"black\" title 'delta (BS)'," <<
       gp.file1d(p_Heston.m_pts_delta) << boost::format("with dots lw 5 lc rgb \"red\" title 'delta (Heston), rho=%1%'") %rho << std::endl;
}

void show_both_gamma(Plot &p_GBM, Plot &p_Heston, double rho){
    Gnuplot gp;
    gp << boost::format("set xrange [%1%:%2%]\nset yrange [0:130]\n set xlabel 'Spot'\n") % p_GBM.m_X_min % p_GBM.m_X_max;
    gp << "plot" <<
       gp.file1d(p_GBM.m_pts_gamma) << "with lines lc rgb \"black\" title 'gamma (BS)'," <<
       gp.file1d(p_Heston.m_pts_gamma) << boost::format("with dots lw 5 lc rgb \"red\" title 'gamma (Heston), rho=%1%'") %rho << std::endl;
}

void show_both_vega(Plot &p_GBM, Plot &p_Heston, double rho){
    Gnuplot gp;
    gp << boost::format("set xrange [%1%:%2%]\nset yrange [0:100]\n set xlabel 'Spot'\n") % p_GBM.m_X_min % p_GBM.m_X_max;
    gp << "plot" <<
       gp.file1d(p_GBM.m_pts_vega)  << "with lines lc rgb \"black\" title 'vega (BS)'," <<
       gp.file1d(p_Heston.m_pts_vega)  << boost::format("with dots lw 5 lc rgb \"red\" title 'vega (Heston), rho=%1%'") %rho << std::endl;
}
