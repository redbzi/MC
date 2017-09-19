#include <iostream>
#include <iomanip>

#include "print.h"
#include "params.h"
#include "plot.h"
#include "mc.h"

void print_results_GBM() {
    std::cout << "\n" << std::endl;
    std::cout << std::setw(100) << "****************************************************************" << std::endl;
    std::cout << std::setw(100) << "****************       Black-Scholes model       ***************" << std::endl;
    std::cout << std::setw(100) << "****************************************************************" << std::endl;
    std::cout << "-- Analytic --" << std::endl;
    double a_price = GBM::Analytic::call_price(S_0, K, r, sigma, T);
    double a_delta = GBM::Analytic::call_delta(S_0, K, r, sigma, T);
    double a_gamma = GBM::Analytic::call_gamma(S_0, K, r, sigma, T);
    double a_vega = GBM::Analytic::call_vega(S_0, K, r, sigma, T);
    std::cout << " price: " << a_price <<
              "\n delta: " << a_delta <<
              "\n gamma: " << a_gamma <<
              "\n vega:  " << a_vega <<
              std::endl;


    std::cout << "\n--- MC-Basic Method" << std::endl;
    std::cout << std::setw(40) << "Value" << std::setw(30) << "Relative error" << std::setw(30) << "Half-Width"
              << std::setw(30) << "Computation Time" << std::endl;
    for (auto m : M) {
        std::vector<double> mc_price = GBM::MC::Basic::call_price(S_0, K, r, sigma, T, m);
        std::vector<double> mc_delta = GBM::MC::Basic::call_delta(S_0, K, r, sigma, T, m, delta_S);
        std::vector<double> mc_gamma = GBM::MC::Basic::call_gamma(S_0, K, r, sigma, T, m, delta_S);
        std::vector<double> mc_vega = GBM::MC::Basic::call_vega(S_0, K, r, sigma, T, m, delta_sigma);
        std::cout << "\n-- Simulations: " << m <<
                  "\n   price: " << std::setw(30) << mc_price[0] << std::setw(30) << (mc_price[0] - a_price) / a_price
                  << std::setw(30) << mc_price[1] << std::setw(30) << mc_price[2] <<
                  "\n   delta: " << std::setw(30) << mc_delta[0] << std::setw(30) << (mc_delta[0] - a_delta) / a_delta
                  << std::setw(30) << mc_delta[1] << std::setw(30) << mc_delta[2] <<
                  "\n   gamma: " << std::setw(30) << mc_gamma[0] << std::setw(30) << (mc_gamma[0] - a_gamma) / a_gamma
                  << std::setw(30) << mc_gamma[1] << std::setw(30) << mc_gamma[2] <<
                  "\n   vega:  " << std::setw(30) << mc_vega[0] << std::setw(30) << (mc_vega[0] - a_vega) / a_vega
                  << std::setw(30) << mc_vega[1] << std::setw(30) << mc_vega[2] <<
                  std::endl;
    }

    std::cout << "\n--- MC-Antithetic Method" << std::endl;
    std::cout << std::setw(40) << "Value" << std::setw(30) << "Relative error" << std::setw(30) << "Half-Width"
              << std::setw(30) << "Computation Time" << std::endl;
    for (auto m : M) {
        std::vector<double> mc_price = GBM::MC::Antithetic::call_price(S_0, K, r, sigma, T, m);
        std::vector<double> mc_delta = GBM::MC::Antithetic::call_delta(S_0, K, r, sigma, T, m, delta_S);
        std::vector<double> mc_gamma = GBM::MC::Antithetic::call_gamma(S_0, K, r, sigma, T, m, delta_S);
        std::vector<double> mc_vega = GBM::MC::Antithetic::call_vega(S_0, K, r, sigma, T, m, delta_sigma);
        std::cout << "\n-- Simulations: " << m <<
                  "\n   price: " << std::setw(30) << mc_price[0] << std::setw(30) << (mc_price[0] - a_price) / a_price
                  << std::setw(30) << mc_price[1] << std::setw(30) << mc_price[2] <<
                  "\n   delta: " << std::setw(30) << mc_delta[0] << std::setw(30) << (mc_delta[0] - a_delta) / a_delta
                  << std::setw(30) << mc_delta[1] << std::setw(30) << mc_delta[2] <<
                  "\n   gamma: " << std::setw(30) << mc_gamma[0] << std::setw(30) << (mc_gamma[0] - a_gamma) / a_gamma
                  << std::setw(30) << mc_gamma[1] << std::setw(30) << mc_gamma[2] <<
                  "\n   vega:  " << std::setw(30) << mc_vega[0] << std::setw(30) << (mc_vega[0] - a_vega) / a_vega
                  << std::setw(30) << mc_vega[1] << std::setw(30) << mc_vega[2] <<
                  std::endl;
    }

    std::cout << "\n--- MC-Importance Sampling" << std::endl;
    std::cout << std::setw(40) << "Value" << std::setw(30) << "Relative error" << std::setw(30) << "Half-Width"
              << std::setw(30) << "Computation Time" << std::endl;
    for (auto m : M) {
        std::vector<double> is_price = GBM::MC::ImportantSampling::call_price(S_0, K, r, sigma, T, m);
        std::vector<double> is_delta = GBM::MC::ImportantSampling::call_delta(S_0, K, r, sigma, T, m, delta_S);
        std::vector<double> is_gamma = GBM::MC::ImportantSampling::call_gamma(S_0, K, r, sigma, T, m, delta_S);
        std::vector<double> is_vega = GBM::MC::ImportantSampling::call_vega(S_0, K, r, sigma, T, m, delta_sigma);
        std::cout << "\n-- Simulations: " << m <<
                  "\n   price: " << std::setw(30) << is_price[0] << std::setw(30) << (is_price[0] - a_price) / a_price
                  << std::setw(30) << is_price[1] << std::setw(30) << is_price[2] <<
                  "\n   delta: " << std::setw(30) << is_delta[0] << std::setw(30) << (is_delta[0] - a_delta) / a_delta
                  << std::setw(30) << is_delta[1] << std::setw(30) << is_delta[2] <<
                  "\n   gamma: " << std::setw(30) << is_gamma[0] << std::setw(30) << (is_gamma[0] - a_gamma) / a_gamma
                  << std::setw(30) << is_gamma[1] << std::setw(30) << is_gamma[2] <<
                  "\n   vega:  " << std::setw(30) << is_vega[0] << std::setw(30) << (is_vega[0] - a_vega) / a_vega
                  << std::setw(30) << is_vega[1] << std::setw(30) << is_vega[2] <<
                  std::endl;
    }

    std::cout << "\n--- MC-Control Variate Method" << std::endl;
    std::cout << std::setw(40) << "Value" << std::setw(30) << "Relative error" << std::setw(30) << "Half-Width"
              << std::setw(30) << "Computation Time" << std::setw(30) << "Estimated correlation" << std::endl;
    for (auto m : M) {
        std::vector<double> mc_price = GBM::MC::ControlVariate::call_price(S_0, K, r, sigma, T, m);
        std::cout << "\n-- Simulations: " << m <<
                  "\n   price: " << std::setw(30) << mc_price[0] << std::setw(30) << (mc_price[0] - a_price) / a_price
                  << std::setw(30) << mc_price[1] << std::setw(30) << mc_price[2] << std::setw(30) << mc_price[3] <<
                  std::endl;
    }
}


void print_results_Heston(){
    std::cout << "\n" << std::endl;
    std::cout << std::setw(100) << "****************************************************************" << std::endl;
    std::cout << std::setw(100) << "****************          Heston model         *****************" << std::endl;
    std::cout << std::setw(100) << "****************************************************************" << std::endl;

    std::cout << "-- FFT --" << std::endl;
    std::cout << " price: " << FFT_price <<
              "\n delta: " << FFT_delta <<
              "\n vega_V_0:  " << FFT_vega_V_0 <<
              std::endl;

    std::cout << "\n-- Numerical Methods --" << std::endl;
    std::cout << " price: " << NM_price <<
              std::endl;

    std::cout << "\n--- MC-Basic Method" << std::endl;
    std::cout << std::setw(40) << "Value" << std::setw(30) << "Relative error" << std::setw(30) << "Half-Width" << std::setw(30) << "Computation Time" << std:: endl;
    for (auto m : M) {
        std::vector<double> mc_price = Heston::MC::Direct::call_price(S_0, V_0, K, r, kappa, theta, sigma, rho, T, m);
        std::vector<double> mc_delta = Heston::MC::Direct::call_delta(S_0, V_0, K, r, kappa, theta, sigma, rho, T, m,
                                                                      delta_S);
        std::vector<double> mc_gamma = Heston::MC::Direct::call_gamma(S_0, V_0, K, r, kappa, theta, sigma, rho, T, m,
                                                                      delta_S);
        std::vector<double> mc_vega_V_t = Heston::MC::Direct::call_vega_V_t(S_0, V_0, K, r, kappa, theta, sigma, rho, T,
                                                                            m, delta_sigma);
        std::vector<double> mc_vega_V_0 = Heston::MC::Direct::call_vega_V_0(S_0, V_0, K, r, kappa, theta, sigma, rho, T,
                                                                            m, delta_sigma);
        std::vector<double> mc_vega_theta = Heston::MC::Direct::call_vega_theta(S_0, V_0, K, r, kappa, theta, sigma,
                                                                                rho, T, m, delta_theta);
        std::cout << "\n-- Simulations: " << m <<
                  "\n   price: " << std::setw(30) << mc_price[0] << std::setw(30)
                  << (mc_price[0] - FFT_price) / FFT_price << std::setw(30) << mc_price[1] << std::setw(30)
                  << mc_price[2] <<
                  "\n   delta: " << std::setw(30) << mc_delta[0] << std::setw(30)
                  << (mc_delta[0] - FFT_delta) / FFT_delta << std::setw(30) << mc_delta[1] << std::setw(30)
                  << mc_delta[2] <<
                  "\n   gamma: " << std::setw(30) << mc_gamma[0] << std::setw(30) << "N/A" << std::setw(30)
                  << mc_gamma[1] << std::setw(30) << mc_gamma[2] <<
                  "\n   vega t:  " << std::setw(28) << mc_vega_V_t[0] << std::setw(30) << "N/A" << std::setw(30)
                  << mc_vega_V_t[1] << std::setw(30) << mc_vega_V_t[2] <<
                  "\n   vega 0:  " << std::setw(28) << mc_vega_V_0[0] << std::setw(30)
                  << (mc_vega_V_0[0] - FFT_vega_V_0) / FFT_vega_V_0 << std::setw(30) << mc_vega_V_0[1] << std::setw(30)
                  << mc_vega_V_0[2] <<
                  "\n   vega theta:  " << std::setw(24) << mc_vega_theta[0] << std::setw(30) << "N/A" << std::setw(30)
                  << mc_vega_theta[1] << std::setw(30) << mc_vega_theta[2] <<
                  std::endl;
    }

    std::cout << "\n--- MC-Antithetic Method" << std::endl;
    std::cout << std::setw(40) << "Value" << std::setw(30) << "Relative error" << std::setw(30) << "Half-Width" << std::setw(30) << "Computation Time" << std:: endl;
    for (auto m : M) {
        std::vector<double> mc_price = Heston::MC::Antithetic::call_price(S_0, V_0, K, r, kappa, theta, sigma, rho, T,
                                                                          m);
        std::vector<double> mc_delta = Heston::MC::Antithetic::call_delta(S_0, V_0, K, r, kappa, theta, sigma, rho, T,
                                                                          m, delta_S);
        std::vector<double> mc_gamma = Heston::MC::Antithetic::call_gamma(S_0, V_0, K, r, kappa, theta, sigma, rho, T,
                                                                          m, delta_S);
        std::vector<double> mc_vega = Heston::MC::Antithetic::call_vega_V_t(S_0, V_0, K, r, kappa, theta, sigma, rho, T,
                                                                            m, delta_sigma);
        std::cout << "\n-- Simulations: " << m <<
                  "\n   price: " << std::setw(30) << mc_price[0] << std::setw(30) << (mc_price[0] - FFT_price) / FFT_price << std::setw(30) << mc_price[1] << std::setw(30) << mc_price[2] <<
                  "\n   delta: " << std::setw(30) << mc_delta[0] << std::setw(30) << (mc_delta[0] - FFT_delta) / FFT_delta << std::setw(30) << mc_delta[1] << std::setw(30) << mc_delta[2] <<
                  "\n   gamma: " << std::setw(30) << mc_gamma[0] << std::setw(30) << "N/A" << std::setw(30)
                  << mc_gamma[1] << std::setw(30) << mc_gamma[2] <<
                  "\n   vega t:  " << std::setw(28) << mc_vega[0] << std::setw(30) << "N/A" << std::setw(30)
                  << mc_vega[1] << std::setw(30) << mc_vega[2] <<
                  std::endl;
    }

    std::cout << "\n--- MC-Control Variate Method" << std::endl;
    std::cout << std::setw(40) << "Value" << std::setw(30) << "Relative error" << std::setw(30) << "Half-Width" << std::setw(30) << "Computation Time" << std::setw(30) << "Estimated correlation" << std:: endl;
    for (auto m : M) {
        std::vector<double> mc_price = Heston::MC::ControlVariate::call_price(S_0, V_0, K, r, kappa, theta, sigma, rho,
                                                                              T, m);
        std::cout << "\n-- Simulations: " << m <<
                  "\n   price: " << std::setw(30) << mc_price[0] << std::setw(30) << (mc_price[0] - FFT_price) / FFT_price << std::setw(30) << mc_price[1] << std::setw(30) << mc_price[2] << std::setw(30) << mc_price[3] <<
                  std::endl;
    }
}


void print_graphs_GBM(){
    Plot *plt_1 = new Plot_GBM(Spot, S_range[0], S_range[1], N);
    plt_1->set_values(Analytic, None, S_0, K, r, sigma, T, M[0], delta_S);
    plt_1->show();
    delete plt_1;

    Plot *plt_2 = new Plot_GBM(Volatility, V_range[0], V_range[1], N);
    plt_2->set_values(Analytic, None, S_0, K, r, sigma, T, M[0], delta_S);
    plt_2->show();
    delete plt_2;

    Plot *plt_3 = new Plot_GBM(Maturity, T_range[0], T_range[1], N);
    plt_3->set_values(Analytic, None, S_0, K, r, sigma, T, M[0], delta_S);
    plt_3->show();
    delete plt_3;

    Plot *plt_4 = new Plot_GBM(Spot, S_range[0], S_range[1], N);
    plt_4->set_values(MC, Antithetic, S_0, K, r, sigma, T, M[2], delta_S);
    plt_4->show();
    delete plt_4;

    Plot *plt_5 = new Plot_GBM(Volatility, V_range[0], V_range[1], N);
    plt_5->set_values(MC, Antithetic, S_0, K, r, sigma, T, M[2], delta_S);
    plt_5->show();
    delete plt_5;

    Plot *plt_6 = new Plot_GBM(Maturity, T_range[0], T_range[1], N);
    plt_6->set_values(MC, Antithetic, S_0, K, r, sigma, T, M[2], delta_S);
    plt_6->show();
    delete plt_6;
}

void print_graphs_Heston(){
    Plot *plt_7 = new Plot_Heston(Spot, S_range[0], S_range[1], N);
    plt_7->set_values(MC, Antithetic, S_0, V_0, K, r, kappa, theta, sigma, rho, T, M[1], delta_S, delta_sigma);
    plt_7->show();
    delete plt_7;

    Plot *plt_8 = new Plot_Heston(Maturity, T_range[0], T_range[1], N);
    plt_8->set_values(MC, Antithetic, S_0, V_0, K, r, kappa, theta, sigma, rho, T, M[1], delta_S, delta_sigma);
    plt_8->show();
    delete plt_8;
}

void print_graphs_both(){
    Plot *plt_1 = new Plot_GBM(Spot, S_range[0], S_range[1], N);
    plt_1->set_values(Analytic, None, S_0, K, r, sigma, T, M[0], delta_S);

    for (auto rho : rho_range) {
        Plot *plt_2 = new Plot_Heston(Spot, S_range[0], S_range[1], N);
        plt_2->set_values(MC, Antithetic, S_0, V_0, K, r, kappa, theta, sigma, rho, T, M[1], delta_S, delta_sigma);

        show_both_price(*plt_1, *plt_2, rho);
        show_both_delta(*plt_1, *plt_2, rho);
        show_both_gamma(*plt_1, *plt_2, rho);
        show_both_vega(*plt_1, *plt_2, rho);
        if (rho == rho_range_size) {
            delete plt_2;
        }
    }
    delete plt_1;
}


