#ifndef SM_PLOT_H
#define SM_PLOT_H
#include <vector>
#include <iostream>
#include <boost/tuple/tuple.hpp>
#include <boost/format.hpp>
#include "gnuplot-iostream/gnuplot-iostream.h"
#include "basic.h"
#include "heston.h"
#include "analytic.h"

enum SolvingType{
    Analytic,
    MC
};
enum Method{
    None,
    Direct,
    Antithetic,
    ControlVariate,
};
enum X_axis{
    Spot,
    Volatility,
    Maturity
};

class Plot{
public:
    Plot(X_axis X_type, double X_min, double X_max, int N): m_X_type(X_type), m_X_min(X_min), m_X_max(X_max), m_N(N) {}
    virtual ~Plot(){}
    virtual void set_values(SolvingType RES, Method MTD, double &S, double &K, double &r, double &sigma, double &T, int &M, double &d_S){};
    virtual void set_values(SolvingType RES, Method MTD, double &S, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho, double &T, int &M, double &delta_S, double &delta_sigma){};
    virtual void show() = 0;

    friend void show_both_price(Plot &p_GBM, Plot &p_Heston, double rho);
    friend void show_both_delta(Plot &p_GBM, Plot &p_Heston, double rho);
    friend void show_both_gamma(Plot &p_GBM, Plot &p_Heston, double rho);
    friend void show_both_vega(Plot &p_GBM, Plot &p_Heston, double rho);

protected:
    std::vector<std::pair<double, double> > m_pts_price;
    std::vector<std::pair<double, double> > m_pts_delta;
    std::vector<std::pair<double, double> > m_pts_gamma;
    std::vector<std::pair<double, double> > m_pts_vega;
    X_axis                                  m_X_type;
    double                                  m_X_min;
    double                                  m_X_max;
    int                                     m_N;
};

class Plot_GBM : public Plot {
public:
    Plot_GBM(X_axis X_type, double X_min, double X_max, int N): Plot(X_type, X_min, X_max, N) {}
    ~Plot_GBM(){}
    void set_values(SolvingType RES, Method MTD, double &S, double &K, double &r, double &sigma, double &T, int &M, double &d_S);
    void show();

};

class Plot_Heston : public Plot {
public:
    Plot_Heston(X_axis X_type, double X_min, double X_max, int N): Plot(X_type, X_min, X_max, N) {}
    ~Plot_Heston(){}
    void set_values(SolvingType RES, Method MTD, double &S, double &K, double &r, double &kappa, double &theta, double &sigma, double &rho, double &T, int &M, double &delta_S, double &delta_sigma);
    void show();

};

#endif //SM_PLOT_H