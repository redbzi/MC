#ifndef SM_ANALYTIC_H
#define SM_ANALYTIC_H
#include<cmath>

namespace GBM {
    namespace Analytic {
        double norm_pdf(const double &x);

        double norm_cdf(const double &x);

        double d_j(const int &j, const double &S, const double &K, const double &r, const double &v, const double &T);

        double call_price(const double &S, const double &K, const double &r, const double &v, const double &T);

        double call_delta(const double &S, const double &K, const double &r, const double &v, const double &T);

        double call_gamma(const double &S, const double &K, const double &r, const double &v, const double &T);

        double call_vega(const double &S, const double &K, const double &r, const double &v, const double &T);
    }
}
#endif //SM_ANALYTIC_H
