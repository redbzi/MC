#ifndef SM_INPUT_PARAMETERS_H
#define SM_INPUT_PARAMETERS_H


////////  Financial parameters ///////
double S_0 = 100.;
double K = 100. ;
double r = 0.025;
double sigma = 0.3;
double T = 1.;
double kappa = 1.5;
double theta = 0.04;
double rho = -0.9;

//////// Animation ////////
double rho_range [11] = {-1., -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1.};


//////// Compute the greeks ///////
double delta_S = 0.001 * S_0;
double delta_sigma = 0.001 * sigma;
double delta_theta = 0.001 * theta;

//////// Number of simulations ///////
int M [3] = {100, 10000, 100000};
int N = 100;

//////// solutions with Fast Fourier Transform ///////
double FFT_price    = 8.8948693600540167;
double FFT_delta    = 0.671558511602186;
double FFT_vega_V_0 = 19.59866330837231;

//////// solutions from Numerical Methods ///////
double NM_price = 8.894628;

//////// Range of X-axis for the graphs ///////
double S_range [2] = {1., 200.};
double V_range [2] = {0., 1.};
double T_range [2] = {1., 10.};



#endif //SM_INPUT_PARAMETERS_H
