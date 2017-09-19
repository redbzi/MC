# MC


## Description ##
In this project we have implemented some Monte-Carlo simulations around the Black-Scholes and Heston models.
Given some parameters of an european call option, this code computes for {100, 10000, 100000} simulations the following quantities:
- price
- delta
- gamma
- vega

Some reduction variance methods have been implemented:
- Antithetic Variables
- Control Variate
- Importance Sampling
All of them show great improvements in the estimator accuracy.


We have also plotted the greeks and the price on a same scale against the spot, maturity and volatility (only for the Black-Scholes model in the latter).


## Requirements ##
- gnuplot > 5.2 (tested with qt terminal) http://www.gnuplot.info/download.html
- gnuplot-iostream interface http://stahlke.org/dan/gnuplot-iostream/
- boost > 1.63.0 library http://www.boost.org/users/download/
- cmake > 3.8 https://cmake.org/download/

## Build and run ##
1. clone the project: git clone https://github.com/RedwanBouizi/MC-Heston.git
2. go to MC directory: cd MC
3. cmake
4. make

Both a library and an executable will be build. To run the latter, follow these commands:
1. cd cmake-build-debug/tests/
2. ./run -1 -2 -3

Option 1: performs pricing and greeks computation using both the Black-Scholes and Heston models.

Option 2: outputs several graphs telling about the influence of spot, maturity, volatility over the price and greeks computations.

Option 3: shows how the asymmetry (price and greeks) evolves as the correlation Heston model increases.

