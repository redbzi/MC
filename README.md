# MC-Heston


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
- gnuplot-iostream interface http://stahlke.org/dan/gnuplot-iostream/
- boost > 1.63.0 library http://www.boost.org/users/download/
- cmake > 3.8 https://cmake.org/download/
