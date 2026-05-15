// File: rk4.cpp
// Purpose: Implementation file for the Runge-Kutta 4th methods
// Author: Daniel Bujno

#include "rk4.h"
#include "array3.h"
#include "hindmarsh_rose.h"

#include <vector>

// Runge-Kutta in a state x, time step dt, and dynamics function f.
// Performs autonomous RK4 method and returns next iteration
array3 rk4(array3 x, double dt, array3 (*f)(array3)) {
	array3 k1 = dt * f(x);
	array3 k2 = dt * f(x + k1 / 2);
	array3 k3 = dt * f(x + k2 / 2);
	array3 k4 = dt * f(x + k3);
	return x + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
}

// Perform Henon's trick with RK4 method. xi is the initial state after varying
// array to insert independent variable before the trick. ti is the state of Henon's
// independent variable. dd is the distance to travel and f is the dynamics for Henon's trick.
std::vector<double> rk4_henon(array3 xi, double ti, double dd, array3 (hindmarsh_rose::*f)(array3, double) const, hindmarsh_rose hr) {
	array3 k1 = dd * (hr.*f)(xi         , ti);
	array3 k2 = dd * (hr.*f)(xi + k1 / 2, ti + dd / 2);
	array3 k3 = dd * (hr.*f)(xi + k2 / 2, ti + dd / 2);
	array3 k4 = dd * (hr.*f)(xi + k3    , ti + dd);
	array3 update = xi + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
	std::vector<double> output = {update.get(0), update.get(1), update.get(2), ti + dd};
	return output;
}