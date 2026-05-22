// File: rk4.h
// Purpose: Runge-Kutta 4th header file
// Author: Daniel Bujno

#ifndef __rk4_h_
#define __rk4_h_

#include "array3.h"
#include "hindmarsh_rose.h"

#include <vector>

// Runge-Kutta in a state x, time step dt, and dynamics function f.
// Performs autonomous RK4 method and returns next iteration
array3 rk4(array3 x, double dt, array3 (hindmarsh_rose::*f)(array3) const, hindmarsh_rose &hr);


// Perform Henon's trick with RK4 method. xi is the initial state after varying
// array to insert independent variable before the trick. ti is the state of Henon's
// independent variable. dd is the distance to travel and f is the dynamics for Henon's trick.
std::vector<double> rk4_henon(array3 xi, double ti, double dd, array3 (hindmarsh_rose::*f)(array3, double) const, hindmarsh_rose &hr);

#endif