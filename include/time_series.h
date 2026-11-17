// File: time_series.cpp
// Purpose: Time Series header file
// Author: Daniel Bujno

#ifndef __time_series_h_
#define __time_series_h_

#include "hindmarsh_rose.h"

// Generate time series data of hr with step size dt on the interval [0,t] with initial conditions IC.
// Save to direc/hr_system/hr_time_series.txt
void generate_time_series(hindmarsh_rose hr, array3 IC, const std::string direc, double dt = 1.0/128, double t_max = 10000);

// Generate time series data of hr with step size dt for N iterations with initial conditions IC.
// Save to direc/hr_system/hr_time_series.txt
void generate_time_series(hindmarsh_rose hr, array3 IC, const std::string direc, double dt = 1.0/128, unsigned int N = 1280000);

#endif