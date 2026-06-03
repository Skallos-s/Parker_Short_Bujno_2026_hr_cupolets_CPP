// File: cupolet_search.h
// Purpose: Cupolet search methods header file
// Author: Daniel Bujno

#ifndef __cupolet_search_h_
#define __cupolet_search_h_

#include "hindmarsh_rose.h"

#include <string>
#include <vector>

// Finds all cupolets with given ctrl sequence of length less than limit, that star ton PS1 plane.
// Returns vector of cupolets. Each cupolet is a vector of pairs {index, ps}.
std::vector<std::vector<std::vector<unsigned int>>> find_cupolets(std::vector<unsigned int> &ctrl, std::string bin_rn_direc, unsigned int limit = 600);

// Generates cupolet based off of control. Reads in neuron state, dt, number of bins and
// crossings, control sequence, number of iterations, directory to pull info, and ctrl0,ctrl1 arrays
// that tell how to implement macrocontrol. Returns time series of cupolet.
std::vector<std::vector<double>> cupolet_time_series(hindmarsh_rose &neuron, double dt, unsigned int bins, std::vector<std::vector<unsigned int>> cupolet, const std::string direc, const std::string bin_rn_direc);

#endif