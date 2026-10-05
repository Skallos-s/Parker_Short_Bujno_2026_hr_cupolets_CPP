// File: cupolet_from_impulse.h
// Purpose: Cupolet impulse methods header file
// Author: Daniel Bujno

#ifndef __cupolet_from_impulse_h_
#define __cupolet_from_impulse_h_

#include "hindmarsh_rose.h"

#include <string>
#include <vector>

// Generates impulse function for entire cupolet from control sequence and starting bin.
// Reads in neuron state, dt, number of bins starting bin, control sequence, and save directories.
// Series is time-aligned with multiples of dt.
void create_impulse_function(hindmarsh_rose &neuron, double dt, unsigned int bins, unsigned int start_bin, std::vector<unsigned int> &ctrl, const std::string direc, const std::string bin_rn_direc);

// Generates cupolet from impulse function.// Generates cupolet from starting bin and impulse series.
// Reads in neuron state, dt, number of bins starting bin, control sequence, and save directories.
void create_time_series_from_impulse(hindmarsh_rose &neuron, double dt, unsigned int bins, unsigned int start_bin, std::vector<unsigned int> &ctrl, const std::string direc, const std::string bin_rn_direc);

// Generates cupolet and impulse function from control sequence and starting bin. Runs stabilized cupolet for N times.
// Reads in neuron state, dt, number of bins starting bin, control sequence, and save directories.
// Series is time-aligned with multiples of dt.
void create_cupolet(hindmarsh_rose &neuron, double dt, unsigned int bins, unsigned int start_bin, unsigned int N, std::vector<unsigned int> &ctrl, const std::string direc, const std::string bin_rn_direc);


#endif