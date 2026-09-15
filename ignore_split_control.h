// File: split_control.h
// Purpose: Split control planes transfer header file
// Author: Daniel Bujno

#ifndef __split_contol_h_
#define __split_contol_h_

#include "hindmarsh_rose.h"

#include <string>
#include <vector>

// Reads in neuron object, direc directory, number of bins bins, and time step dt.
// Generates the initial points of the PS2b bins.
void establish_split_control_plane_bins(hindmarsh_rose &neuron, const std::string direc, const std::string bin_rn_direc, double dt = 1.0/128, unsigned int bins = 1600);

#endif