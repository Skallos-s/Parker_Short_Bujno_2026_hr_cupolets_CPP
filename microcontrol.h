// File: microcontrol.h
// Purpose: Microcontrol header file
// Author: Daniel Bujno

#ifndef __microcontrol_h_
#define __microcontrol_h_

#include "hindmarsh_rose.h"

#include <string>
#include <vector>

// Reads in neuron object, direc directory, number of bins bins, and time step dt.
// Generates ps0_micro and ps1_micro that represents mapping for each bin and save to direc.
void establish_microcontrol(hindmarsh_rose &neuron, const std::string direc, const std::string bin_rn_direc, double dt = 1.0/128, unsigned int bins = 1600);

// Reads in neuron object, dt time step, number of bins bins, ps control plane, and direc directory. Returns corresponding PS control plane microcontrol mapping.
std::vector<std::vector<unsigned int>> micro_control_map(hindmarsh_rose &neuron, const std::string direc, const std::string bin_rn_direc, double dt, unsigned int bins, unsigned int ps);

#endif