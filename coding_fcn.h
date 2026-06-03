// File: coding_fcn.h
// Purpose: Coding Functions methods header file
// Author: Daniel Bujno

#ifndef __coding_fcn_h_
#define __coding_fcn_h_


#include "hindmarsh_rose.h"

#include <string>
#include <vector>

// Read in HR neuron, dt, direc directory to store files, the number of iterations for the coding
// function to store and the number of bins to parse each poincare surface
void coding_fcn(hindmarsh_rose &neuron, const std::string direc, const std::string bin_rn_direc, double dt = 1.0/128, unsigned int coding_fcn_N = 16, unsigned int bins = 1600, bool generate_poly = false);

// For surface PS reads in x vertex values psx, y vertex values psy and polynmoial fit poly_coords,
// and bins to generate the bin spacing, midpoints, and initial data for evaluation 
void gen_bin_data(std::vector<double> &rbin_endpts, std::vector<double> &rbin_mids, std::vector<array3> &rbin_inits, std::vector<double> &psx, std::vector<double> &psy, std::vector<double> &poly_coords, unsigned int bins, int ps);

// Reads in list of initial bins x, crossings N, number of bins, dt, vertex data of PS0 and PS1 (ps0x, ps0y, ps1x, ps1y)
// Assigns decimal value of binary visitation sequence and corresponding symbolic dynamics for each bin to map and rmap
// map is a vector of length (bins) containing vectors of length (N)
// rmap is a vector of length (bins)
// map will contain the binary sequence of the next N control planes
// rmap will contain the binary decimal that map produces
void gen_crossing_sequence(std::vector<double> &rmap, std::vector<std::vector<unsigned int>> &map, std::vector<array3> &x, hindmarsh_rose &neuron, std::vector<double> &ps0x, std::vector<double> &ps0y, std::vector<double> &ps1x, std::vector<double> &ps1y, unsigned int N = 16, unsigned int bins = 1600, double dt = 1.0/128);

// Returns polynomial coefficients for quadratic/cubic fit
// 3x3 and 4x4 matrix inverses are specialized
std::vector<double> quadratic_regression(std::vector<double> &ind, std::vector<double> &dep);
std::vector<double> cubic_regression(std::vector<double> &ind, std::vector<double> &dep);




#endif