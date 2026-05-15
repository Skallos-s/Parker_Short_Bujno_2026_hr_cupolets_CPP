// File: helper.h
// Purpose: Numpy and Scipy methods header file
// Author: Daniel Bujno

#ifndef __helper_h_
#define __helper_h_

#include "array3.h"

#include <string>
#include <vector>

// Load data from text file
void loadtxt(std::vector<double> &vec_t, std::vector<double> &vec_x, std::vector<double> &vec_y, std::vector<double> &vec_z, const std::string file_path);

// Find local maxima of vector
// Returns vector of indices
std::vector<unsigned int> find_peaks(std::vector<double> &vec);

// Like find_peaks, but negates input and finds peaks within [0,1]
std::vector<unsigned int> find_mid_peaks(std::vector<double> &vec);

// Write 2D array to file
// std::vector<array3> data has only 2 array3s
// std::vector<std::vector<double>> data contains vector<double> of length 5 (t x y z ps)
void save_array3_vector(std::vector<array3> data, const std::string header, const std::string location);
void save_vector_vector(std::vector<std::vector<double>> data, const std::string header, const std::string location);

// Signum function
// From https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

#endif