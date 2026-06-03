// File: helper.h
// Purpose: Numpy and Scipy methods header file
// Author: Daniel Bujno

#ifndef __helper_h_
#define __helper_h_

#include "array3.h"

#include <string>
#include <vector>

// Load data from text file
void loadtxt_4(std::vector<double> &vec_t, std::vector<double> &vec_x, std::vector<double> &vec_y, std::vector<double> &vec_z, const std::string file_path);
void loadtxt_3(std::vector<double> &psx, std::vector<double> &psy, std::vector<double> &psz, const std::string file_path);
void loadtxt_2(std::vector<unsigned int> &index, std::vector<unsigned int> &ps, const std::string file_path);
void loadtxt_1(std::vector<std::vector<double>> &pts, const std::string file_path);
void loadtxt_1(std::vector<double> &pts, const std::string file_path);
void loadtxt_1(std::vector<array3> &pts, const std::string file_path);
void loadtxt_1(std::vector<unsigned int> &data, const std::string file_path);

// Find local maxima of vector
// Returns vector of indices
std::vector<unsigned int> find_peaks(std::vector<double> &vec);

// Like find_peaks, but negates input and finds peaks within [0,1]
std::vector<unsigned int> find_mid_peaks(std::vector<double> &vec);

// Write 2D array to file
// std::vector<std::vector<double>> data contains vector<double> of length 5 (t x y z ps)
void save_array3_vector(std::vector<array3> &data, const std::string header, const std::string location);
void save_vector_vector(std::vector<std::vector<double>> &data, const std::string header, const std::string location);
void save_vector_unsigned_int(std::vector<std::vector<unsigned int>> &data, const std::string header, const std::string location);
void save_double_vector(std::vector<double> &data, const std::string header, const std::string location);
void save_uint_vector(std::vector<unsigned int> &data, const std::string header, const std::string location);
void save_binary(std::vector<std::vector<unsigned int>> &map, const std::string header, const std::string location);

// Signum function
// From https://stackoverflow.com/questions/1903954/is-there-a-standard-sign-function-signum-sgn-in-c-c
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

// Evaluate a polynomial with ax^n+...+b where poly_coords=<a,...,b>
double eval_poly(std::vector<double> &poly_coords, double x);

// Find first index of array greater than val.
// array must be sorted
template <typename T> unsigned int bisect(std::vector<T> array, T val) {
	// Initialize min and max indices
	unsigned int min = 0;
	unsigned int max = array.size();
	unsigned int mid;
	
	// Perform binary search
	while (min < max) {
		// Midpoint index
		mid = (min + max) / 2;
		
		// Index in lower half
		if (array[mid] > val) {
			max = mid;
		} else { // Index is in upper half
			min = mid + 1; // Index guaranteed larger than mid. Needed so that while loop terminates
		}
	}
	
	return min;
}

// Find absolute value of difference of two values
template <typename T> T abs_diff(T val1, T val2) {
	if (val1 < val2) {return val2 - val1;}
	else {return val1 - val2;}
}

#endif