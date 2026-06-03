// File: macrocontrol.cpp
// Purpose: macrocontrol implementation file
// Author: Daniel Bujno


#include "macrocontrol.h"
#include "keep_data.h"
#include "helper.h"

#include <cmath>
#include <cfloat>
#include <string>
#include <vector>
#include <climits>
#include <iostream>

// Reads in the top level directory name. Generates the macrocontrol sequence that corresponds
// to the number of bins and crossings and saves in directory /direc/macrocontrol.
void establish_macrocontrol(const std::string bin_rn_direc) {
	
	std::cout << "Creating macrocontrol sequence based on cupolet control..." << std::endl;
	
	// Directory where to store macrocontrol info
	const std::string store_direc = bin_rn_direc + "/macrocontrol";
	
	// Create the directory if it does not exist
	check_direc(store_direc);
	
	// Read in the coding function for each control plane
	std::vector<double> rn0, rn1;
	
	loadtxt_1(rn0, bin_rn_direc + "/coding_fcn/ps0_rn_fcn.txt");
	loadtxt_1(rn1, bin_rn_direc + "/coding_fcn/ps1_rn_fcn.txt");
	
	// Find the macrocontrol map for each control plane, and the difference map
	std::vector<unsigned int> ctrl0, ctrl1;
	std::vector<double>       diff0, diff1;
	
	macro_control_map(ctrl0, diff0, rn0);
	macro_control_map(ctrl1, diff1, rn1);
	
	// Save the macrocontrol maps for PS0 and PS1
	save_uint_vector(ctrl0, "PS0 macrocontrol map.", store_direc + "/ps0_macrocontrol.txt");
	save_uint_vector(ctrl1, "PS1 macrocontrol map.", store_direc + "/ps1_macrocontrol.txt");
	
	// Save the macrocontrol difference maps for PS0 and PS1
	save_double_vector(diff0, "PS0 macrocontrol difference map.", store_direc + "/ps0_diff_macrocontrol.txt");
	save_double_vector(diff1, "PS1 macrocontrol difference map.", store_direc + "/ps1_diff_macrocontrol.txt");
	
	std::cout << "Saved macrocontrol map to " + store_direc << std::endl;
}

// Reads in a coding sequence rn. Write to an array ctrl that maps bin B to the nearest
// smallest difference of the coding sequence and an array diff that says the size of that difference.
void macro_control_map(std::vector<unsigned int> &ctrl, std::vector<double> &diff, std::vector<double> rn) {
	// Define the number of bins
	unsigned int bins = rn.size();
	
	// Iterate through each bin and find the macrocontrol destination
	for (unsigned int i = 0; i < bins; i++) {
		// Minimum positive difference between ith bin and any other bin
		double min_diff = DBL_MAX;
		
		// Index of right-most closest bin with difference min_diff
		unsigned int index = 0;
		
		// Helper variable to find closest bin index with difference min_diff
		unsigned int index_offset = UINT_MAX;
		
		// Loop over all bins
		for (unsigned int j = 0; j < bins; j++) {
			// Difference between jth and ith bins
			double ndiff = std::abs(rn[j] - rn[i]);
			
			// Smallest positive difference found
			if ((ndiff > 0) and (ndiff < min_diff)) {
				// Save minimum difference
				min_diff = ndiff;
				
				// Save index of minimum
				index = j;
				
				// Save offset
				index_offset = abs_diff(i, j);
			} else if (ndiff == min_diff) {
				unsigned int d = abs_diff(i, j);
				if (d <= index_offset) {
					index_offset = d;
					index = j;
				}
			}
		}
		
		// Push to ctrl and diff
		ctrl.push_back(index);
		diff.push_back(min_diff);
	}
	
	return;
}