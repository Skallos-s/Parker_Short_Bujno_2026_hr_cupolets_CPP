// File: TEST_MAIN.cpp
// Purpose: Testing file
// Author: Daniel Bujno

#include <cassert>
#include <algorithm>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iostream>
#include <vector>

// #include <bits/stdc++.h>

#include "microcontrol.h"
#include "macrocontrol.h"
#include "cupolet_search.h"
#include "control_planes.h"
#include "hindmarsh_rose.h"
#include "array3.h"
#include "rk4.h"
#include "keep_data.h"
#include "helper.h"
#include "coding_fcn.h"
#include "time_series.h"
#include "cupolet_from_impulse.h"

int main(int argc, char* argv[]) {
	
	// GENERATION FLAGS
	bool GENERATE_TIME_SERIES_DATA  = false;
	bool GENERATE_CONTROL_PLANE     = false;
	bool GENERATE_PS2_CONTROL_PLANE = false;
	bool GENERATE_CODING_FUNCTIONS  = false;
	bool DO_POLYNOMIAL_REGRESSION   = false;
	bool GENERATE_PS2_INITIAL_PTS   = false;
	bool GENERATE_PS2_PREIMAGE      = false;
	bool ESTABLISH_MICRO_MACRO      = false;
	
	// SIMULATION FLAGS
	bool RUN_ONE_CUPOLET  = false;
	bool RUN_ALL_CUPOLETS = true;
	
	// Hindmarsh Rose attractor
	hindmarsh_rose hr = hindmarsh_rose();
	
	// Bin and crossing data
	unsigned int bins = 1600;
	unsigned int crossings = 16;
	
	// Save directories
	std::string direc = "save_direct";
	std::string bin_rn_direc = direc + "/bins_" + std::to_string(bins) + "_rN_" + std::to_string(crossings);
	
	// Time series data initial conditions
	array3 IC = array3(1,1,3);
	
	unsigned int N = 1280000;
	//double t_max = 10000;
	double dt = 1.0/128;
	
	// Generate time series data
	if (GENERATE_TIME_SERIES_DATA) {generate_time_series(hr, IC, direc, dt, N);}
	
	// Generate control plane data
	if (GENERATE_CONTROL_PLANE || GENERATE_PS2_CONTROL_PLANE) {
		find_surfaces(hr, GENERATE_CONTROL_PLANE, GENERATE_PS2_CONTROL_PLANE);
	}
	
	// Generate coding functions
	if (GENERATE_CODING_FUNCTIONS){
		coding_fcn(hr, direc, bin_rn_direc, dt, crossings, bins, DO_POLYNOMIAL_REGRESSION);
	}
	
	// Generate  PS2 initial points data
	if (GENERATE_PS2_INITIAL_PTS) {
		establish_split_control_plane_bins(hr, direc, bin_rn_direc, dt, bins);
	}
	
	// Generate preimage of PS1 inits onto PS2a
	if (GENERATE_PS2_PREIMAGE) {
		establish_ps1_preimage(hr, direc, bin_rn_direc, dt, bins);
	}
	
	// Establish micro and macro controls
	if (ESTABLISH_MICRO_MACRO) {
		establish_microcontrol(hr, direc, bin_rn_direc, dt, bins);
		establish_macrocontrol(bin_rn_direc);
	}

	
	
	
	if (RUN_ONE_CUPOLET) {
		// Single cupolet
		std::vector<unsigned int> ctrl{0,0,1};
		
		// Starting bin locations for cupolets
		std::vector<unsigned int> bin_start = find_cupolets_start_bins(ctrl, bin_rn_direc);
		
		std::cout << "Number of Cupolets Found: " << bin_start.size() << std::endl;
		
		// Cupolet choice
		unsigned int CUPOLET_INDEX = 0;
		
		// Bound index
		CUPOLET_INDEX = std::min(CUPOLET_INDEX, static_cast<unsigned int>(bin_start.size()));
		
		// How many periods to record
		unsigned int PERIOD_COUNT = 50;
		
		create_cupolet(hr, dt, bins, bin_start[CUPOLET_INDEX], PERIOD_COUNT, ctrl, direc, bin_rn_direc);
		create_time_series_from_impulse(hr, dt, bins, bin_start[CUPOLET_INDEX], ctrl, direc, bin_rn_direc);
	}
	
	if (RUN_ALL_CUPOLETS) {
		// Create save directory for all cupolets
		check_direc(bin_rn_direc + "/primitive_cupolets");
		
		// Loop over control sequence lengths
		for (int length = 1; length <= 6; length++) {
			// Loop over control sequences of set length
			// Encoded as binary
			for (int ctrl_binary = 0; ctrl_binary < 1 << length; ctrl_binary++) {
				// Flag for non-primitivity detection
				bool skip_ctrl = false;
				
				// Loop over all non-trivial cyclic permutations to check primitivity of control sequence
				for (int cyc_perm = 1; cyc_perm < length; cyc_perm++) {
					int val1 = ctrl_binary >> cyc_perm;               // left-most (length-cyc_perm) bits of ctrl_binary
					int val2 = ctrl_binary & ((1 << cyc_perm) - 1);     // right-most cyc_perm bits of ctrl_binary
					int val3 = val1 + (val2 << (length - cyc_perm)); // cyclic permutation of ctrl_binary
					
					// If a non-trivial cyclic permutation is equal to the original ctrl sequence, it is not primitive
					// If not primitive, break out of loop
					skip_ctrl = val3 == ctrl_binary;
					if (skip_ctrl) {break;}
				}
				
				// If not primitive, skip loop iteration
				if (skip_ctrl) {continue;}
				
				// Cupolet name string and array
				std::string cupolet_name = "C";
				std::vector<unsigned int> ctrl;
				for (unsigned int i = length; i > 0; i--) {
					cupolet_name += (ctrl_binary >> (i - 1) & 1) ? "1" : "0";
					ctrl.push_back(ctrl_binary >> (i - 1) & 1);
				}
				
				std::cout << "Finding cupolets " << cupolet_name << std::endl;
				
				// Control sequence folder file path
				std::string ctrl_direc = bin_rn_direc + "/primitive_cupolets/" + cupolet_name;
				
				// Create save directory for given control sequence
				check_direc(ctrl_direc);
				
				// Starting bin locations for cupolets
				std::vector<unsigned int> bin_start = find_cupolets_start_bins(ctrl, bin_rn_direc);
				
				// Loop over bin starts
				for (unsigned int bin_index = 0; bin_index < bin_start.size(); bin_index++) {
					// Cupolet folder file path
					std::string cupolet_direc = ctrl_direc + "/bin_" + std::to_string(bin_start[bin_index]);
					
					std::cout << "  Bin " << bin_start[bin_index] << std::endl;
					
					// Create save directory for given cupolet
					check_direc(cupolet_direc);
					check_direc(cupolet_direc + "/text_data");
					check_direc(cupolet_direc + "/image_data");
					
					create_cupolet(hr, dt, bins, bin_start[bin_index], 1, ctrl, direc, bin_rn_direc, true);
					create_cupolet(hr, dt, bins, bin_start[bin_index], 50, ctrl, direc, bin_rn_direc, true);
				}
			}
		}
		
	}
	
	
	//create_linked_cupolet(hr, dt, bins, bin_start[k], bin_start[k]+20, 0.0001, ctrl, direc, bin_rn_direc);
	
	//save_time_series(hr, dt, bins, bin_start[k], ctrl, direc, bin_rn_direc);
	
	//for (unsigned int i = 0; i < bins; i++) {
	//	save_impulse_function(hr, dt, bins, i, 1, direc, bin_rn_direc);
	//}
	
	//std::vector<std::vector<std::vector<unsigned int>>> var = find_cupolets(ctrl, bin_rn_direc);
	//std::cout << var.size() << " " << var[k].size() << std::endl;
	//
	//for (unsigned int i = 0; i < var[k].size(); i++) {
	//	std::cout << var[k][i][0] << " " << var[k][i][1] << " " << var[k][i][2] << std::endl;
	//}
	
	//std::vector<std::vector<double>> cts = cupolet_time_series(hr, 1.0/128, bins, var[k], direc, bin_rn_direc, true, 0.5);
	//save_data(cts, "Time series data of cupolet C001", bin_rn_direc + "/cupolet_time_series_data.txt");
	
	
	
	return 0;
}
