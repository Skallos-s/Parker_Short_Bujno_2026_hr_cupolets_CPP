// File: TEST_MAIN.cpp
// Purpose: Testing file
// Author: Daniel Bujno

#include <cassert>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iostream>
#include <vector>

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

int main(int argc, char* argv[]) {
	
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
	generate_time_series(hr, IC, direc, dt, N);
	
	// Generate control plane data
	//find_surfaces(hr);
	
	// Generate coding functions
	//coding_fcn(hr, direc, bin_rn_direc, dt, crossings, bins, false);
	
	// Establish micro and macro controls
	//establish_microcontrol(hr, direc, bin_rn_direc, 1.0/128, bins);
	//establish_macrocontrol(bin_rn_direc);
	
	// Single cupolet
	std::vector<unsigned int> ctrl{0,1};
	unsigned int k = 0;
	
	std::vector<std::vector<std::vector<unsigned int>>> var = find_cupolets(ctrl, bin_rn_direc);
	std::cout << var.size() << " " << var[k].size() << std::endl;
	
	for (unsigned int i = 0; i < var[k].size(); i++) {
		std::cout << var[k][i][0] << " " << var[k][i][1] << " " << var[k][i][2] << std::endl;
	}
	
	std::vector<std::vector<double>> cts = cupolet_time_series(hr, 1.0/128, bins, var[k], direc, bin_rn_direc);
	save_vector_vector(cts, "Time series data of cupolet 001", bin_rn_direc + "/cupolet_time_series_data.txt");
	
	
	return 0;
}
