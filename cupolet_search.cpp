// File: cupolet_search.cpp
// Purpose: Cupolet search methods implementation file
// Author: Daniel Bujno

#include "cupolet_search.h"
#include "control_planes.h"
#include "hindmarsh_rose.h"
#include "helper.h"
#include "array3.h"
#include "rk4.h"

#include <string>
#include <vector>

#include <iostream>

// Finds all cupolets with given ctrl sequence of length less than limit, that star ton PS1 plane.
// Returns vector of cupolets. Each cupolet is a vector of pairs {index, ps}.
std::vector<std::vector<std::vector<unsigned int>>> find_cupolets(std::vector<unsigned int> &ctrl, std::string bin_rn_direc, unsigned int limit) {
	// Read in PS0/PS1 micro/macro control maps.
	// ps0ps and ps1ps are ps columns in ps0micro and ps1micro respectively.
	std::vector<unsigned int> ps0micro, ps1micro, ps0macro, ps1macro, ps0ps, ps1ps;
	
	loadtxt_2(ps0micro, ps0ps, bin_rn_direc + "/microcontrol/ps0_microcontrol.txt");
	loadtxt_2(ps1micro, ps1ps, bin_rn_direc + "/microcontrol/ps1_microcontrol.txt");
	loadtxt_1(ps0macro,        bin_rn_direc + "/macrocontrol/ps0_macrocontrol.txt");
	loadtxt_1(ps1macro,        bin_rn_direc + "/macrocontrol/ps1_macrocontrol.txt");
	
	std::vector<std::vector<std::vector<unsigned int>>> found_cupolets;

	// Iterate through all possible initial conditions on PS1 to find a cupolet
	for (unsigned int i = 0; i < ps0micro.size(); i++) {
		// Boolean that determines if cupolet has been found
		bool cupolet = false;
		
		// Index of control
		unsigned int ci = 0;
		
		// Initial starting position
		std::vector<std::vector<unsigned int>> cbins;
		std::vector<unsigned int> entry{i,1,2}; // Also last element in cbins
		cbins.push_back(entry);
		
		// Find copulets until limit is reached
		unsigned int checks = 0;
		
		// While cupolet is not found, map the bin connections
		while ((!cupolet) and (checks < limit)) {
			// If the current control is a 0, apply a microcontrol
			if (ctrl[ci] == 0) {
				// No push applied
				entry[2] = 0;
				
				// If the current bin is on PS1 apply PS1 micrcontrol
				if (entry[1] == 1) {
					// next entry
					entry[1] = ps1ps[entry[0]];
					entry[0] = ps1micro[entry[0]];
					
					// save entry
					cbins.push_back(entry);
					
					// PS1 plane intersected
					checks++;
					
					// If the current bin is on PS0 apply PS0 micrcontrol
				} else {
					// next entry
					entry[1] = ps0ps[entry[0]];
					entry[0] = ps0micro[entry[0]];
					
					// save entry
					cbins.push_back(entry);
				}
				
				// If the curret control is a 1, apply a macrocontrol
			} else {
				// Push applied
				entry[2] = 1;
				
				// If the current bin is on PS1 apply PS1 micrcontrol
				if (entry[1] == 1) {
					// perturbed entry
					entry[0] = ps1macro[entry[0]];
					
					// save entry
					cbins.push_back(entry);
					
					// next entry
					entry[2] = 3;
					entry[1] = ps1ps[entry[0]];
					entry[0] = ps1micro[entry[0]];
					
					// save entry
					cbins.push_back(entry);
					
					// PS1 plane intersected
					checks++;
					
					// If the current bin is on PS0 apply PS0 micrcontrol
				} else {
					// perturbed entry
					entry[0] = ps0macro[entry[0]];
					
					// save entry
					cbins.push_back(entry);
					
					// next entry
					entry[2] = 3;
					entry[1] = ps0ps[entry[0]];
					entry[0] = ps0micro[entry[0]];
					
					// save entry
					cbins.push_back(entry);
				}
			}
			
			// Tick control index
			ci++;
			if (ci == ctrl.size()) {ci = 0;}
			
			// Check if cupolet has been found
			// First entry must equal last entry
			if ((cbins[0][0] == entry[0]) and (cbins[0][1] == entry[1]) and (ci == 0)) {
				cupolet = true;
				found_cupolets.push_back(cbins);
			}
		}
	}
	
	// return found cupolets
	return found_cupolets;
}

// Generates cupolet based off of control. Reads in neuron state, dt, number of bins and
// crossings, control sequence, number of iterations, directory to pull info, and ctrl0,ctrl1 arrays
// that tell how to implement macrocontrol. Returns time series of cupolet.
std::vector<std::vector<double>> cupolet_time_series(hindmarsh_rose &neuron, double dt, unsigned int bins, std::vector<std::vector<unsigned int>> cupolet, const std::string direc, const std::string bin_rn_direc) {
	// Read in the control plane initial conditions
	std::vector<array3> ps0inits, ps1inits;
	
	loadtxt_1(ps0inits, bin_rn_direc + "/coding_fcn/ps0_bin_inits.txt");
	loadtxt_1(ps1inits, bin_rn_direc + "/coding_fcn/ps1_bin_inits.txt");
	
	// Read in the control plane endpoints for the polynomial fits
	std::vector<double> ps0endpts, ps1endpts;
	
	loadtxt_1(ps0endpts, bin_rn_direc + "/coding_fcn/ps0_bin_endpoints.txt");
	loadtxt_1(ps1endpts, bin_rn_direc + "/coding_fcn/ps1_bin_endpoints.txt");
	
	// Read in the vertices of each control plane
	std::vector<double> ps0x, ps0y, ps0z;
	std::vector<double> ps1x, ps1y, ps1z;
	
	loadtxt_3(ps0x, ps0y, ps0z, direc + "/control_planes/ps0_vertices.txt");
	loadtxt_3(ps1x, ps1y, ps1z, direc + "/control_planes/ps1_vertices.txt");
	
	// Temporary initialization
	std::vector<array3> verts0;
	std::vector<array3> verts1;
		
	// Fill temp vector with values
	verts0.push_back(array3(ps0x[0], ps0y[0], 0));
	verts0.push_back(array3(ps0x[1], ps0y[1], 0));
	verts1.push_back(array3(ps1x[0], ps1y[0], 0));
	verts1.push_back(array3(ps1x[1], ps1y[1], 0));
	
	// Initialize time series array (t,x,y,z)
	std::vector<std::vector<double>> time_series;
	
	// cupolet array index
	unsigned int i = 0;
	
	// Start with t = 0
	double t = 0;
	
	while (i < cupolet.size() - 1) {
		// If push occurs, skip next entry
		if (cupolet[i + 1][2] == 1) {i++;}
		
		// Bin and plane
		unsigned int bn = cupolet[i][0];
		unsigned int ps = cupolet[i][1];
		
		// Initial point dependent on which plane to start
		array3 curr = ps == 0 ? ps0inits[bn] : ps1inits[bn];
		
		// Contiue to integrate system until control plane is crossed
		bool crossed_bool = false;
		
		// Store starting position into time series
		std::vector<double> row{t, curr.get(0), curr.get(1), curr.get(2)};
		time_series.push_back(row);
		
		while (not crossed_bool) {
			// Integrate one step forward
			array3 next = rk4(curr, dt, &hindmarsh_rose::hr_dynamics, neuron);
			
			// Check if PS0 has been crossed
			if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 0, verts0)) {
				// Find where trajectory intersects control plane
				std::vector<double> tyzxp = rk4_henon(array3(t, curr.get(1), curr.get(2)), curr.get(0), -(curr.get(0)-ps0x[0]), &hindmarsh_rose::hr_dx_dynamics, neuron);
				
				// Mark as crossed to move to next bin
				crossed_bool = true;
				
				// Move time forward
				t = tyzxp[0];
				
				// Store ending position into time series
				std::vector<double> row{t, next.get(0), next.get(1), next.get(2)};
				time_series.push_back(row);
				
				// Check if PS1 has been crossed
			} else if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 1, verts1)) {
				// Find where trajectory intersects control plane
				std::vector<double> txzyp = rk4_henon(array3(t, curr.get(0), curr.get(2)), curr.get(1), -(curr.get(1)-ps1y[0]), &hindmarsh_rose::hr_dy_dynamics, neuron);
				
				// Mark as crossed to move to next bin
				crossed_bool = true;
				
				// Move time forward
				t = txzyp[0];
				
				// Store ending position into time series
				std::vector<double> row{t, next.get(0), next.get(1), next.get(2)};
				time_series.push_back(row);
			} else {
				// Move time forward
				t += dt;
				
				// Store position into time series
				std::vector<double> row{t, next.get(0), next.get(1), next.get(2)};
				time_series.push_back(row);
			}
			
			// Set new values to old values
			curr = next;
		}
		
		// Increment i
		i++;
	}
	
	return time_series;
}