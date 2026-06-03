// File: microcontrol.cpp
// Purpose: Microcontrol implementation file
// Author: Daniel Bujno


#include "hindmarsh_rose.h"
#include "control_planes.h"
#include "microcontrol.h"
#include "keep_data.h"
#include "array3.h"
#include "helper.h"
#include "rk4.h"

#include <string>
#include <vector>
#include <iostream>

// Reads in neuron object, direc directory, number of bins bins, and time step dt.
// Generates ps0_micro and ps1_micro that represents mapping for each bin and save to direc.
void establish_microcontrol(hindmarsh_rose &neuron, const std::string direc, const std::string bin_rn_direc, double dt, unsigned int bins) {
	
	std::cout << "Creating microcontrol sequence based on cupolet control..." << std::endl;
	
	// Directory where to store microcontrol info
	const std::string store_direc = bin_rn_direc + "/microcontrol";
	
	// Create the directory if it does not exist
	check_direc(store_direc);
	
	// Create the microcontrol map for each control plane
	std::vector<std::vector<unsigned int>> ps0_micro = micro_control_map(neuron, direc, bin_rn_direc, dt, bins, 0);
	std::vector<std::vector<unsigned int>> ps1_micro = micro_control_map(neuron, direc, bin_rn_direc, dt, bins, 1);
	
	save_vector_unsigned_int(ps0_micro, "Microcontrol map for PS0", store_direc + "/ps0_microcontrol.txt");
	save_vector_unsigned_int(ps1_micro, "Microcontrol map for PS1", store_direc + "/ps1_microcontrol.txt");
	
	std::cout << "Saved microcontrol map to " + store_direc << std::endl;
}

// Reads in neuron object, dt time step, number of bins bins, ps control plane, and direc directory. Returns corresponding PS control plane microcontrol mapping.
std::vector<std::vector<unsigned int>> micro_control_map(hindmarsh_rose &neuron, const std::string direc, const std::string bin_rn_direc, double dt, unsigned int bins, unsigned int ps) {
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
	
	// Create empty array to store microcontrol map
	std::vector<std::vector<unsigned int>> psmicro_map;
	
	// Temporary initialization
	std::vector<array3> verts0;
	std::vector<array3> verts1;
		
	// Fill temp vector with values
	verts0.push_back(array3(ps0x[0], ps0y[0], 0));
	verts0.push_back(array3(ps0x[1], ps0y[1], 0));
	verts1.push_back(array3(ps1x[0], ps1y[0], 0));
	verts1.push_back(array3(ps1x[1], ps1y[1], 0));
	
	// Iterate through each bin and find the microcontrol map
	for (unsigned int i = 0; i < bins; i++) {
		// Initial point dependent on which plane to start
		array3 curr = ps == 0 ? ps0inits[i] : ps1inits[i];
		
		// Start with t = 0
		double t = 0;
		
		// Contiue to integrate system until control plane is crossed
		bool crossed_bool = false;
		
		// Push empty vector row to psmicro_map
		std::vector<unsigned int> row{0,0};
		psmicro_map.push_back(row);
		
		while (not crossed_bool) {
			// Integrate one step forward
			array3 next = rk4(curr, dt, &hindmarsh_rose::hr_dynamics, neuron);
			
			// Check if PS0 has been crossed
			if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 0, verts0)) {
				// Find where trajectory intersects control plane
				std::vector<double> tyzxp = rk4_henon(array3(t, curr.get(1), curr.get(2)), curr.get(0), -(curr.get(0)-ps0x[0]), &hindmarsh_rose::hr_dx_dynamics, neuron);
                
				// Find the left side of the intersecting bin endpoint
				unsigned int ii = bisect(ps0endpts, tyzxp[1]) - 1;
				
				// Store the values
				psmicro_map[i][0] = ii;
				psmicro_map[i][1] = 0;
				
				// Mark as crossed to move to next bin
				crossed_bool = true;
				
				// Check if PS1 has been crossed
			} else if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 1, verts1)) {
				// Find where trajectory intersects control plane
				std::vector<double> txzyp = rk4_henon(array3(t, curr.get(0), curr.get(2)), curr.get(1), -(curr.get(1)-ps1y[0]), &hindmarsh_rose::hr_dy_dynamics, neuron);
                
				// Find the left side of the intersecting bin endpoint
				unsigned int ii = bisect(ps1endpts, txzyp[1]) - 1;
				
				// Store the values
				psmicro_map[i][0] = ii;
				psmicro_map[i][1] = 1;
				
				// Mark as crossed to move to next bin
				crossed_bool = true;
			}
			
			// Set new values to old values
			curr = next;
			
			// Move time forward
			t += dt;
		}
	}
	
	// Return the microcontrol map
	return psmicro_map;
}
