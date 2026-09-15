// File: split_control.cpp
// Purpose: Split control planes transfer implementation file
// Author: Daniel Bujno

#include "split_control.h"
#include "control_planes.h"
#include "hindmarsh_rose.h"
#include "helper.h"
#include "array3.h"
#include "rk4.h"

#include <string>
#include <vector>
#include <iostream>

// Reads in neuron object, direc directory, number of bins bins, and time step dt.
// Generates the initial points of the PS2b bins.
void establish_split_control_plane_bins(hindmarsh_rose &neuron, const std::string direc, const std::string bin_rn_direc, double dt, unsigned int bins) {
	
	std::cout << "Creating PS2 bin data..." << std::endl;
	
	// Directory where to store PS2 bin data
	const std::string store_direc = bin_rn_direc + "/coding_fcn";
	
	// Create the directory if it does not exist
	check_direc(store_direc);
	
	// Read in the PS1 control plane initial conditions
	std::vector<array3> ps1inits;
	loadtxt_1(ps1inits, bin_rn_direc + "/coding_fcn/ps1_bin_inits.txt");
	
	// Read in the PS1 control plane endpoints
	std::vector<double> ps1endpts;
	loadtxt_1(ps1endpts, bin_rn_direc + "/coding_fcn/ps1_bin_endpoints.txt");
	
	// Read in the vertices of the PS1 and PS2 control planes
	std::vector<double> ps1x, ps1y, ps1z;
	std::vector<double> ps2x, ps2y, ps2z;
	
	loadtxt_3(ps1x, ps1y, ps1z, direc + "/control_planes/ps1_vertices.txt");
	loadtxt_3(ps2x, ps2y, ps2z, direc + "/control_planes/ps2_vertices.txt");
	
	// Populate p2ax and co
	//p2ax.push_back(ps2x[0]); p2ax.push_back(ps2x[1]);
	//p2ay.push_back(ps2y[0]); p2ay.push_back(ps2y[1]);
	//p2az.push_back(ps2z[0]); p2az.push_back(ps2z[1]);
	//p2bx.push_back(ps2x[2]); p2bx.push_back(ps2x[3]);
	//p2by.push_back(ps2y[2]); p2by.push_back(ps2y[3]);
	//p2bz.push_back(ps2z[2]); p2bz.push_back(ps2z[3]);
	
	// Fill temp vector with values
	verts1.push_back(array3(ps1x[0], ps1y[0], 0));
	verts1.push_back(array3(ps1x[1], ps1y[1], 0));
	
	verts2a.push_back(array3(ps2x[0], ps2y[0], 0));
	verts2a.push_back(array3(ps2x[1], ps2y[1], 0));
	verts2b.push_back(array3(ps2x[3], ps2y[3], 0));
	verts2b.push_back(array3(ps2x[4], ps2y[4], 0));
	
	
	// PS2b bin initial points
	std::vector<array3> ps2b_inits;
	
	// Step forward ps1inits to find ps2b_inits
	for (unsigned int i = 0; i < bins; i++) {
		// Initial point 
		array3 curr = ps1inits[i];
		
		// Contiue to integrate system until control plane is crossed
		bool crossed_bool = false;
		
		while (not crossed_bool) {
			// Integrate one step forward
			array3 next = rk4(curr, dt, &hindmarsh_rose::hr_dynamics, neuron);
			
			// Check if PS1 has been crossed
			if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 1, verts1)) {
				// Find where trajectory intersects control plane
				std::vector<double> txzyp = rk4_henon(array3(t, curr.get(0), curr.get(2)), curr.get(1), -(curr.get(1)-ps1y[0]), &hindmarsh_rose::hr_dy_dynamics, neuron);
				
				// Mark as crossed to move to next bin
				crossed_bool = true;
				
				// Store ending position into ps2binits
				ps2b_inits.push_back(array3(txzyp[1], txzyp[3], txzyp[2]));
			}
			
		}
	}
	
	return;
}
