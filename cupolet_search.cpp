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


// Finds all cupolets with given ctrl sequence of length less than limit, that star ton PS1 plane.
// Returns vector of cupolets. Each cupolet is a vector of pairs {index, ps}.
std::vector<std::vector<std::vector<unsigned int>>> find_cupolets(std::vector<unsigned int> &ctrl, const std::string bin_rn_direc, unsigned int limit) {
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
		std::vector<unsigned int> entry{i,1,ctrl[0]}; // Also last element in cbins
		cbins.push_back(entry);
		
		// entry[0] is bin #
		// entry[1] is PS  #
		// entry[2] is ctrl value of previous entry
		
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
			
			// When all starting positions are obtained and cupolets cyclically sorted
			// to start at the earliest bin at the start of the control sequence,
			// a cupolet will never go to a lower bin at the start of the sequence
			// than the starting (after sorting) value.
			// Hence, if a cupolet reaches a lower bin than starting at the start
			// of the control sequence, it has been already found and can be discarded.
			if ((cbins[0][0] < entry[0]) and (cbins[0][1] == entry[1]) and (ci == 0)) {
				cupolet = true;
				// Skip cupolet
			}
			
			
			// Check if cupolet has been found
			// First entry must equal last entry
			 else if ((cbins[0][0] == entry[0]) and (cbins[0][1] == entry[1]) and (ci == 0)) {
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
std::vector<std::vector<double>> cupolet_time_series(hindmarsh_rose &neuron, double dt, unsigned int bins, std::vector<std::vector<unsigned int>> cupolet, const std::string direc, const std::string bin_rn_direc, bool split_planes, double BEZIER_STRENGTH) {
	
	// Read in the control plane initial conditions
	std::vector<array3> ps0inits, ps1inits, ps2inits;
	
	loadtxt_1(ps0inits, bin_rn_direc + "/coding_fcn/ps0_bin_inits.txt");
	loadtxt_1(ps1inits, bin_rn_direc + "/coding_fcn/ps1_bin_inits.txt");
	loadtxt_1(ps2inits, bin_rn_direc + "/coding_fcn/ps2_bin_inits.txt");
	
	// Read in the travel time between PS1 and PS2b
	std::vector<double> ps2b_time;
	
	loadtxt_1(ps2b_time, bin_rn_direc + "/coding_fcn/ps2_travel_time.txt");
	
	// Read in the previous positions of initial conditions of PS2b
	std::vector<array3> ps2b_prev;
	
	loadtxt_1(ps2b_prev, bin_rn_direc + "/coding_fcn/ps2_prev_position.txt");
	
	// Read in the vertices of each control plane
	std::vector<double> ps0x, ps0y, ps0z;
	std::vector<double> ps1x, ps1y, ps1z;
	std::vector<double> ps2x, ps2y, ps2z;
	
	loadtxt_3(ps0x, ps0y, ps0z, direc + "/control_planes/ps0_vertices.txt");
	loadtxt_3(ps1x, ps1y, ps1z, direc + "/control_planes/ps1_vertices.txt");
	loadtxt_3(ps2x, ps2y, ps2z, direc + "/control_planes/ps2_vertices.txt");
	
	// Plane data used for checking crossing
	std::vector<array3> verts0;
	std::vector<array3> verts1;
	std::vector<array3> verts2a;
	std::vector<array3> verts2b;
	
	// Fill plane data with values
	verts0.push_back(array3(ps0x[0], ps0y[0], 0));
	verts0.push_back(array3(ps0x[1], ps0y[1], 0));
	verts1.push_back(array3(ps1x[0], ps1y[0], 0));
	verts1.push_back(array3(ps1x[1], ps1y[1], 0));
	
	verts2a.push_back(array3(ps2x[0], ps2y[0], 0));
	verts2a.push_back(array3(ps2x[1], ps2y[1], 0));
	verts2b.push_back(array3(ps2x[3], ps2y[3], 0));
	verts2b.push_back(array3(ps2x[4], ps2y[4], 0));
	
	// Initialize time series array (t,x,y,z)
	std::vector<std::vector<double>> time_series;
	
	// cupolet array index
	unsigned int i = 0;
	
	// Start with t = 0
	double t = 0;
	
	while (i + 1 < cupolet.size()) {
		// If push occurs, skip next entry
		if (cupolet[i + 1][2] == 1) {i++;}
		
		// Bin and plane
		unsigned int bn = cupolet[i][0];
		unsigned int ps = cupolet[i][1];
		
		// Increment travel time if PS1 has been passed and planes are split
		// First entry (i=0) could not have passed PS1
		//if (split_planes && i != 0 && ps != 0) {t += ps2b_time[bn];}
		
		// Initial point dependent on which plane to start
		array3 curr = array3(0,0,0);
		if (ps == 0) {curr = ps0inits[bn];}
		else {
			if (split_planes) {
				curr = ps2inits[bn];
			} else {
				curr = ps1inits[bn];
			}
		}
		
		// Contiue to integrate system until control plane is crossed
		bool crossed_bool = false;
		
		// Don't save output when PS2a is crossed
		bool PS2a_crossed_bool = false;
		
		// Store starting position into time series
		// Skip if point is on PS2b and is not the first point
		std::vector<double> row{t, curr.get(0), curr.get(1), curr.get(2)};
		if (ps == 0 || !split_planes || i == 0) {time_series.push_back(row);}
		
		// Used for split planes bezier curve
		double t0, t1;
		array3 P0, Q1;
		
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
				std::vector<double> row{t, tyzxp[3], tyzxp[1], tyzxp[2]};
				time_series.push_back(row);
				
				// Check if PS1 has been crossed
			} else if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 1, verts1)) {
				// Find where trajectory intersects control plane
				std::vector<double> txzyp = rk4_henon(array3(t, curr.get(0), curr.get(2)), curr.get(1), -(curr.get(1)-ps1y[0]), &hindmarsh_rose::hr_dy_dynamics, neuron);
				
				// Mark as crossed to move to next bin
				crossed_bool = true;
				
				// Move time forward
				t = txzyp[0];
				
				// Store ending position into time series if planes are not split
				if (!split_planes) {
					std::vector<double> row{t, txzyp[1], txzyp[3], txzyp[2]};
					time_series.push_back(row);
					
					// Fill space between PS2a and PS2b
				} else {
					// Next index
					unsigned int j = i + 1;
					
					// If push occurs, skip next entry
					if (j + 1 < cupolet.size() && cupolet[j + 1][2] == 1) {j++;}
		
					// Next bin
					unsigned int bn_j = cupolet[j][0];
					
					// End time step for bezier curve
					t += ps2b_time[bn_j];
					t1 = t;
					
					// Do bezier splining
					array3 Q2 = ps2b_prev[bn_j];
					array3 P3 = ps2inits[bn_j];
					
					// Scaling factor for bezier curve
					double SCALE1 = BEZIER_STRENGTH / (P0.get(1) - Q1.get(1));
					double SCALE2 = BEZIER_STRENGTH / (Q2.get(1) - P3.get(1));
					
					// Finish doing bezier splining
					array3 P1 = 3 * (P0 + SCALE1 * (Q1 - P0));
					array3 P2 = 3 * (P3 + SCALE2 * (Q2 - P3));
					
					// Number of time steps between P0 and P3
					unsigned int STEP_COUNT = 20;
					
					// Step size for bezier curve
					double dtbc1 = (t1 - t0) / STEP_COUNT;
					
					// Bezier curve
					for (unsigned int k = 0; k <= STEP_COUNT; k++) {
						// Time step along unit interval
						double dtui = 1.0 / STEP_COUNT * k;
						
						// t^2 and (1-t)^2
						double tsq = dtui * dtui;
						double rts = (1 - dtui) * (1 - dtui);
						
						// Polynomials (1-t)^3, t(1-t)^2, t^2(1-t), and t^3
						double P0C = rts * (1 - dtui);
						double P1C = rts * dtui;
						double P2C = tsq * (1 - dtui);
						double P3C = tsq * dtui;
						
						// Point on bezier curve
						array3 BC = P0 * P0C + P1 * P1C + P2 * P2C + P3 * P3C;
						
						// Intermediate time step
						double dtbc = t0 + dtbc1 * k;
						
						// Save point
						std::vector<double> row{dtbc, BC.get(0), BC.get(1), BC.get(2)};
						time_series.push_back(row);
					}
				}
				
				// Check if PS2a has been crossed if control planes are split
				// Short-circuiting used 
			} else if (split_planes && crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 1, verts2a)) {
				// Find where trajectory intersects control plane
				std::vector<double> txzyp = rk4_henon(array3(t, curr.get(0), curr.get(2)), curr.get(1), -(curr.get(1)-ps2y[0]), &hindmarsh_rose::hr_dy_dynamics, neuron);
				
				// Mark PS2a as crossed to ignore saving of unperturbed path
				PS2a_crossed_bool = true;
				
				// Move time forward
				t = txzyp[0];
				
				// Initial time step for bezier curve
				t0 = t;
				
				// Control points for bezier curve
				P0 = array3(txzyp[1], txzyp[3], txzyp[2]);
				Q1 = next;
				
				// Store ending position into time series
				std::vector<double> row{t, txzyp[1], txzyp[3], txzyp[2]};
				time_series.push_back(row);
			} else {
				// Move time forward
				t += dt;
				
				// Store position into time series if PS2a has not been crossed
				if (!PS2a_crossed_bool) {
					std::vector<double> row{t, next.get(0), next.get(1), next.get(2)};
					time_series.push_back(row);
				}
			}
			
			// Set new values to old values
			curr = next;
		}
		
		// Increment i
		i++;
	}
	
	return time_series;
}