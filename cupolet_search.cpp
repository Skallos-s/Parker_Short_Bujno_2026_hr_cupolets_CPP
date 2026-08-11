// File: cupolet_search.cpp
// Purpose: Cupolet search methods implementation file
// Author: Daniel Bujno

#include "cupolet_search.h"
#include "control_planes.h"
#include "hindmarsh_rose.h"
#include "keep_data.h"
#include "helper.h"
#include "array3.h"
#include "rk4.h"

#include <string>
#include <vector>
#include <iostream>


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
			ci = (ci + 1) % ctrl.size();
			
			// When all starting positions are obtained and cupolets cyclically sorted
			// to start at the earliest bin at the start of the control sequence,
			// a cupolet will never go to a lower bin at the start of the sequence
			// than the starting (after sorting) value.
			// Hence, if a cupolet reaches a lower bin than starting at the start
			// of the control sequence, it has been already found and can be discarded.
			if ((cbins[0][0] > entry[0]) and (cbins[0][1] == entry[1]) and (ci == 0)) {
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
	//std::vector<array3> ps2b_prev;
	//loadtxt_1(ps2b_prev, bin_rn_direc + "/coding_fcn/ps2_prev_position.txt");
	
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
					array3 P3 = ps2inits[bn_j];
					array3 Q2 = neuron.hr_dynamics(P3);
					
					// Scaling factor for bezier curve
					double SCALE1 = std::abs(BEZIER_STRENGTH / Q1.get(1));
					double SCALE2 = std::abs(BEZIER_STRENGTH / Q2.get(1));
					
					// Finish doing bezier splining
					array3 P1 = 3 * (P0 + SCALE1 * Q1);
					array3 P2 = 3 * (P3 - SCALE2 * Q2);
					
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
				Q1 = neuron.hr_dynamics(P0);
				
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

// Finds all cupolets with given ctrl sequence of length less than limit, that star ton PS1 plane.
// Returns start bins of each cupolet
std::vector<unsigned int> find_cupolets_start_bins(std::vector<unsigned int> &ctrl, const std::string bin_rn_direc, unsigned int limit) {
	// Read in PS0/PS1 micro/macro control maps.
	// ps0ps and ps1ps are ps columns in ps0micro and ps1micro respectively.
	std::vector<unsigned int> ps0micro, ps1micro, ps0macro, ps1macro, ps0ps, ps1ps;
	
	loadtxt_2(ps0micro, ps0ps, bin_rn_direc + "/microcontrol/ps0_microcontrol.txt");
	loadtxt_2(ps1micro, ps1ps, bin_rn_direc + "/microcontrol/ps1_microcontrol.txt");
	loadtxt_1(ps0macro,        bin_rn_direc + "/macrocontrol/ps0_macrocontrol.txt");
	loadtxt_1(ps1macro,        bin_rn_direc + "/macrocontrol/ps1_macrocontrol.txt");
	
	std::vector<unsigned int> found_cupolets;

	// Iterate through all possible initial conditions on PS1 to find a cupolet
	for (unsigned int i = 0; i < ps0micro.size(); i++) {
		// Boolean that determines if cupolet has been found
		bool cupolet = false;
		
		// Index of control
		unsigned int ci = 0;
		
		// Initial position
		unsigned int bin_init = i;
		unsigned int ps_init  = 1;
		
		// Current position
		unsigned int bin_curr = i;
		unsigned int ps_curr  = 1;
		
		// Find copulets until limit is reached
		unsigned int checks = 0;
		
		// While cupolet is not found, map the bin connections
		while ((!cupolet) and (checks < limit)) {
			// If the current control is a 0, apply a microcontrol
			if (ctrl[ci] == 0) {
				// If the current bin is on PS1 apply PS1 micrcontrol
				if (ps_curr == 1) {
					// next entry
					ps_curr  = ps1ps[bin_curr];
					bin_curr = ps1micro[bin_curr];
					
					// PS1 plane intersected
					checks++;
					
					// If the current bin is on PS0 apply PS0 micrcontrol
				} else {
					// next entry
					ps_curr  = ps0ps[bin_curr];
					bin_curr = ps0micro[bin_curr];
				}
				
				// If the curret control is a 1, apply a macrocontrol
			} else {
				// If the current bin is on PS1 apply PS1 micrcontrol
				if (ps_curr == 1) {
					// perturbed entry
					bin_curr = ps1macro[bin_curr];
					
					// next entry
					ps_curr  = ps1ps[bin_curr];
					bin_curr = ps1micro[bin_curr];
					
					// PS1 plane intersected
					checks++;
					
					// If the current bin is on PS0 apply PS0 micrcontrol
				} else {
					// perturbed entry
					bin_curr = ps0macro[bin_curr];
					
					// next entry
					ps_curr  = ps0ps[bin_curr];
					bin_curr = ps0micro[bin_curr];
				}
			}
			
			// Tick control index
			ci = (ci + 1) % ctrl.size();
			
			// When all starting positions are obtained and cupolets cyclically sorted
			// to start at the earliest bin at the start of the control sequence,
			// a cupolet will never go to a lower bin at the start of the sequence
			// than the starting (after sorting) value.
			// Hence, if a cupolet reaches a lower bin than starting at the start
			// of the control sequence, it has been already found and can be discarded.
			if ((bin_init > bin_curr) and (ps_init == ps_curr) and (ci == 0)) {
				cupolet = true;
				// Skip cupolet
			}
			
			
			// Check if cupolet has been found
			// First entry must equal last entry
			 else if ((bin_init == bin_curr) and (ps_init == ps_curr)  and (ci == 0)) {
				cupolet = true;
				found_cupolets.push_back(bin_init);
			}
		}
	}
	
	// return found cupolets
	return found_cupolets;
}

// Generates cupolet from control sequence and starting bin. Reads in neuron state, dt, number of bins
// starting bin, control sequence, and save directories. Generates and saves time series of cupolet
void save_time_series(hindmarsh_rose &neuron, double dt, unsigned int bins, unsigned int start_bin, std::vector<unsigned int> &ctrl, const std::string direc, const std::string bin_rn_direc) {
	// Read in the control plane initial conditions
	std::vector<array3> ps0inits, ps1inits, ps2inits;
	
	loadtxt_1(ps0inits, bin_rn_direc + "/coding_fcn/ps0_bin_inits.txt");
	loadtxt_1(ps1inits, bin_rn_direc + "/coding_fcn/ps1_bin_inits.txt");
	loadtxt_1(ps2inits, bin_rn_direc + "/coding_fcn/ps2_bin_inits.txt");
	
	// Read in the vertices of each control plane
	std::vector<double> ps0x, ps0y, ps0z;
	std::vector<double> ps1x, ps1y, ps1z;
	std::vector<double> ps2x, ps2y, ps2z;
	
	loadtxt_3(ps0x, ps0y, ps0z, direc + "/control_planes/ps0_vertices.txt");
	loadtxt_3(ps1x, ps1y, ps1z, direc + "/control_planes/ps1_vertices.txt");
	loadtxt_3(ps2x, ps2y, ps2z, direc + "/control_planes/ps2_vertices.txt");
	
	// Read in the control plane endpoints for the polynomial fits
	std::vector<double> ps0endpts, ps1endpts;
	
	loadtxt_1(ps0endpts, bin_rn_direc + "/coding_fcn/ps0_bin_endpoints.txt");
	loadtxt_1(ps1endpts, bin_rn_direc + "/coding_fcn/ps1_bin_endpoints.txt");
	
	// Read in PS0/PS1 macrocontrol maps.
	std::vector<unsigned int> ps0macro, ps1macro;
	
	loadtxt_1(ps0macro, bin_rn_direc + "/macrocontrol/ps0_macrocontrol.txt");
	loadtxt_1(ps1macro, bin_rn_direc + "/macrocontrol/ps1_macrocontrol.txt");
	
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
	verts2b.push_back(array3(ps2x[2], ps2y[2], 0));
	verts2b.push_back(array3(ps2x[3], ps2y[3], 0));
	
	// Initialize time series arrays (t,x,y,z)
	std::vector<std::vector<double>> true_time_series;
	std::vector<std::vector<double>> gray_time_series;
	
	// Delta between Bezier curve interprolation direction vector and HR derivative
	std::vector<std::vector<double>> HR_bezier_delta;
	
	// Start and end bins of bezier curve
	// Used for file naming when saving HR_bezier_delta;
	unsigned int bezier_start_bin, bezier_end_bin;

	// Index of control
	unsigned int ci = 0;
	
	// Start with t = 0
	double t = 0;
	
	// Integrate until cupolet is fully mapped
	bool CONTINUE_INTEGRATING  = true;
	bool STOP_INTEGRATING_SOON = false;
	
	// Check if PS2a has been crossed
	bool PS2A_CROSSED = false;
	
	// Used for split planes bezier curve
	double t0, t1;
	array3 P0, Q1;
	
	// Initial point dependent on which plane to start
	array3 curr = ps2inits[start_bin];
	
	// Store starting position into time series
	std::vector<double> row{0,curr.get(0), curr.get(1), curr.get(2)};
	true_time_series.push_back(row);
	gray_time_series.push_back(row);
	
	// Indicate which of the PS1 bins have been reached to detect if the cupolet has looped onto itself
	std::vector<bool> PS0_reached, PS1_reached;
	PS0_reached.resize(bins, false);
	PS1_reached.resize(bins, false);
	
	// Indicate that the starting bin on PS1 has been reached
	PS1_reached[start_bin] = true;
	
	// Cupolet name string
	std::string CUPOLET_NAME = "C";
	for (unsigned int i = 0; i < ctrl.size(); i++) {
		CUPOLET_NAME += (ctrl[i]) ? "!" : "0";
	}
	
	// Cupolet generation loop
	while (CONTINUE_INTEGRATING) {
		// Integrate one step forward
		array3 next = rk4(curr, dt, &hindmarsh_rose::hr_dynamics, neuron);
		
		// Check if PS0 has been crossed
		if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 0, verts0)) {
			// Find where trajectory intersects control plane
			std::vector<double> tyzxp = rk4_henon(array3(t, curr.get(1), curr.get(2)), curr.get(0), -(curr.get(0)-ps0x[0]), &hindmarsh_rose::hr_dx_dynamics, neuron);
			
			// Move time forward
			t = tyzxp[0];
			
			// Store position into time series
			row = {t, tyzxp[3], tyzxp[1], tyzxp[2]};
			true_time_series.push_back(row);
			gray_time_series.push_back(row);
			
			// Tick control index
			ci = (ci + 1) % ctrl.size();
			
			// Find the left side of the intersecting bin endpoint
			unsigned int ii = bisect(ps0endpts, tyzxp[1]) - 1;
			
			// Apply macrocontrol
			if (ctrl[ci]) {ii = ps0macro[ii];}
			
			// Center to bin
			next = ps0inits[ii];
			
			// If bin has been reached at the start of control sequence already,
			// set flag to false. Then mark bin as having been reached
			// If flag is false, then integration can stop
			CONTINUE_INTEGRATING = !PS0_reached[ii] || (ci > 0);
			PS0_reached[ii] = ci == 0;
			
			// Store position into time series
			row = {t, next.get(0), next.get(1), next.get(2)};
			true_time_series.push_back(row);
			gray_time_series.push_back(row);
			
		// Check if PS2a has been crossed
		} else if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 1, verts2a)) {
			// Find where trajectory intersects control plane
			std::vector<double> txzyp = rk4_henon(array3(t, curr.get(0), curr.get(2)), curr.get(1), -(curr.get(1)-ps2y[0]), &hindmarsh_rose::hr_dy_dynamics, neuron);
			
			// Move time forward
			t = txzyp[0];
			
			// Initial time step for bezier curve
			t0 = t;
			
			// PS2a has been crossed
			PS2A_CROSSED = true;
			
			// Control points for bezier curve
			P0 = array3(txzyp[1], txzyp[3], txzyp[2]);
			Q1 = neuron.hr_dynamics(P0);
			
			// Snap point to PS2a
			next = P0;
			
			// Store position into gray time series
			row = {t, next.get(0), next.get(1), next.get(2)};
			gray_time_series.push_back(row);
			
			// Reset HR_bezier_delta vector
			HR_bezier_delta = {};
			
		// Check if PS1 has been crossed
		} else if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 1, verts1)) {
			// Find where trajectory intersects control plane
			std::vector<double> txzyp = rk4_henon(array3(t, curr.get(0), curr.get(2)), curr.get(1), -(curr.get(1)-ps1y[0]), &hindmarsh_rose::hr_dy_dynamics, neuron);
			
			// Move time forward
			t = txzyp[0];
			
			// Store position into only gray time series
			row = {t, txzyp[1], txzyp[3], txzyp[2]};
			gray_time_series.push_back(row);
			
			// Tick control index
			ci = (ci + 1) % ctrl.size();
			
			// Find the left side of the intersecting bin endpoint
			unsigned int ii = bisect(ps1endpts, txzyp[1]) - 1;
			
			// Save start bin
			bezier_start_bin = ii;
			
			// Apply macrocontrol
			if (ctrl[ci]) {ii = ps1macro[ii];}
			
			// Save end bin
			bezier_end_bin = ii;
			
			// Center to bin
			next = ps1inits[ii];
			
			// If bin has been reached at the start of control sequence already,
			// set flag to true. Then mark bin as having been reached
			// If flag is true, integration must continue until PS2b
			STOP_INTEGRATING_SOON = PS1_reached[ii] & (ci == 0);
			PS1_reached[ii] = ci == 0;
			
			// Store position into gray time series
			row = {t, next.get(0), next.get(1), next.get(2)};
			gray_time_series.push_back(row);
		
		// Check if PS2b has been crossed
		} else if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 1, verts2b)) {
			// Find where trajectory intersects control plane
			std::vector<double> txzyp = rk4_henon(array3(t, curr.get(0), curr.get(2)), curr.get(1), -(curr.get(1)-ps2y[2]), &hindmarsh_rose::hr_dy_dynamics, neuron);
			
			// Move time forward
			t = txzyp[0];
			
			// Final time step for bezier curve
			t1 = t;
			
			// If bin has been reached in PS1, stop integrating
			CONTINUE_INTEGRATING = !STOP_INTEGRATING_SOON;
			
			// Do bezier splining
			array3 P3 = array3(txzyp[1], txzyp[3], txzyp[2]);
			array3 Q2 = neuron.hr_dynamics(P3);
			
			// Finish doing bezier splining
			array3 P1 = P0 + Q1 / 3 * (t1 - t0);
			array3 P2 = P3 - Q2 / 3 * (t1 - t0);
			
			// Snap point to PS2b
			next = P3;
			
			// PS2b has been crossed
			PS2A_CROSSED = false;
			
			// Store position into gray time series
			row = {t, next.get(0), next.get(1), next.get(2)};
			gray_time_series.push_back(row);
			
			// Number of time steps between P0 and P3
			unsigned int STEP_COUNT = 20;
			
			// Step size for bezier curve
			double dtbc1 = (t1 - t0) / STEP_COUNT;
			
			// Bezier curve
			for (unsigned int k = 0; k <= STEP_COUNT; k++) {
				// Time step along unit interval
				double dtui = 1.0 / STEP_COUNT * k;
				double rvdt = 1.0 - dtui;

				// Polynomials (1-t)^3, t(1-t)^2, t^2(1-t), and t^3
				double P0C = rvdt * rvdt * rvdt;
				double P1C = rvdt * rvdt * dtui;
				double P2C = rvdt * dtui * dtui;
				double P3C = dtui * dtui * dtui;
				
				// Polynomials (1-t)^2, t(1-t), and t^2
				double D0C = rvdt * rvdt;
				double D1C = rvdt * dtui;
				double D2C = dtui * dtui;
				
				// Point on bezier curve
				array3 PBC = P0 * P0C + 3 * P1 * P1C + 3 * P2 * P2C + P3 * P3C;
				
				// Derivative of bezier curve
				// Normalized across unit interval, not time interval
				array3 DBC = 3 * D0C * (P1 - P0) + 6 * D1C * (P2 - P1) + 3 * D2C * (P3 - P2);
				
				// Hindmarsh Rose derivative
				array3 DHR = neuron.hr_dynamics(PBC);
				
				// Delta
				// DBC normalized across time interval
				array3 delta = DBC / (t1 - t0) - DHR;
				
				// Intermediate time step
				double dtbc = t0 + dtbc1 * k;
				
				// Store position into true time series
				row = {dtbc, PBC.get(0), PBC.get(1), PBC.get(2)};
				true_time_series.push_back(row);
				
				// Store delta between HR and Bezier
				row = {dtbc, delta.get(0), delta.get(1), delta.get(2)};
				HR_bezier_delta.push_back(row);
			}
			
			// Save delta time series
			save_data(HR_bezier_delta, "Delta difference between Hindmarsh Rose derivative and Bezier curve derivative",
			          bin_rn_direc + "/delta_" + std::to_string(bezier_start_bin) + "_" + std::to_string(bezier_end_bin) + ".txt");
			
			
		// No plane has been crossed
		} else {
			// Move time forward
			t += dt;
			
			// Store position into time series
			std::vector<double> row{t, next.get(0), next.get(1), next.get(2)};
			if (!PS2A_CROSSED) {true_time_series.push_back(row);}
			gray_time_series.push_back(row);
		}
		
		// Set new values to old values
		curr = next;
	}
	
	save_data(true_time_series, "True time series data of cupolet " + CUPOLET_NAME, bin_rn_direc + "/" + CUPOLET_NAME + "_cupolet_true_time_series_data.txt");
	save_data(gray_time_series, "Gray time series data of cupolet " + CUPOLET_NAME, bin_rn_direc + "/" + CUPOLET_NAME + "_cupolet_gray_time_series_data.txt");
	
	return;
}

// Generate impulse function for difference between true path and Bezier curve from PS2a to PS2b
// Reads in neuron stte, dt, number of bins, starting bin, whether a push is done or not, save directories.
void save_impulse_function(hindmarsh_rose &neuron, double dt, unsigned int bins, unsigned int start_bin, unsigned int ctrl, const std::string direc, const std::string bin_rn_direc) {
	// Read in the control plane initial conditions
	std::vector<array3> ps1inits, ps2ainits, ps2binits;
	
	loadtxt_1(ps1inits , bin_rn_direc + "/coding_fcn/ps1_bin_inits.txt" );
	loadtxt_1(ps2ainits, bin_rn_direc + "/coding_fcn/ps2a_bin_inits.txt");
	loadtxt_1(ps2binits, bin_rn_direc + "/coding_fcn/ps2b_bin_inits.txt");
	
	// Read in the vertices of each control plane
	std::vector<double> ps1x, ps1y, ps1z;
	std::vector<double> ps2x, ps2y, ps2z;
	
	loadtxt_3(ps1x, ps1y, ps1z, direc + "/control_planes/ps1_vertices.txt");
	loadtxt_3(ps2x, ps2y, ps2z, direc + "/control_planes/ps2_vertices.txt");
	
	// Read in travel times from PS2a to PS1 and from PS1 to PS2b
	std::vector<double> travel_time_1, travel_time_2;
	
	loadtxt_1(travel_time_1, bin_rn_direc + "/coding_fcn/ps2a_ps1_travel_times.txt");
	loadtxt_1(travel_time_2, bin_rn_direc + "/coding_fcn/ps1_ps2b_travel_times.txt");
	
	// Read in the control plane endpoints
	std::vector<double> ps1endpts;
	
	loadtxt_1(ps1endpts, bin_rn_direc + "/coding_fcn/ps1_bin_endpoints.txt");
	
	// Read in PS1 macrocontrol maps.
	std::vector<unsigned int> ps1macro;
	
	loadtxt_1(ps1macro, bin_rn_direc + "/macrocontrol/ps1_macrocontrol.txt");
	
	// Delta between Bezier curve interprolation direction vector and HR derivative
	std::vector<std::vector<double>> HR_bezier_delta;
	
	// Ending bin
	unsigned int end_bin = start_bin;
	
	// Apply kick
	if (ctrl) {end_bin = ps1macro[start_bin];}
	
	// Total travel time
	double total_travel_time = travel_time_1[start_bin] + travel_time_2[end_bin];
	
	// Bezier curve control points
	array3 P0 = ps2ainits[start_bin];
	array3 P3 = ps2binits[end_bin];
	
	array3 Q1 = neuron.hr_dynamics(P0);
	array3 Q2 = neuron.hr_dynamics(P3);
	
	array3 P1 = P0 + Q1 / 3 * total_travel_time;
	array3 P2 = P3 - Q2 / 3 * total_travel_time;
	
	// Number of time steps between P0 and P3
	unsigned int STEP_COUNT = 20;
	
	// Step size for bezier curve
	double dtbc1 = total_travel_time / STEP_COUNT;
	
	// Bezier curve
	for (unsigned int k = 0; k <= STEP_COUNT; k++) {
		// Time step along unit interval
		double dtui = 1.0 / STEP_COUNT * k;
		double rvdt = 1.0 - dtui;

		// Polynomials (1-t)^3, t(1-t)^2, t^2(1-t), and t^3
		double P0C = rvdt * rvdt * rvdt;
		double P1C = rvdt * rvdt * dtui;
		double P2C = rvdt * dtui * dtui;
		double P3C = dtui * dtui * dtui;
		
		// Polynomials (1-t)^2, t(1-t), and t^2
		double D0C = rvdt * rvdt;
		double D1C = rvdt * dtui;
		double D2C = dtui * dtui;
		
		// Point on bezier curve
		array3 PBC = P0 * P0C + 3 * P1 * P1C + 3 * P2 * P2C + P3 * P3C;
		
		// Derivative of bezier curve
		// Normalized across unit interval, not time interval
		array3 DBC = 3 * D0C * (P1 - P0) + 6 * D1C * (P2 - P1) + 3 * D2C * (P3 - P2);
		
		// Hindmarsh Rose derivative
		array3 DHR = neuron.hr_dynamics(PBC);
		
		// Delta
		// DBC normalized across time interval
		array3 delta = DBC / total_travel_time - DHR;
		
		// Intermediate time step
		double dtbc = dtbc1 * k;
		
		// Store delta between HR and Bezier
		std::vector<double> row = {dtbc, delta.get(0), delta.get(1), delta.get(2)};
		HR_bezier_delta.push_back(row);
	}
	
	// Create impulse function directory if it does not exist
	check_direc(bin_rn_direc + "/impulse_functions");
	
	// Save delta time series
	save_data(HR_bezier_delta, "Impulse function for HR",
			  bin_rn_direc + "/impulse_functions/delta_" + std::to_string(start_bin) + ".txt");

	return;
}


// Generates impulse function for entire cupolet from control sequence and starting bin.
// Reads in neuron state, dt, number of bins starting bin, control sequence, and save directories.
void save_impulse_series(hindmarsh_rose &neuron, double dt, unsigned int bins, unsigned int start_bin, std::vector<unsigned int> &ctrl, const std::string direc, const std::string bin_rn_direc) {
	// Read in the control plane initial conditions
	std::vector<array3> ps0inits, ps1inits, ps2inits;
	
	loadtxt_1(ps0inits, bin_rn_direc + "/coding_fcn/ps0_bin_inits.txt");
	loadtxt_1(ps1inits, bin_rn_direc + "/coding_fcn/ps1_bin_inits.txt");
	loadtxt_1(ps2inits, bin_rn_direc + "/coding_fcn/ps2b_bin_inits.txt");
	
	// Read in the vertices of each control plane
	std::vector<double> ps0x, ps0y, ps0z;
	std::vector<double> ps1x, ps1y, ps1z;
	std::vector<double> ps2x, ps2y, ps2z;
	
	loadtxt_3(ps0x, ps0y, ps0z, direc + "/control_planes/ps0_vertices.txt");
	loadtxt_3(ps1x, ps1y, ps1z, direc + "/control_planes/ps1_vertices.txt");
	loadtxt_3(ps2x, ps2y, ps2z, direc + "/control_planes/ps2_vertices.txt");
	
	// Read in the control plane endpoints for the polynomial fits
	std::vector<double> ps0endpts, ps1endpts;
	
	loadtxt_1(ps0endpts, bin_rn_direc + "/coding_fcn/ps0_bin_endpoints.txt");
	loadtxt_1(ps1endpts, bin_rn_direc + "/coding_fcn/ps1_bin_endpoints.txt");
	
	// Read in PS0/PS1 macrocontrol maps.
	std::vector<unsigned int> ps0macro, ps1macro;
	
	loadtxt_1(ps0macro, bin_rn_direc + "/macrocontrol/ps0_macrocontrol.txt");
	loadtxt_1(ps1macro, bin_rn_direc + "/macrocontrol/ps1_macrocontrol.txt");
	
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
	verts2b.push_back(array3(ps2x[2], ps2y[2], 0));
	verts2b.push_back(array3(ps2x[3], ps2y[3], 0));
	
	// Initialize impulse time series array (t,dx,dy,dz)
	std::vector<std::vector<double>> impulse_series;

	// Index of control
	unsigned int ci = 0;
	
	// Start with t = 0
	double t = 0;
	
	// Integrate until cupolet is fully mapped
	bool CONTINUE_INTEGRATING  = true;
	bool STOP_INTEGRATING_SOON = false;
	bool SAVE_IMPULSE_FUNCTION = true;
	
	// Used for split planes bezier curve
	double t0, t1;
	array3 P0, Q1;
	
	// Initial point dependent on which plane to start
	array3 curr = ps2inits[start_bin];
	
	// Store first entry into impulse series
	std::vector<double> row{0, 0, 0, 0};
	impulse_series.push_back(row);
	
	// Indicate which of the PS1 bins have been reached to detect if the cupolet has looped onto itself
	std::vector<bool> PS0_reached, PS1_reached;
	PS0_reached.resize(bins, false);
	PS1_reached.resize(bins, false);
	
	// Indicate that the starting bin on PS1 has been reached
	PS1_reached[start_bin] = true;
	
	// Cupolet name string
	std::string CUPOLET_NAME = "C";
	for (unsigned int i = 0; i < ctrl.size(); i++) {
		CUPOLET_NAME += (ctrl[i]) ? "1" : "0";
	}
	
	// Cupolet generation loop
	while (CONTINUE_INTEGRATING) {
		// Integrate one step forward
		array3 next = rk4(curr, dt, &hindmarsh_rose::hr_dynamics, neuron);
		
		// Check if PS0 has been crossed
		if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 0, verts0)) {
			// Find where trajectory intersects control plane
			std::vector<double> tyzxp = rk4_henon(array3(t, curr.get(1), curr.get(2)), curr.get(0), -(curr.get(0)-ps0x[0]), &hindmarsh_rose::hr_dx_dynamics, neuron);
			
			// Move time forward
			t = tyzxp[0];
			
			// Tick control index
			ci = (ci + 1) % ctrl.size();
			
			// Find the left side of the intersecting bin endpoint
			unsigned int ii = bisect(ps0endpts, tyzxp[1]) - 1;
			
			// Apply macrocontrol
			if (ctrl[ci]) {ii = ps0macro[ii];}
			
			
			// Store impulse into time series
			array3 impulse_array = ps0inits[ii] - next;
			row = {t, impulse_array.get(0), impulse_array.get(1), impulse_array.get(2)};
			impulse_series.push_back(row);
			
			// Center to bin
			next = ps0inits[ii];
			
			// If bin has been reached at the start of control sequence already,
			// set flag to false. Then mark bin as having been reached
			// If flag is false, then integration can stop
			CONTINUE_INTEGRATING = !PS0_reached[ii] || (ci > 0);
			PS0_reached[ii] = ci == 0;
			
		// Check if PS2a has been crossed
		} else if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 1, verts2a)) {
			// Find where trajectory intersects control plane
			std::vector<double> txzyp = rk4_henon(array3(t, curr.get(0), curr.get(2)), curr.get(1), -(curr.get(1)-ps2y[0]), &hindmarsh_rose::hr_dy_dynamics, neuron);
			
			// Move time forward
			t = txzyp[0];
			
			// Initial time step for bezier curve
			t0 = t;
			
			// Control points for bezier curve
			P0 = array3(txzyp[1], txzyp[3], txzyp[2]);
			Q1 = neuron.hr_dynamics(P0);
			
			// Snap point to PS2a
			next = P0;
			
			// Store impulse into time series
			//row = {0, 0, 0, 0};
			//impulse_series.push_back(row);
			
			// Impulse function will be replaced with Bezier curve,
			SAVE_IMPULSE_FUNCTION = false;
			
		// Check if PS1 has been crossed
		} else if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 1, verts1)) {
			// Find where trajectory intersects control plane
			std::vector<double> txzyp = rk4_henon(array3(t, curr.get(0), curr.get(2)), curr.get(1), -(curr.get(1)-ps1y[0]), &hindmarsh_rose::hr_dy_dynamics, neuron);
			
			// Move time forward
			t = txzyp[0];
			
			// Tick control index
			ci = (ci + 1) % ctrl.size();
			
			// Find the left side of the intersecting bin endpoint
			unsigned int ii = bisect(ps1endpts, txzyp[1]) - 1;
			
			// Apply macrocontrol
			if (ctrl[ci]) {ii = ps1macro[ii];}
			
			// Center to bin
			next = ps1inits[ii];
			
			// If bin has been reached at the start of control sequence already,
			// set flag to true. Then mark bin as having been reached
			// If flag is true, integration must continue until PS2b
			STOP_INTEGRATING_SOON = PS1_reached[ii] & (ci == 0);
			PS1_reached[ii] = ci == 0;
		
		// Check if PS2b has been crossed
		} else if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 1, verts2b)) {
			// Find where trajectory intersects control plane
			std::vector<double> txzyp = rk4_henon(array3(t, curr.get(0), curr.get(2)), curr.get(1), -(curr.get(1)-ps2y[2]), &hindmarsh_rose::hr_dy_dynamics, neuron);
			
			// Move time forward
			t = txzyp[0];
			
			// Final time step for bezier curve
			t1 = t;
			
			std::cout << t0 << std::endl;
			std::cout << t1 << std::endl;
			
			// If bin has been reached in PS1, stop integrating
			CONTINUE_INTEGRATING = !STOP_INTEGRATING_SOON;
			
			// Do bezier splining
			array3 P3 = array3(txzyp[1], txzyp[3], txzyp[2]);
			array3 Q2 = neuron.hr_dynamics(P3);
			
			// Finish doing bezier splining
			array3 P1 = P0 + Q1 / 3 * (t1 - t0);
			array3 P2 = P3 - Q2 / 3 * (t1 - t0);
			
			// Snap point to PS2b
			next = P3;
			
			// Initialize temporary curr array3 and time for Bezier curve
			array3 curr_bez = P0;
			double t_bez = t0;
			
			// Integrate from PS2a to PS2b and snap to Bezier curve at each step
			while (true) {
				// Integrate one step forward
				array3 next_bez = rk4(curr_bez, dt, &hindmarsh_rose::hr_dynamics, neuron);
				std::cout << next_bez << std::endl;
				
				// Check if PS2b has been crossed
				if (crossed(curr_bez.get(0), next_bez.get(0), curr_bez.get(1), next_bez.get(1), 1, verts2b)) {
					// Find where trajectory intersects control plane
					std::vector<double> txzyp_bez = rk4_henon(array3(t_bez, curr_bez.get(0), curr_bez.get(2)), curr_bez.get(1), -(curr_bez.get(1)-ps2y[2]), &hindmarsh_rose::hr_dy_dynamics, neuron);
					
					// Store impulse into time series
					row = {t1 - txzyp_bez[0], P3.get(0) - txzyp_bez[1], P3.get(1) - txzyp_bez[3], P3.get(2) - txzyp_bez[2]};
					impulse_series.push_back(row);
					
					// Save impulse function now
					SAVE_IMPULSE_FUNCTION = true;
					
					// Break out of loop
					break;
				
				// No plane has been crossed
				} else {
					// Move time forward
					t_bez += dt;
					
					// Normalized Bezier time
					double t_bez_norm = (t_bez - t0) / (t1 - t0);
					double t_bez_nrev = 1 - t_bez_norm;
					// Maybe normalize time so that PBC and next have the same y-coord?
					
					// Polynomials (1-t)^3, t(1-t)^2, t^2(1-t), and t^3
					double P0C = t_bez_nrev * t_bez_nrev * t_bez_nrev;
					double P1C = t_bez_nrev * t_bez_nrev * t_bez_norm;
					double P2C = t_bez_nrev * t_bez_norm * t_bez_norm;
					double P3C = t_bez_norm * t_bez_norm * t_bez_norm;
					
					// Corresponding point on Bezier curve
					array3 PBC = P0 * P0C + 3 * P1 * P1C + 3 * P2 * P2C + P3 * P3C;
					
					// Store impulse into time series
					row = {0, PBC.get(0) - next_bez.get(0), PBC.get(1) - next_bez.get(1), PBC.get(2) - next_bez.get(2)};
					impulse_series.push_back(row);
					
					// Set new values to old values
					curr_bez = PBC;
				}
			}
			
		// No plane has been crossed
		} else {
			// Move time forward
			t += dt;
			
			if (SAVE_IMPULSE_FUNCTION) {
				// Store impulse into time series
				std::vector<double> row{0, 0, 0, 0};
				impulse_series.push_back(row);
			}
		}
		
		// Set new values to old values
		curr = next;
	}
	
	save_data(impulse_series, "Impulse series data of cupolet " + CUPOLET_NAME, bin_rn_direc + "/" + CUPOLET_NAME + "_cupolet_impulse_series_data.txt");
	
	return;
}


// Generates impulse function for entire cupolet from control sequence and starting bin.
// Reads in neuron state, dt, number of bins starting bin, control sequence, and save directories.
// Generates cupolet from starting bin and impulse series.
void time_series_from_impulse(hindmarsh_rose &neuron, double dt, unsigned int bins, unsigned int start_bin, std::vector<unsigned int> &ctrl, const std::string direc, const std::string bin_rn_direc) {
	// Read in the control plane initial conditions
	std::vector<array3> ps2inits;
	
	loadtxt_1(ps2inits, bin_rn_direc + "/coding_fcn/ps2b_bin_inits.txt");
	
	// Read in the vertices of each control plane
	std::vector<double> ps0x, ps0y, ps0z;
	std::vector<double> ps2x, ps2y, ps2z;
	
	loadtxt_3(ps0x, ps0y, ps0z, direc + "/control_planes/ps0_vertices.txt");
	loadtxt_3(ps2x, ps2y, ps2z, direc + "/control_planes/ps2_vertices.txt");
	
	// Plane data used for checking crossing
	std::vector<array3> verts0;
	std::vector<array3> verts2a;
	std::vector<array3> verts2b;
	
	// Fill plane data with values
	verts0.push_back(array3(ps0x[0], ps0y[0], 0));
	verts0.push_back(array3(ps0x[1], ps0y[1], 0));
	
	verts2a.push_back(array3(ps2x[0], ps2y[0], 0));
	verts2a.push_back(array3(ps2x[1], ps2y[1], 0));
	verts2b.push_back(array3(ps2x[2], ps2y[2], 0));
	verts2b.push_back(array3(ps2x[3], ps2y[3], 0));
	
	// Initialize cupolet time series array (t,x,y,z)
	std::vector<std::vector<double>> time_series;
	
	// Start with t = 0
	double t = 0;
	
	// Initial point dependent on which plane to start
	array3 curr = ps2inits[start_bin];
	
	// Store first entry into time series
	std::vector<double> row{0, curr.get(0), curr.get(1), curr.get(2)};
	time_series.push_back(row);
	
	// Cupolet name string
	std::string CUPOLET_NAME = "C";
	for (unsigned int i = 0; i < ctrl.size(); i++) {
		CUPOLET_NAME += (ctrl[i]) ? "1" : "0";
	}
	
	// Impulse Seires Counter
	unsigned int ISC = 0;
	
	// Load impulse time series
	std::vector<double> time_impulse;
	std::vector<array3> xyz_impulse;
	
	loadtxt_2(time_impulse, xyz_impulse, bin_rn_direc + "/" + CUPOLET_NAME + "_cupolet_impulse_series_data.txt");
	
	// Cupolet generation loop
	while (ISC < time_impulse.size()) {
		// Integrate one step forward
		array3 next = rk4(curr, dt, &hindmarsh_rose::hr_dynamics, neuron);
		
		// Check if PS0 has been crossed
		if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 0, verts0)) {
			// Find where trajectory intersects control plane
			std::vector<double> tyzxp = rk4_henon(array3(t, curr.get(1), curr.get(2)), curr.get(0), -(curr.get(0)-ps0x[0]), &hindmarsh_rose::hr_dx_dynamics, neuron);
			
			// Move time forward
			t = tyzxp[0];
			
			// Snap point to PS0
			next = array3(tyzxp[3], tyzxp[1], tyzxp[2]);
			
			// Store position into time series
			row = {t, next.get(0), next.get(1), next.get(2)};
			time_series.push_back(row);
			
			// Apply impulse
			t += time_impulse[ISC];
			next += xyz_impulse[ISC];
			
			// Store position into time series
			row = {t, next.get(0), next.get(1), next.get(2)};
			time_series.push_back(row);
			
		// Check if PS2a has been crossed
		} else if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 1, verts2a)) {
			// Find where trajectory intersects control plane
			std::vector<double> txzyp = rk4_henon(array3(t, curr.get(0), curr.get(2)), curr.get(1), -(curr.get(1)-ps2y[0]), &hindmarsh_rose::hr_dy_dynamics, neuron);
			
			// Move time forward
			t = txzyp[0];
			
			// Snap point to PS2a
			next = array3(txzyp[1], txzyp[3], txzyp[2]);
			
			// Apply impulse
			t += time_impulse[ISC];
			next += xyz_impulse[ISC];
			
			// Store position into time series
			row = {t, next.get(0), next.get(1), next.get(2)};
			time_series.push_back(row);
			
		// Check if PS2b has been crossed
		} else if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 1, verts2b)) {
			// Find where trajectory intersects control plane
			std::vector<double> txzyp = rk4_henon(array3(t, curr.get(0), curr.get(2)), curr.get(1), -(curr.get(1)-ps2y[2]), &hindmarsh_rose::hr_dy_dynamics, neuron);
			
			// Move time forward
			t = txzyp[0];
			
			// Snap point to PS2b
			next = array3(txzyp[1], txzyp[3], txzyp[2]);
			
			// Apply impulse
			t += time_impulse[ISC];
			next += xyz_impulse[ISC];
			
			// Store position into time series
			row = {t, next.get(0), next.get(1), next.get(2)};
			time_series.push_back(row);
			
		// No plane has been crossed
		} else {
			// Move time forward
			t += dt;
			
			// Apply impulse
			t += time_impulse[ISC];
			next += xyz_impulse[ISC];
			
			// Store position into time series
			row = {t, next.get(0), next.get(1), next.get(2)};
			time_series.push_back(row);
		}
		
		// Set new values to old values
		curr = next;
		
		// Increment impulse counter
		ISC++;
	}
	
	save_data(time_series, "Time series data of cupolet " + CUPOLET_NAME + " generated by an impulse function", bin_rn_direc + "/" + CUPOLET_NAME + "_cupolet_from_impulse_data.txt");
	
	return;
}