// File: cupolet_from_impulse.cpp
// Purpose: Cupolet impulse methods implementation file
// Author: Daniel Bujno

#include "cupolet_from_impulse.h"
#include "control_planes.h"
#include "hindmarsh_rose.h"
#include "keep_data.h"
#include "helper.h"
#include "array3.h"
#include "rk4.h"

#include <string>
#include <vector>
#include <fstream>
#include <iostream>


// Generates impulse function for entire cupolet from control sequence and starting bin.
// Reads in neuron state, dt, number of bins starting bin, control sequence, and save directories.
// Series is time-aligned with multiples of dt.
void create_impulse_function(hindmarsh_rose &neuron, double dt, unsigned int bins, unsigned int start_bin, std::vector<unsigned int> &ctrl, const std::string direc, const std::string bin_rn_direc) {
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
	
	// Initialize impulse time series array (dx,dy,dz)
	std::vector<array3> impulse_series;

	// Index of control
	unsigned int ci = 0;
	
	// Start with t = 0
	double t = 0;
	
	// Integrate flags
	bool CONTINUE_INTEGRATING  = true;  // While loop flag
	bool STOP_INTEGRATING_SOON = false; // Indicate if last Bezier curve is being saved
	bool SAVE_IMPULSE_FUNCTION = true;  // Don't save impulse during Bezier transition
	bool SKIP_PS2B_CTRL_PLANE  = true;  // Skip considering first PS2b crossing (Crossing implies Bezier transition)
	
	// Used for split planes bezier curve
	// t0 and t1 are double-precise times when cupolet intersects PS2a and PS2b
	// curr_pre_bezier is the curent position immediately before PS2a
	// t0a is the time at which curr_pre_bezier is saved
	double t0, t1, t0a;
	array3 P0, Q1, curr_pre_bezier;
	
	// Initial point dependent on which plane to start
	array3 curr = ps1inits[start_bin];
	
	// Store first entry into impulse series
	impulse_series.push_back(array3(0,0,0));
	
	// Indicate which of the PS1 bins have been reached to detect if the cupolet has looped onto itself
	std::vector<bool> PS0_reached, PS1_reached;
	PS0_reached.resize(bins, false);
	PS1_reached.resize(bins, false);
	
	// Indicate that the starting bin on PS1 has been reached
	PS1_reached[start_bin] = true;
	
	// Cupolet name string
	std::string cupolet_name = "C";
	for (unsigned int i = 0; i < ctrl.size(); i++) {
		cupolet_name += (ctrl[i]) ? "1" : "0";
	}
	
	// Cupolet generation loop
	while (CONTINUE_INTEGRATING) {
		// Integrate one step forward
		array3 next = rk4(curr, dt, &hindmarsh_rose::hr_dynamics, neuron);
		
		// Check if PS0 has been crossed
		if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 0, verts0)) {
			// Find where trajectory intersects control plane
			std::vector<double> tyzxp = rk4_henon(array3(t, curr.get(1), curr.get(2)), curr.get(0), -(curr.get(0)-ps0x[0]), &hindmarsh_rose::hr_dx_dynamics, neuron);
			
			// Tick control index
			ci = (ci + 1) % ctrl.size();
			
			// Find the left side of the intersecting bin endpoint
			unsigned int ii = bisect(ps0endpts, tyzxp[1]) - 1;
			
			// Apply macrocontrol
			if (ctrl[ci]) {ii = ps0macro[ii];}
			
			// Integrate from center of bin
			next = rk4(ps0inits[ii], t - tyzxp[0] + dt, &hindmarsh_rose::hr_dynamics, neuron);
			
			// Store impulse into time series
			array3 impulse_array = rk4_reverse(curr, next, dt, &hindmarsh_rose::hr_dynamics, neuron);
			impulse_series.push_back(impulse_array);
			
			// If bin has been reached at the start of control sequence already,
			// set flag to false. Then mark bin as having been reached
			// If flag is false, then integration can stop
			CONTINUE_INTEGRATING = !PS0_reached[ii] || (ci > 0);
			PS0_reached[ii] = ci == 0;
			
		// Check if PS2a has been crossed
		} else if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 1, verts2a)) {
			// Find where trajectory intersects control plane
			std::vector<double> txzyp = rk4_henon(array3(t, curr.get(0), curr.get(2)), curr.get(1), -(curr.get(1)-ps2y[0]), &hindmarsh_rose::hr_dy_dynamics, neuron);
			
			// Initial time step for bezier curve
			t0 = txzyp[0];
			
			// Control points and pseudo control points for bezier curve
			P0 = array3(txzyp[1], txzyp[3], txzyp[2]);
			Q1 = neuron.hr_dynamics(P0);
			
			// Save curr for use in Bezier curve
			curr_pre_bezier = curr;
			t0a = t;
			
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
			// Skip first crossing of PS2b
			if (SKIP_PS2B_CTRL_PLANE) {
				SKIP_PS2B_CTRL_PLANE = false;
				
				// Set new values to old values
				curr = next;
				
				continue;
			}
			
			// Find where trajectory intersects control plane
			std::vector<double> txzyp = rk4_henon(array3(t, curr.get(0), curr.get(2)), curr.get(1), -(curr.get(1)-ps2y[2]), &hindmarsh_rose::hr_dy_dynamics, neuron);
			
			// Final time step for bezier curve
			t1 = txzyp[0];
			
			// If bin has been reached in PS1, stop integrating
			CONTINUE_INTEGRATING = !STOP_INTEGRATING_SOON;
			
			// Control points and pseudo control points for bezier curve
			array3 P3 = array3(txzyp[1], txzyp[3], txzyp[2]);
			array3 Q2 = neuron.hr_dynamics(P3);
			array3 P1 = P0 + Q1 / 3 * (t1 - t0);
			array3 P2 = P3 - Q2 / 3 * (t1 - t0);
			
			// Calculate time until dt alignment
			double ta = t1;
			while (ta > dt) {ta -= dt;}
			
			// Next dt-alligned position
			next = rk4(P3, dt - ta, &hindmarsh_rose::hr_dynamics, neuron);
			
			// Initialize temporary curr array3 and time for Bezier curve
			array3 curr_bez = curr_pre_bezier;
			double t_bez = t0a;
			
			// Integrate from PS2a to PS2b and snap to Bezier curve at each step
			while (true) {
				// Integrate one step forward
				array3 next_bez = rk4(curr_bez, dt, &hindmarsh_rose::hr_dynamics, neuron);
				
				// Check if PS2b has been crossed
				if (crossed(curr_bez.get(0), next_bez.get(0), curr_bez.get(1), next_bez.get(1), 1, verts2b)) {
					// Store impulse into time series
					array3 impulse_array = rk4_reverse(curr_bez, next, dt, &hindmarsh_rose::hr_dynamics, neuron);
					impulse_series.push_back(impulse_array);
					
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
					array3 impulse_array = rk4_reverse(curr_bez, PBC, dt, &hindmarsh_rose::hr_dynamics, neuron);
					impulse_series.push_back(impulse_array);
					
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
				impulse_series.push_back(array3(0,0,0));
			}
		}
		
		// Set new values to old values
		curr = next;
	}
	
	save_data(impulse_series, "Impulse series data of cupolet " + cupolet_name, bin_rn_direc + "/" + cupolet_name + "_cupolet_impulse_series_data.txt");
	
	return;
}

// Generates cupolet from impulse function.// Generates cupolet from starting bin and impulse series.
// Reads in neuron state, dt, number of bins starting bin, control sequence, and save directories.
void create_time_series_from_impulse(hindmarsh_rose &neuron, double dt, unsigned int bins, unsigned int start_bin, std::vector<unsigned int> &ctrl, const std::string direc, const std::string bin_rn_direc) {
	// Read in the control plane initial conditions
	std::vector<array3> ps1inits;
	
	loadtxt_1(ps1inits, bin_rn_direc + "/coding_fcn/ps1_bin_inits.txt");

	// Initialize cupolet time series array (t,x,y,z)
	std::vector<std::vector<double>> time_series;
	
	// Start with t = 0
	double t = 0;
	
	// Initial point dependent on which plane to start
	array3 curr = ps1inits[start_bin];
	
	// Store first entry into time series
	std::vector<double> row{0, curr.get(0), curr.get(1), curr.get(2)};
	time_series.push_back(row);
	
	// Cupolet name string
	std::string cupolet_name = "C";
	for (unsigned int i = 0; i < ctrl.size(); i++) {
		cupolet_name += (ctrl[i]) ? "1" : "0";
	}
	
	// Impulse Seires Counter
	unsigned int ISC = 1; // Off by 1 error
	
	// Load impulse time series
	std::vector<array3> xyz_impulse;
	
	loadtxt_1(xyz_impulse, bin_rn_direc + "/" + cupolet_name + "_cupolet_impulse_series_data.txt");
	
	// Cupolet generation loop
	while (ISC < xyz_impulse.size()) {
		// Integrate one step forward and apply impulse
		array3 next = rk4_impulse(curr, xyz_impulse[ISC], dt, &hindmarsh_rose::hr_dynamics, neuron);
		
		// Move time forward
		t += dt;
		
		// Store position into time series
		row = {t, next.get(0), next.get(1), next.get(2)};
		time_series.push_back(row);
		
		// Set new values to old values
		curr = next;
		
		// Increment impulse counter
		ISC++;
	}
	
	save_data(time_series, "Time series data of cupolet " + cupolet_name + " generated by an impulse function", bin_rn_direc + "/" + cupolet_name + "_cupolet_from_impulse_data.txt");
	
	return;
}

// Generates cupolet and impulse function from control sequence and starting bin. Runs stabilized cupolet for N times.
// Reads in neuron state, dt, number of bins starting bin, control sequence, and save directories.
// Series is time-aligned with multiples of dt.
void create_cupolet(hindmarsh_rose &neuron, double dt, unsigned int bins, unsigned int start_bin, unsigned int N, std::vector<unsigned int> &ctrl, const std::string direc, const std::string bin_rn_direc, bool save_mode) {
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
	
	// Cupolet name string
	std::string cupolet_name = "C";
	for (unsigned int i = 0; i < ctrl.size(); i++) {
		cupolet_name += (ctrl[i]) ? "1" : "0";
	}
	
	// Files for impulse, offsets, and time series
	std::ofstream impulse_file, offsets_file, cupolet_file;
	
	// Whether to use modified save directory or original save directory.
	if (save_mode) {
		// Save directory prefixes and postfixes
		std::string save_dir = bin_rn_direc + "/fundamental_cupolets/" + cupolet_name + "/bin_" + std::to_string(start_bin) + "/text_data/";
		std::string prefix   = save_dir + cupolet_name + "_bin_" + std::to_string(start_bin) + "_";
		std::string postfix  = "_";
		if (N < 10) {postfix += "0";}
		postfix += std::to_string(N) + "_Period.txt";
		
		// Save locations for impulse, offsets, and time series
		impulse_file.open(prefix + "Impulse" + postfix);
		offsets_file.open(prefix + "Offsets" + postfix);
		cupolet_file.open(prefix + "Cupolet" + postfix);
	} else {
		// Save locations for impulse, offsets, and time series
		impulse_file.open(bin_rn_direc + "/cupolet_impulse.txt");
		offsets_file.open(bin_rn_direc + "/cupolet_offsets.txt");
		cupolet_file.open(bin_rn_direc + "/cupolet_time_series.txt");
	}
	
	// Save header for each file
	impulse_file << "# " << "Impulse series data of cupolet " + cupolet_name << std::endl;
	offsets_file << "# " << "Impulse offsets of cupolet " + cupolet_name << std::endl;
	cupolet_file << "# " << "Time series data of cupolet " + cupolet_name << std::endl;
	
	// Index of control
	unsigned int ci = 0;
	
	// Start with t = 0
	double t = 0;
	
	// Integrate flags
	bool CONTINUE_INTEGRATING  = true;  // While loop flag
	bool STOP_INTEGRATING_SOON = false; // Indicate if last Bezier curve is being saved
	bool SAVE_IMPULSE_CUPOLET  = true;  // Don't save impulse during Bezier transition
	bool SKIP_PS2B_CTRL_PLANE  = true;  // Skip considering first PS2b crossing (Crossing implies Bezier transition)
	
	// Used for split planes bezier curve
	// t0 and t1 are double-precise times when cupolet intersects PS2a and PS2b
	// curr_pre_bezier is the curent position immediately before PS2a
	// t0a is the time at which curr_pre_bezier is saved
	double t0, t1, t0a;
	array3 P0, Q1, curr_pre_bezier;
	
	// Initial point dependent on which plane to start
	array3 curr = ps1inits[start_bin];
	
	// Store first entry into cupolet and impulse series
	save_datum(0, curr, cupolet_file);
	save_datum(0, 0, 0, impulse_file);
	
	// Indicate how much each of the PS1 bins have been reached to detect if the cupolet has looped onto itself N times
	std::vector<unsigned int> PS0_reached, PS1_reached;
	PS0_reached.resize(bins, 0);
	PS1_reached.resize(bins, 0);
	
	// Indicate that the starting bin on PS1 has been reached
	PS1_reached[start_bin] = 1;
	
	// Cupolet generation loop
	while (CONTINUE_INTEGRATING) {
		
		// Integrate one step forward
		array3 next = rk4(curr, dt, &hindmarsh_rose::hr_dynamics, neuron);
		
		// Check if PS0 has been crossed
		if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 0, verts0)) {
			// Find where trajectory intersects control plane
			std::vector<double> tyzxp = rk4_henon(array3(t, curr.get(1), curr.get(2)), curr.get(0), -(curr.get(0)-ps0x[0]), &hindmarsh_rose::hr_dx_dynamics, neuron);
			
			// Tick control index
			ci = (ci + 1) % ctrl.size();
			
			// Find the left side of the intersecting bin endpoint
			unsigned int ii = bisect(ps0endpts, tyzxp[1]) - 1;
			
			// Apply macrocontrol
			if (ctrl[ci]) {ii = ps0macro[ii];}
			
			// Move time forward
			t += dt;
			
			// Integrate from center of bin
			next = rk4(ps0inits[ii], t - tyzxp[0], &hindmarsh_rose::hr_dynamics, neuron);
			
			// Store offsets, impulse and cupolet position
			array3 impulse_array = rk4_reverse(curr, next, dt, &hindmarsh_rose::hr_dynamics, neuron);
			save_datum(t, next, cupolet_file);
			save_datum(impulse_array, impulse_file);
			save_datum(dt + tyzxp[0] - t, offsets_file);
			
			
			// If bin has been reached at the start of control sequence already,
			// set flag to false. Then increment bin counter
			// If flag is false, then integration can stop
			CONTINUE_INTEGRATING = (PS0_reached[ii] < N) || (ci > 0);
			PS0_reached[ii] += ci == 0;
			
		// Check if PS2a has been crossed
		} else if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 1, verts2a)) {
			// Find where trajectory intersects control plane
			std::vector<double> txzyp = rk4_henon(array3(t, curr.get(0), curr.get(2)), curr.get(1), -(curr.get(1)-ps2y[0]), &hindmarsh_rose::hr_dy_dynamics, neuron);
			
			// Initial time step for bezier curve
			t0 = txzyp[0];
			
			// Control points and pseudo control points for bezier curve
			P0 = array3(txzyp[1], txzyp[3], txzyp[2]);
			Q1 = neuron.hr_dynamics(P0);
			
			// Save curr for use in Bezier curve
			curr_pre_bezier = curr;
			t0a = t;
			
			// Save Bezier impulse offsets
			save_datum(t0 - t, offsets_file);
			
			// cupolet will be replaced with Bezier curve,
			SAVE_IMPULSE_CUPOLET = false;
			
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
			STOP_INTEGRATING_SOON = (PS1_reached[ii] >= N) & (ci == 0);
			PS1_reached[ii] += ci == 0;
		
		// Check if PS2b has been crossed
		} else if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 1, verts2b)) {
			// Skip first crossing of PS2b
			if (SKIP_PS2B_CTRL_PLANE) {
				SKIP_PS2B_CTRL_PLANE = false;
				
				// Set new values to old values
				curr = next;
				
				// Move time forward
				t += dt;
				
				// Store impulse and cupolet position
				save_datum(t, next, cupolet_file);
				save_datum(0, 0, 0, impulse_file);
				
				continue;
			}
			
			// Find where trajectory intersects control plane
			std::vector<double> txzyp = rk4_henon(array3(t, curr.get(0), curr.get(2)), curr.get(1), -(curr.get(1)-ps2y[2]), &hindmarsh_rose::hr_dy_dynamics, neuron);
			
			// Final time step for bezier curve
			t1 = txzyp[0];
			
			// If bin has been reached in PS1, stop integrating
			CONTINUE_INTEGRATING = !STOP_INTEGRATING_SOON;
			
			// Control points and pseudo control points for bezier curve
			array3 P3 = array3(txzyp[1], txzyp[3], txzyp[2]);
			array3 Q2 = neuron.hr_dynamics(P3);
			array3 P1 = P0 + Q1 / 3 * (t1 - t0);
			array3 P2 = P3 - Q2 / 3 * (t1 - t0);
			
			// Calculate time until dt alignment
			double ta = t1;
			while (ta > dt) {ta -= dt;}
			
			// Next dt-alligned position
			next = rk4(P3, dt - ta, &hindmarsh_rose::hr_dynamics, neuron);
			
			// Align t
			t = t1 + dt - ta;
			
			// Initialize temporary curr array3 and time for Bezier curve
			array3 curr_bez = curr_pre_bezier;
			double t_bez = t0a;
			
			// Integrate from PS2a to PS2b and snap to Bezier curve at each step
			while (true) {
				// Integrate one step forward
				array3 next_bez = rk4(curr_bez, dt, &hindmarsh_rose::hr_dynamics, neuron);
				if (curr_bez.get(0) > 3) {break;}
				// Move time forward
				t_bez += dt;
				
				// Check if PS2b has been crossed
				if (t_bez >= t || crossed(curr_bez.get(0), next_bez.get(0), curr_bez.get(1), next_bez.get(1), 1, verts2b)) {
					// Store impulse and cupolet position
					array3 impulse_array = rk4_reverse(curr_bez, next, dt, &hindmarsh_rose::hr_dynamics, neuron);
					save_datum(t_bez, next, cupolet_file);
					save_datum(impulse_array, impulse_file);
					
					// Save impulse and cupolet now
					SAVE_IMPULSE_CUPOLET = true;
					
					// Break out of loop
					break;
				
				// No plane has been crossed
				} else {
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
					
					// Store impulse and cupolet position
					array3 impulse_array = rk4_reverse(curr_bez, PBC, dt, &hindmarsh_rose::hr_dynamics, neuron);
					save_datum(t_bez, PBC, cupolet_file);
					save_datum(impulse_array, impulse_file);
					
					// Set new values to old values
					curr_bez = PBC;
				}
			}
			
		// No plane has been crossed
		} else {
			// Move time forward
			t += dt;
			
			if (SAVE_IMPULSE_CUPOLET) {
				// Store impulse and cupolet position
				save_datum(t, next, cupolet_file);
				save_datum(0, 0, 0, impulse_file);
			}
		}
		
		// Set new values to old values
		curr = next;
	}
	
	// Close files
	impulse_file.close();
	offsets_file.close();
	cupolet_file.close();
	
	return;
}

// Generates impulse curves along a given bin preimage from PS1 onto PS2a. N+1 points are sampled between PS2a and PS2b inclusive
// Reads in neuron state, dt, number of bins, starting bin, whether to kick, and save directories.
void bezier_impulse(hindmarsh_rose &neuron, double dt, int bins, unsigned int start_bin, unsigned int N, bool kick, const std::string direc, const std::string bin_rn_direc) {
	// Read in the PS2a control plane initial conditions
	std::vector<array3> ps2ainits;
	loadtxt_1(ps2ainits, bin_rn_direc + "/coding_fcn/ps2a_bin_inits.txt");
	
	// Read in the vertices of each control plane
	std::vector<double> ps2x, ps2y, ps2z;
	loadtxt_3(ps2x, ps2y, ps2z, direc + "/control_planes/ps2_vertices.txt");
	
	// Read in PS0/PS1 macrocontrol maps.
	std::vector<unsigned int> ps1macro;
	
	loadtxt_1(ps1macro, bin_rn_direc + "/macrocontrol/ps1_macrocontrol.txt");
	
	// Files for impulse, offsets, and time series
	std::ofstream impulse_file, impulse_file;
	
	// Is kick performed?
	if (kick) {
		// Open file
		impulse_file.open(bin_rn_direc + "/bezier_impulses/kick/bin_" + std::to_string(start_bin) + "_kick.txt");
		
		// Save header for each file
		impulse_file << "# " << "Bezier impulse with kick starting on bin" + cupolet_name << std::endl;
		
	} else {
		// Open file
		impulse_file.open(bin_rn_direc + "/bezier_impulses/no_kick/bin_" + std::to_string(start_bin) + "_no_kick.txt");
		
		// Save header for each file
		impulse_file << "# " << "Bezier impulse without kick starting on bin" + cupolet_name << std::endl;
	}
	
	// Start with t = 0
	double t = 0;
	
	// Bezier curve data
	double t_delta;
	array3 P1, P2, P3, Q2;
	array3 P0 = ps2ainits[start_bin];
	array3 Q1 = neuron.hr_dynamics(P0);
	
	// Initial point on PS2a
	// Preimages integrate onto PS1 inits, at which point a kick can be performed
	// Kick can be replicated by kicking the starting bin
	array3 curr = kick ? ps2ainits[ps1macro[start_bin]] : ps2ainits[start_bin];
	
	// Integrate until corresponding point on PS2b is reached.
	while (true) {
		// Integrate one step forward
		array3 next = rk4(curr, dt, &hindmarsh_rose::hr_dynamics, neuron);
		
		// Check if PS2b has been crossed
		if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 1, verts2b)) {
			// Find where trajectory intersects control plane
			std::vector<double> txzyp = rk4_henon(array3(t, curr.get(0), curr.get(2)), curr.get(1), -(curr.get(1)-ps2y[2]), &hindmarsh_rose::hr_dy_dynamics, neuron);
			
			// Save finishing time
			t_delta = txzyp[0];
			
			// Control points and pseudo control points for bezier curve
			P3 = array3(txzyp[1], txzyp[3], txzyp[2]);
			Q2 = neuron.hr_dynamics(P3);
			P1 = P0 + Q1 / 3 * t_delta;
			P2 = P3 - Q2 / 3 * t_delta;
			
			break;
		}
		
		// Increment time
		t += dt;
	}
	
	// Starting point for computing impulse
	array3 start = rk4_undo(P0, dt, &hindmarsh_rose::hr_dynamics, neuron);
	
	// Time difference between starting point and ending point
	double t_span = t_delta + dt;
	
	// Compute impulse function
	for (unsigned int k = 0; k <= N; k++) {
		// Time along bezier transition
		// A negative time is a point before the bezier transition
		// A time within dt of t_delta will integrate past the bezier transition
		double t_rel = t_span * k / N - dt;
		
		// Start and end point for determining impulse
		array3 R0, R1;
		
		// If t_rel is negative, R0 is not on bezier transition
		if (t_rel < 0) {
			R0 = rk4(start, dt + t_rel, &hindmarsh_rose::hr_dynamics, neuron);
		} else {
			// Normalized Bezier time
			double t_bez_norm = t_rel / t_delta;
			double t_bez_nrev = 1 - t_bez_norm;
			// Maybe normalize time so that PBC and next have the same y-coord?
			
			// Polynomials (1-t)^3, t(1-t)^2, t^2(1-t), and t^3
			double P0C = t_bez_nrev * t_bez_nrev * t_bez_nrev;
			double P1C = t_bez_nrev * t_bez_nrev * t_bez_norm;
			double P2C = t_bez_nrev * t_bez_norm * t_bez_norm;
			double P3C = t_bez_norm * t_bez_norm * t_bez_norm;
			
			// Corresponding point on Bezier curve
			R0 = P0 * P0C + 3 * P1 * P1C + 3 * P2 * P2C + P3 * P3C;
		}
		
		// If t_rel is within dt of t_delta, R1 is not on bezier transition
		if (t_rel + dt > t_delta) {
			R1 = rk4(P3, t_rel + dt - t_delta, &hindmarsh_rose::hr_dynamics, neuron);
		} else {
			// Normalized Bezier time
			double t_bez_norm = (t_rel + dt) / t_delta;
			double t_bez_nrev = 1 - t_bez_norm;
			// Maybe normalize time so that PBC and next have the same y-coord?
			
			// Polynomials (1-t)^3, t(1-t)^2, t^2(1-t), and t^3
			double P0C = t_bez_nrev * t_bez_nrev * t_bez_nrev;
			double P1C = t_bez_nrev * t_bez_nrev * t_bez_norm;
			double P2C = t_bez_nrev * t_bez_norm * t_bez_norm;
			double P3C = t_bez_norm * t_bez_norm * t_bez_norm;
			
			// Corresponding point on Bezier curve
			R1 = P0 * P0C + 3 * P1 * P1C + 3 * P2 * P2C + P3 * P3C;
		}
		
		// Store impulse
		array3 impulse = rk4_reverse(R0, R1, dt, &hindmarsh_rose::hr_dynamics, neuron);
		save_datum(t_rel, impulse_array, impulse_file);
	}
	
	// Close files
	impulse_file.close();
	
	return;
}