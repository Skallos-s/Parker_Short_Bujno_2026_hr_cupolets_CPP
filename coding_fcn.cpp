// File: coding_fcn.cpp
// Purpose: Implementation file for Coding Functions methods
// Author: Daniel Bujno


#include "hindmarsh_rose.h"
#include "control_planes.h"
#include "coding_fcn.h"
#include "keep_data.h"
#include "helper.h"
#include "rk4.h"

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>

// Read in HR neuron, dt, direc directory to store files, the number of iterations for the coding
// function to store and the number of bins to parse each poincare surface
void coding_fcn(hindmarsh_rose &neuron, const std::string direc, const std::string bin_rn_direc, double dt, unsigned int coding_fcn_N, unsigned int bins, bool generate_poly) {
	
	std::cout << "Generating coding functions..." << std::endl;
	
	// Directory to store all coding function data
	const std::string store_direc = bin_rn_direc + "/coding_fcn";
	
	// Create coding function directory if it does not exist
	check_direc(store_direc);
	
	// Read in PS1 and PS0 vertices
	std::vector<double> ps0x, ps0y, ps0z;
	std::vector<double> ps1x, ps1y, ps1z;
	
	loadtxt_3(ps0x, ps0y, ps0z, direc + "/control_planes/ps0_vertices.txt");
	loadtxt_3(ps1x, ps1y, ps1z, direc + "/control_planes/ps1_vertices.txt");
	
	// Read in PS1 and PS0 poincare points
	std::vector<double> ps0pts_t, ps0pts_x, ps0pts_y, ps0pts_z;
	std::vector<double> ps1pts_t, ps1pts_x, ps1pts_y, ps1pts_z;
	
	loadtxt_4(ps0pts_t, ps0pts_x, ps0pts_y, ps0pts_z, direc + "/control_planes/ps0_pts.txt");
	loadtxt_4(ps1pts_t, ps1pts_x, ps1pts_y, ps1pts_z, direc + "/control_planes/ps1_pts.txt");
	
	// Create polynomial fits for PS1 and PS0 pts
	std::vector<double> ps0_p;
	std::vector<double> ps1_p;
	if (generate_poly) {
		ps0_p = quadratic_regression(ps0pts_y, ps0pts_z);
		ps1_p = quadratic_regression(ps1pts_y, ps1pts_z); // Implement cubic regression and replace quadratic with cubic
		
		// Save the polynomial coefficients for each plane
		save_double_vector(ps0_p, "Coefficients of approximating polynomial for PS0 (first is coefficient of highest degree).", store_direc + "/ps0_poly_coeffs.txt");
		save_double_vector(ps1_p, "Coefficients of approximating polynomial for PS1 (first is coefficient of highest degree).", store_direc + "/ps1_poly_coeffs.txt");
	} else {
		loadtxt_1(ps0_p, store_direc + "/ps0_poly_coeffs.txt");
		loadtxt_1(ps1_p, store_direc + "/ps1_poly_coeffs.txt");
	}
	
	// Generate the end points, middle, and correspoding x y z values of each midpoint for both planes
	std::vector<double> bin0, bin1, bin0_mids, bin1_mids;
	std::vector<array3> bin0_inits, bin1_inits;
	gen_bin_data(bin0, bin0_mids, bin0_inits, ps0x, ps0y, ps0_p, bins, 0);
	gen_bin_data(bin1, bin1_mids, bin1_inits, ps1x, ps1y, ps1_p, bins, 1);
	
	// Save the end points of the bins for PS0 and PS1
	save_double_vector(bin0, "End points of each bin for PS0 (y values).", store_direc + "/ps0_bin_endpoints.txt");
	save_double_vector(bin1, "End points of each bin for PS1 (x values).", store_direc + "/ps1_bin_endpoints.txt");
	
	// Save the mid points of each of the bins for PS0 and PS1
	save_double_vector(bin0_mids, "Mid points of each bin for PS0 (y values).", store_direc + "/ps0_bin_midpoints.txt");
	save_double_vector(bin1_mids, "Mid points of each bin for PS1 (x values).", store_direc + "/ps1_bin_midpoints.txt");
	
	// Save the initial points of the bins for PS0 and PS1 (x y z values corresponding to the middle of the bin on the planes)
	save_array3_vector(bin0_inits, "Initial conditions of each bin for PS0 (x y z values).", store_direc + "/ps0_bin_inits.txt");
	save_array3_vector(bin1_inits, "Initial conditions of each bin for PS1 (x y z values).", store_direc + "/ps1_bin_inits.txt");
	
	std::vector<double> rn0, rn1;
	std::vector<std::vector<unsigned int>> map0, map1;
	
	// Determine the coding sequence and corresponding binary representation of middle point on each bin
	gen_crossing_sequence(rn0, map0, bin0_inits, neuron, ps0x, ps0y, ps1x, ps1y, coding_fcn_N, bins, dt);
	gen_crossing_sequence(rn1, map1, bin1_inits, neuron, ps0x, ps0y, ps1x, ps1y, coding_fcn_N, bins, dt);
	
	// Save binary decmial for each coding map for each control plane
	save_double_vector(rn0, "Rn function for PS0 (binary decimal).", store_direc + "/ps0_rn_fcn.txt");
	save_double_vector(rn1, "Rn function for PS1 (binary decimal).", store_direc + "/ps1_rn_fcn.txt");
	
	// Save each coding map for each control plane
	save_binary(map0, "Symbolic future crossings of PS0 for each bin.", store_direc + "/ps0_rn_map_fcn.txt");
	save_binary(map1, "Symbolic future crossings of PS1 for each bin.", store_direc + "/ps1_rn_map_fcn.txt");
	
	return;
}

// For surface PS reads in x vertex values psx, y vertex values psy and polynmoial fit poly_coords,
// and bins to generate the bin spacing, midpoints, and initial data for evaluation 
// write output to bin, bin_mids, and bin_inits
void gen_bin_data(std::vector<double> &rbin_endpts, std::vector<double> &rbin_mids, std::vector<array3> &rbin_inits, std::vector<double> &psx, std::vector<double> &psy, std::vector<double> &poly_coords, unsigned int bins, int ps) {
	if (ps == 0) {
		
		// The edges of all the bins
		double DELTA = (psy[1] - psy[0]) / bins;
		for (unsigned int i = 0; i <= bins; i++) {
			rbin_endpts.push_back(psy[0] + DELTA * i);
		}
		
		// Half the distance between two bins
		double mid_db = (rbin_endpts[1] - rbin_endpts[0]) / 2;
		
		// All the middles of the bins
		DELTA = (rbin_endpts[rbin_endpts.size() - 1] - rbin_endpts[0] - 2 * mid_db) / (bins - 1);
		for (unsigned int i = 0; i < bins; i++) {
			rbin_mids.push_back(rbin_endpts[0] + mid_db + DELTA * i);
		}
		
		// Initial conditions (x y z) corresponding to each bin
		for (unsigned int i = 0; i < bins; i++) {
			rbin_inits.push_back(array3(psx[0], rbin_mids[i], eval_poly(poly_coords,rbin_mids[i])));
		}
		
	} else {
		// The edges of all the bins
		double DELTA = (psx[1] - psx[0]) / bins;
		for (unsigned int i = 0; i <= bins; i++) {
			rbin_endpts.push_back(psx[0] + DELTA * i);
		}
		
		// Half the distance between two bins
		double mid_db = (rbin_endpts[1] - rbin_endpts[0]) / 2;
		
		// All the middles of the bins
		DELTA = (rbin_endpts[rbin_endpts.size() - 1] - rbin_endpts[0] - 2 * mid_db) / (bins - 1);
		for (unsigned int i = 0; i < bins; i++) {
			rbin_mids.push_back(rbin_endpts[0] + mid_db + DELTA * i);
		}
		
		// Initial conditions (x y z) corresponding to each bin
		for (unsigned int i = 0; i < bins; i++) {
			rbin_inits.push_back(array3(rbin_mids[i], psy[0], eval_poly(poly_coords,rbin_mids[i])));
		}
	}
	
	return;
}

// Reads in list of initial bins x, crossings N, number of bins, dt, vertex data of PS0 and PS1 (ps0x, ps0y, ps1x, ps1y)
// Assigns decimal value of binary visitation sequence and corresponding symbolic dynamics for each bin to map and rmap
// map is a vector of length (bins) containing vectors of length (N)
// rmap is a vector of length (bins)
// map will contain the binary sequence of the next N control planes
// rmap will contain the binary decimal that map produces
void gen_crossing_sequence(std::vector<double> &rmap, std::vector<std::vector<unsigned int>> &map, std::vector<array3> &x, hindmarsh_rose &neuron, std::vector<double> &ps0x, std::vector<double> &ps0y, std::vector<double> &ps1x, std::vector<double> &ps1y, unsigned int N, unsigned int bins, double dt) {
	// Initialize rmap and map to 0
	for (unsigned int i = 0; i < bins; i++) {
		// Temporary row for map
		std::vector<unsigned int> row;
		
		// Initialize row
		for (unsigned int j = 0; j < N; j++) {
			row.push_back(0);
		}
		
		// Initialize rmap
		rmap.push_back(0);
		
		// Initialize map
		map.push_back(row);
	}
	
	// Temporary initialization
	std::vector<array3> verts0;
	std::vector<array3> verts1;
		
	// Fill temp vector with values
	verts0.push_back(array3(ps0x[0], ps0y[0], 0));
	verts0.push_back(array3(ps0x[1], ps0y[1], 0));
	verts1.push_back(array3(ps1x[0], ps1y[0], 0));
	verts1.push_back(array3(ps1x[1], ps1y[1], 0));
	
	// Loop through each initial condition of each bin
	for (unsigned int i = 0; i < bins; i++) {
		// Start with 0 crossings and the initial condition is the middle of the bin
		unsigned int cross = 0;
		array3 curr = x[i];
		
		// Continue collecting crossings until N are found
		while (cross < N) {
			// Integrate one step forward
			array3 next = rk4(curr, dt, &hindmarsh_rose::hr_dynamics, neuron);
			
			// Check if PS0 has been crossed
			if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 0, verts0)) {
				// Add a 0 to the symbolic dynamics
				map[i][cross] = 0;
				
				// Increase the number of crossings
				cross += 1;
				
				// Check if PS1 has been crossed
			} else if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 1, verts1)) {
				// Add a 1 to the symbolic dynamics
				map[i][cross] = 1;
				
				// Increase the number of crossings
				cross += 1;
			}
			
			// New data becomes old data
			curr = next;
		}
		
		// Set rmap value
		rmap[i] = 0;
		for (unsigned int j = 0; j < N; j++) {
			double val = map[i][j];
			rmap[i] += val / (1 << (j + 1));
		}
	}
	
	return;
}


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
	std::vector<double> ps2x, ps2y, ps2z;
	loadtxt_3(ps2x, ps2y, ps2z, direc + "/control_planes/ps2_vertices.txt");
	
	// Plane data used for checking crossing
	std::vector<array3> verts2b;
	
	// Fill plane data with values
	verts2b.push_back(array3(ps2x[2], ps2y[2], ps2z[2]));
	verts2b.push_back(array3(ps2x[3], ps2y[3], ps2z[3]));
	
	// PS2b bin initial points
	std::vector<array3> ps2b_inits;
	
	// Travel time between PS1 and PS2b
	std::vector<double> ps2b_time;
	
	// Previous point prior to crossing PS2b
	std::vector<array3> ps2b_prev;
	
	// Step forward ps1inits to find ps2b_inits
	for (unsigned int i = 0; i < bins; i++) {
		// Initial point 
		array3 curr = ps1inits[i];
		
		// Contiue to integrate system until control plane is crossed
		bool crossed_bool = false;
		
		// Local time variable
		double t = 0;
		
		while (not crossed_bool) {
			// Integrate one step forward
			array3 next = rk4(curr, dt, &hindmarsh_rose::hr_dynamics, neuron);
			
			// Check if PS1 has been crossed
			if (crossed(curr.get(0), next.get(0), curr.get(1), next.get(1), 1, verts2b)) {
				// Find where trajectory intersects control plane
				std::vector<double> txzyp = rk4_henon(array3(t, curr.get(0), curr.get(2)), curr.get(1), -(curr.get(1)-ps2y[2]), &hindmarsh_rose::hr_dy_dynamics, neuron);
				
				// Mark as crossed to move to next bin
				crossed_bool = true;
				
				// Move time forward
				t = txzyp[0];
				
				// Store ending position into ps2b_inits
				ps2b_inits.push_back(array3(txzyp[1], txzyp[3], txzyp[2]));
				
				// Store travel time into ps2b_time
				ps2b_time.push_back(t);
				
				// Store previous position into ps2b_prev;
				ps2b_prev.push_back(curr);
			} else {
				// Increment time
				t += dt;
			}
			
			// Set new values to old values
			curr = next;
			
		}
	}
	
	// Save the initial points of the PS2b virtual bins (x y z values corresponding to the middle of the bin on the planes)
	save_array3_vector(ps2b_inits, "Initial conditions of each bin for PS2b (x y z values).", store_direc + "/ps2_bin_inits.txt");
	
	// Save the travel time between PS1 and PS2b
	save_double_vector(ps2b_time, "Travel time from each bin in PS1 to PS2b.", store_direc + "/ps2_travel_time.txt");
	
	// Save the previous positions from each initial position
	save_array3_vector(ps2b_prev, "Previous point from each initial condition of PS2b.", store_direc + "/ps2_prev_position.txt");
	
	return;
}



// Returns polynomial coefficients for quadratic fit
// 3x3 matrix inverse is specialized
std::vector<double> quadratic_regression(std::vector<double> &ind, std::vector<double> &dep) {
	
	double p1    , p2    , p3    , p4    ; // Single data point powers
	double e1 = 0, e2 = 0, e3 = 0, e4 = 0; // Running sum of data point powers
	double e0 = static_cast<int>(ind.size());
	
	// e_i are entries in the design matrix X^TX
	// Matrix X has rows [1 x_i x_i^2 ... x_i^m] for each data point in ind
	for (unsigned int i = 0; i < ind.size(); i++) {
		// Compute powers of data point
		p1 = ind[i];
		p2 = p1 * p1;
		p3 = p2 * p1;
		p4 = p2 * p2;
		
		// Accumulate powers
		e1 += p1;
		e2 += p2;
		e3 += p3;
		e4 += p4;
	}
	
	// Determinant of matrix X^TX
	double detXTX = e0*e2*e4 + 2*e1*e2*e3 - e0*e3*e3 - e4*e1*e1 - e2*e2*e2;
	
	// Matrix inverse entries
	// Adjoint matrix entries divided by detXTX
	double a00 = (e2*e4 - e3*e3) / detXTX;
	double a01 = (e2*e3 - e1*e4) / detXTX;
	double a02 = (e1*e3 - e2*e2) / detXTX;
	double a11 = (e0*e4 - e2*e2) / detXTX;
	double a12 = (e1*e2 - e0*e3) / detXTX;
	double a22 = (e0*e2 - e1*e1) / detXTX;
	
	// Regression coefficients
	std::vector<double> beta = {0,0,0};
	
	// beta = (X^TX)^-1 X^T dep
	for (unsigned int i = 0; i < ind.size(); i++) {
		beta[0] += (a02 + a12*ind[i] + a22*ind[i]*ind[i]) * dep[i];
		beta[1] += (a01 + a11*ind[i] + a12*ind[i]*ind[i]) * dep[i];
		beta[2] += (a00 + a01*ind[i] + a02*ind[i]*ind[i]) * dep[i];
	}
	
	return beta;
}



// Returns polynomial coefficients for cubic fit
// 4x4 matrix inverse is specialized
std::vector<double> cubic_regression(std::vector<double> &ind, std::vector<double> &dep) {
	
	double p1    , p2    , p3    , p4    , p5    , p6     ; // Single data point powers
	double e1 = 0, e2 = 0, e3 = 0, e4 = 0, e5 = 0, e6 = 0; // Running sum of data point powers
	double e0 = static_cast<int>(ind.size());
	
	// e_i are entries in the design matrix X^TX
	// Matrix X has rows [1 x_i x_i^2 ... x_i^m] for each data point in ind
	for (unsigned int i = 0; i < ind.size(); i++) {
		// Compute powers of data point
		p1 = ind[i];
		p2 = p1 * p1;
		p3 = p2 * p1;
		p4 = p2 * p2;
		p5 = p2 * p3;
		p6 = p3 * p3;
		
		// Accumulate powers
		e1 += p1;
		e2 += p2;
		e3 += p3;
		e4 += p4;
		e5 += p5;
		e6 += p6;
	}
	
	// Determinant of matrix X^TX
	double detXTX = e0*e2*e4 + 2*e1*e2*e3 - e0*e3*e3 - e4*e1*e1 - e2*e2*e2;
	
	// Matrix inverse entries
	// Adjoint matrix entries divided by detXTX
	double a00 = (e2*e4 - e3*e3) / detXTX;
	double a01 = (e2*e3 - e1*e4) / detXTX;
	double a02 = (e1*e3 - e2*e2) / detXTX;
	double a11 = (e0*e4 - e2*e2) / detXTX;
	double a12 = (e1*e2 - e0*e3) / detXTX;
	double a22 = (e0*e2 - e1*e1) / detXTX;
	std::cout << a00 << " " << a01 << " " << a02 << " " << a11 << " " << a12 << " " << a22 << std::endl;
	std::cout << detXTX << std::endl;
	// Regression coefficients
	std::vector<double> beta = {0,0,0};
	
	// beta = (X^TX)^-1 X^T dep
	for (unsigned int i = 0; i < ind.size(); i++) {
		beta[0] += (a02 + a12*ind[i] + a22*ind[i]*ind[i]) * dep[i];
		beta[1] += (a01 + a11*ind[i] + a12*ind[i]*ind[i]) * dep[i];
		beta[2] += (a00 + a01*ind[i] + a02*ind[i]*ind[i]) * dep[i];
	}
	
	return beta;
}
