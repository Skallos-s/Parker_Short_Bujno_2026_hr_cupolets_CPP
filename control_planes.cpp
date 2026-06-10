// File: control_plane.cpp
// Purpose: Implementation file for Control Plane methods
// Author: Daniel Bujno

#include "hindmarsh_rose.h"
#include "control_planes.h"
#include "keep_data.h"
#include "helper.h"
#include "array3.h"
#include "rk4.h"

#include <cmath>
#include <cfloat>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <experimental/filesystem>

#include <iomanip>

// Reads in neuron, dt, and directory direc.
// Generates all related Poincare surface information
// through various functions and plots that are stored in direc.
void find_surfaces(hindmarsh_rose &neuron, double percent, const std::string direc) {
	
	std::cout << "Generating Poincare surfaces and corresponding intersections..." << std::endl;
	
	// Define directory name for storing control plane data
	const std::string store_direc = direc + "/control_planes";
	
	// Create control plane directory if it does not exist
	check_direc(store_direc);
	
	// Read in data
	std::vector<double> t, x, y, z;
	loadtxt_4(t, x, y, z, direc + "/hr_system/hr_time_series.txt");
	
	// Find index for percent data remaining
	int start = static_cast<int>(round((1 - percent / 100) * (t.size())));
	
	// Redefine data for percent data
	if (start > 0) {
		t.erase(t.begin(), t.begin() + start);
		x.erase(x.begin(), x.begin() + start);
		y.erase(y.begin(), y.begin() + start);
		z.erase(z.begin(), z.begin() + start);
	}
	
	// Find indicies of peaks for x time series - used for PS1
	std::vector<unsigned int> peaks = find_peaks(x);
	std::vector<unsigned int> mid_peaks = find_mid_peaks(x);
	
	// Find the vertices of the surface for PS0 and PS1
	std::vector<array3> verts0 = ps_verts(x,y,z,mid_peaks,0);
	std::vector<array3> verts1 = ps_verts(x,y,z,peaks,1);
	
	
	// Find all the points on the control planes
	std::vector<std::vector<double>> pts = ps_points(t, x, y, z, verts0, verts1, neuron);
	
	// Points on PS0 and PS1
	std::vector<std::vector<double>> pts0;
	std::vector<std::vector<double>> pts1;
	
	// y values of PS0 points on plane
	std::vector<double> yp0;
	
	
	// These lines are used to shrink PS0 to a better fit for the points on the plane
	for (unsigned int i = 0; i < pts.size(); i++) {
		if (pts[i][4] == 0) {
			// Grab all the y values for the PS0 points on plane
			yp0.push_back(pts[i][2]);
			
			// PS0 points
			pts0.push_back(pts[i]);
		} else {
			// PS1 points
			pts1.push_back(pts[i]);
		}
	}
	
	// Initialize minimum and maximum double values
	double ypmin = DBL_MAX, ypmax = -DBL_MAX;
	
	// Find the min and max of the y values
	for (unsigned int i = 0; i < yp0.size(); i++) {
		ypmin = std::min(ypmin, yp0[i]);
		ypmax = std::max(ypmax, yp0[i]);
	}
    
	// Slightly extend the y range by 50%
    double ypinc = abs_diff(ypmax , ypmin) * 0.50;
	
	// Replace y values of PS0 with augmented ones
    verts0[0].set(1, ypmin - ypinc);
	verts0[1].set(1, ypmax + ypinc);
	
	// Save PS0 and PS1 data
	save_array3_vector(verts0, "X Y Z limits of PS0 plane", store_direc + "/ps0_vertices.txt");
	save_array3_vector(verts1, "X Y Z limits of PS1 plane", store_direc + "/ps1_vertices.txt");
	
	// Save all the PS0 and PS1 points
	save_vector_vector(pts0, "PS0 points on plane."        , store_direc + "/ps0_pts.txt");
	save_vector_vector(pts1, "PS1 points on plane."        , store_direc + "/ps1_pts.txt");
	save_vector_vector(pts, "All points on control planes.", store_direc + "/ps_all_pts.txt");
	
	std::cout << "Completed integration of Poincare surfaces and established control planes. Saved to " + direc + "/control_planes" << std::endl;
	
	return;
}

// Function that reads in x y z values, a peak array, and plane indicator ps.
// Returns a matrix containing the x y z values of either PS1 values (ps = 1) or PS0 values (ps = 0)
std::vector<array3> ps_verts(std::vector<double> &vec_x, std::vector<double> &vec_y, std::vector<double> &vec_z, std::vector<unsigned int> &peaks, int ps) {
	
	std::vector<array3> verts;
	
	// Initialize minimum and maximum double values
	double xpmin = DBL_MAX, xpmax = -DBL_MAX;
	double ypmin = DBL_MAX, ypmax = -DBL_MAX;
	double zpmin = DBL_MAX, zpmax = -DBL_MAX;
	double yavg = 0;
	
	// Min and max of x peaks
	for (unsigned int i = 0; i < peaks.size(); i++) {
		xpmin = std::min(xpmin, vec_x[peaks[i]]);
		xpmax = std::max(xpmax, vec_x[peaks[i]]);
	}
	
	// Min and max of y values
	for (unsigned int i = 0; i < vec_y.size(); i++) {
		ypmin = std::min(ypmin, vec_y[i]);
		ypmax = std::max(ypmax, vec_y[i]);
		yavg += (vec_y[i] - yavg) / (i+1); // Iterative (accurate) average
	}
	
	// Min and max of z values
	for (unsigned int i = 0; i < vec_z.size(); i++) {
		zpmin = std::min(zpmin, vec_z[i]);
		zpmax = std::max(zpmax, vec_z[i]);
	}
	
	// Increase x y z endpoints
	double x_ext = 0.05 * abs_diff(xpmax , xpmin);
	double y_ext = 0.05 * abs_diff(ypmax , ypmin);
	double z_ext = 0.01 * abs_diff(zpmax , zpmin);
	
	if (ps == 0) {
		// Establish PS0 values
		verts.push_back(array3(xpmin, ypmin-y_ext, zpmin-z_ext));
		verts.push_back(array3(xpmin, ypmax+y_ext, zpmax+z_ext));
	} else { // ps == 1
		// Establish PS1 values
		verts.push_back(array3(xpmin-x_ext, yavg, zpmin-z_ext));
		verts.push_back(array3(xpmax+x_ext, yavg, zpmax+z_ext));
	}
	return verts;
}

// Reads in time series t,x,y,z, vertex data for surfaces 0 and 1, as well as neural system neuron.
// Returns list of point of all crossings of PS0 and PS1 as t,x,y,z,PS where PS = 1 for PS1 or 0 for
// PS0. Uses Henon's trick to integrate exactly onto surface.
std::vector<std::vector<double>> ps_points(std::vector<double> &t, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z, std::vector<array3> &verts0, std::vector<array3> &verts1, hindmarsh_rose &neuron) {
	
	// Array to be filled with pts on contro plane
	std::vector<std::vector<double>> pts;
	
	// Iterate through all points
	for (unsigned int i = 0; i + 1 < x.size(); i++) {
		double xn = x[i], xn1 = x[i+1]; // Grab consecutive x points
		double yn = y[i], yn1 = y[i+1]; // Grab consecutive y points
		double zn = z[i];               // Grab z point
		
		// Check to see PS1 has been crossed
		if (crossed(xn, xn1, yn, yn1, 1, verts1)) {
			// If crossed, use Henon's trick to integrate onto surface
			std::vector<double> txzyp = rk4_henon(array3(t[i], xn, zn), yn, -(yn-verts1[0].get(1)), &hindmarsh_rose::hr_dy_dynamics, neuron);
			std::vector<double> output = {txzyp[0], txzyp[1], txzyp[3], txzyp[2], 1};
			pts.push_back(output);
		}
		
		// Check to see PS0 has been crossed
		if (crossed(xn, xn1, yn, yn1, 0, verts0)) {
			// If crossed, use Henon's trick to integrate onto surface
			std::vector<double> tyzxp = rk4_henon(array3(t[i], yn, zn), xn, -(xn-verts0[0].get(0)), &hindmarsh_rose::hr_dx_dynamics, neuron);
			std::vector<double> output = {tyzxp[0], tyzxp[3], tyzxp[1], tyzxp[2], 0};
			pts.push_back(output);
		}
	}
	
	return pts;
}

// Reads in points x_n, x_n1, y_n, y_n1 and using surface ps (indicates which plane, 0 or 1) with
// vertices verts, checks to see if pts are on opposite sides of the surface
bool crossed(double x_n, double x_n1, double y_n, double y_n1, int ps, std::vector<array3> &verts) {
	// Check for crossing PS1
	if (ps == 1) {
		// Change of signs in y direction indicates crossing PS1
		double sp = y_n - verts[0].get(1), sp1 = y_n1 - verts[0].get(1);
		if (x_n <= verts[1].get(0) && x_n1 >= verts[0].get(0) && sgn(sp1) != sgn(sp) && sgn(sp) == 1) {
			// PS1 crossed if yn and yn1 signs differ, and the xn and xn1 are within the PS1 boundary
			return true;
		}
	}
	
	// Check for crossing PS0
	if (ps == 0) {
		// Change of signs in x direction indicates crossing PS0
		double sp = x_n - verts[0].get(0), sp1 = x_n1 - verts[0].get(0);
		if (y_n <= verts[1].get(1) && y_n1 >= verts[0].get(1) && sgn(sp1) != sgn(sp) && sgn(sp) == 1) {
			// PS0 crossed if xn and xn1 signs differ, and the yn and yn1 are within the PS0 boundary
			return true;
		}
	}
	
	// False if neither plane crossed
	return false;
}