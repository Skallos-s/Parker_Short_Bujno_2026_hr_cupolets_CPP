// File: control_planes.h
// Purpose: Control planes header file
// Author: Daniel Bujno

#ifndef __control_planes_h_
#define __control_planes_h_

#include <string>
#include <vector>

#include "hindmarsh_rose.h"
#include "array3.h"

// Shorthand namespace
// namespace vecfl = std::vector<double>;
// namespace veca3 = std::vector<array3>;

// Reads in neuron, dt, and directory direc.
// Generates all related Poincare surface information
// through various functions and plots that are stored in direc.
void find_surfaces(hindmarsh_rose &neuron, bool GENERATE_CONTROL_PLANE = true, bool GENERATE_PS2_CONTROL_PLANE = true, double percent = 75.0, const std::string direc = "save_direct");

// Function that reads in x y z values, a peak array, and plane indicator ps.
// Returns a matrix containing the x y z values of either PS1 values (ps = 1) or PS0 values (ps = 0)
std::vector<array3> ps_verts(std::vector<double> &vec_x, std::vector<double> &vec_y, std::vector<double> &vec_z, std::vector<unsigned int> &peaks, int ps);

// Like ps_verts, but PS1 is split into PS2a and PS2b. The planes are shifted horizontally along y. Moreover the input is assumed to be ps1.
// Returns a matrix containing the x y z values of PS2a, followed by he x y z values of PS2b
// The y positions of the plane have average yavg (calculated) and distance y_width (given).
std::vector<array3> ps2_verts(std::vector<double> &vec_x, std::vector<double> &vec_y, std::vector<double> &vec_z, double y_width);

// Reads in time series t,x,y,z, vertex data for surfaces 0 and 1, as well as neural system neuron.
// Returns list of point of all crossings of PS0 and PS1 as t,x,y,z,PS where PS = 1 for PS1 or 0 for
// PS0. Uses Henon's trick to integrate exactly onto surface.
// If verts1 has 4 entries, assumed to be split. 1 is PS1a, 2 is PS1b.
std::vector<std::vector<double>> ps_points(std::vector<double> &t, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z, std::vector<array3> &verts0, std::vector<array3> &verts1, hindmarsh_rose &neuron);

// Reads in points x_n, x_n1, y_n, y_n1 and using surface ps (indicates which plane, 0 for PS0, 1 for PS1 and PS1a, and 2 for PS1b) with
// vertices verts, checks to see if pts are on opposite sides of the surface
bool crossed(double x_n, double x_n1, double y_n, double y_n1, int ps, std::vector<array3> &verts);

#endif