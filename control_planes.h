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
void find_surfaces(hindmarsh_rose &neuron, double dt = 1.0/128, double percent = 75.0, const std::string direc = "save_direct");

// Function that reads in x y z values, a peak array, and plane indicator ps.
// Returns a matrix containing the x y z values of either PS1 values (ps = 1) or PS0 values (ps = 0)
std::vector<array3> ps_verts(std::vector<double> &vec_x, std::vector<double> &vec_y, std::vector<double> &vec_z, std::vector<unsigned int> &peaks, int ps);

// Reads in time series t,x,y,z, vertex data for surfaces 0 and 1, as well as neural system neuron.
// Returns list of point of all crossings of PS0 and PS1 as t,x,y,z,PS where PS = 1 for PS1 or 0 for
// PS0. Uses Henon's trick to integrate exactly onto surface.
std::vector<std::vector<double>> ps_points(std::vector<double> &t, std::vector<double> &x, std::vector<double> &y, std::vector<double> &z, std::vector<array3> &verts0, std::vector<array3> &verts1, hindmarsh_rose &neuron);

// Reads in points x_n, x_n1, y_n, y_n1 and using surface ps (indicates which plane, 0 or 1) with
// vertices verts, checks to see if pts are on opposite sides of the surface
bool crossed(double x_n, double x_n1, double y_n, double y_n1, int ps, std::vector<array3> &verts);

#endif