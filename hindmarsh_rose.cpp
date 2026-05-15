// File: hindmarsh_rose.cpp
// Purpose: Implementation file for the hindmarsh_rose class
// Author: Daniel Bujno

#include "array3.h"
#include "hindmarsh_rose.h"

// Default Constructor
// Initialize Hindmarsh-Rose neuron with specified parameters
hindmarsh_rose::hindmarsh_rose() {
	a  = 1;
	b  = 3;
	c  = 1;
	d  = 5;
	s  = 4;
	xr = -8.0/5;
	r  = 0.006;
	I  = 3.25;
}
	
// Parameter Constructor
// Initialize Hindmarsh-Rose neuron with specified parameters
hindmarsh_rose::hindmarsh_rose(double Ca, double Cb, double Cc,
				               double Cd, double Cs, double Cxr,
				               double Cr, double CI) {
	a  = Ca;
	b  = Cb;
	c  = Cc;
	d  = Cd;
	s  = Cs;
	xr = Cxr;
	r  = Cr;
	I  = CI;
}

// Set Parameters
void hindmarsh_rose::set_a(double Ca)   {a  = Ca;}
void hindmarsh_rose::set_b(double Cb)   {b  = Cb;}
void hindmarsh_rose::set_c(double Cc)   {c  = Cc;}
void hindmarsh_rose::set_d(double Cd)   {d  = Cd;}
void hindmarsh_rose::set_s(double Cs)   {s  = Cs;}
void hindmarsh_rose::set_xr(double Cxr) {xr = Cxr;}
void hindmarsh_rose::set_r(double Cr)   {r  = Cr;}
void hindmarsh_rose::set_I(double CI)   {I  = CI;}

// Get Parameters
double hindmarsh_rose::get_a()  {return a;}
double hindmarsh_rose::get_b()  {return b;}
double hindmarsh_rose::get_c()  {return c;}
double hindmarsh_rose::get_d()  {return d;}
double hindmarsh_rose::get_s()  {return s;}
double hindmarsh_rose::get_xr() {return xr;}
double hindmarsh_rose::get_r()  {return r;}
double hindmarsh_rose::get_I()  {return I;}


// Generates the Hindmarsh-Rose dynamical state based off state (vector of x y z values)
array3 hindmarsh_rose::hr_dynamics(array3 state) const {
	double x,y,z;
	state.unzip(x, y, z); // Redefine values from state for readability
	
	double nx,ny,nz;
	nx = y - a*x*x*x + b*x*x - z + I;
	ny = c - d*x*x - y;
	nz = r * (s * (x - xr) - z);
	return array3(nx, ny, nz);
}

// Generates the Hindmarsh-Rose dynamical state based off state (vector of x y z values)
array3 hindmarsh_rose::hr_dy_dynamics(array3 state, double y) const {
	double t,x,z;
	state.unzip(t, x, z); // Redefine values from state for readability
	
	double dy,nt,nx,nz;
	dy = c - d*x*x - y;
	nt = 1 / dy;
	nx = (y - a*x*x*x + b*x*x - z + I) / dy;
	nz = r * (s * (x - xr) - z) / dy;
	return array3(nt, nx, nz);
}

// Generates the Hindmarsh-Rose dynamical state based off state (vector of x y z values)
array3 hindmarsh_rose::hr_dx_dynamics(array3 state, double x) const {
	double t,y,z;
	state.unzip(t, y, z); // Redefine values from state for readability
	
	double dx,nt,ny,nz;
	dx = y - a*x*x*x + b*x*x - z + I;
	nt = 1 / dx;
	ny = (c - d*x*x - y) / dx;
	nz = r * (s * (x - xr) - z) / dx;
	return array3(nt, ny, nz);
}