// File: hindmarsh_rose.h
// Purpose: Hindmarsh Rose header file
// Author: Daniel Bujno

#ifndef __hindmarsh_rose_h_
#define __hindmarsh_rose_h_

#include "array3.h"

class hindmarsh_rose {
public:
	// Default Constructor
	hindmarsh_rose();
	
	// Parameter Constructor
	hindmarsh_rose(double Ca, double Cb, double Cc,
				   double Cd, double Cs, double Cxr,
				   double Cr, double CI);
	
	// Copy
	//hindmarsh_rose(const hindmarsh_rose &h);
	
	// Destructor
	// NOTE: Only needed when class has pointer data
	// ~hindmarsh_rose();
	
	// Set Parameters
	void set_a (double Ca);
	void set_b (double Cb);
	void set_c (double Cc);
	void set_d (double Cd);
	void set_s (double Cs);
	void set_xr(double Cxr);
	void set_r (double Cr);
	void set_I (double CI);
	
	// Get Parameters
	double get_a ();
	double get_b ();
	double get_c ();
	double get_d ();
	double get_s ();
	double get_xr();
	double get_r ();
	double get_I ();
	
	// Dynamics
	array3 hr_dynamics(array3 state) const;
	array3 hr_dy_dynamics(array3 state, double y) const;
	array3 hr_dx_dynamics(array3 state, double x) const;
	//array3 hr_jacobian(double t, array3 state);
	
private:
	// Parameters
	double a; // Default parameter value in Hindmarsh & Rose 1984
	double b; // ~~~
	double c; // ~~~
	double d; // ~~~
	double s; // ~~~
	double xr; // ~~
	double r; // Bifurcation parameter, default chosen to be in chaotic regime with default I
	double I; // Represents incoming current, default chosen to be in chaotic regime with default r
};

#endif