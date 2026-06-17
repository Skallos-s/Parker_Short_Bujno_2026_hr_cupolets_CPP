// File: array.h
// Purpose: 3D array header file
// Author: Daniel Bujno

#ifndef __array3_h_
#define __array3_h_

#include <ostream>

class array3 {
public:
	// Constructor
	array3();
	array3(double xp, double yp, double zp);
	
	// Getter
	double get(int n) const;
	
	// Setter
	void set(int n, double kp);
	
	// Unzip
	void unzip(double &vx, double &vy, double &vz);
	
	// Operator overloading
	
	// Arithmetic operator friends
	friend array3 operator-(const array3 &a1);
	friend array3 operator+(const array3 &a1, const array3 &a2);
	friend array3 operator-(const array3 &a1, const array3 &a2);
	friend array3 operator*(const array3 &a1, const double f);
	friend array3 operator*(const double f, const array3 &a2);
	friend array3 operator/(const array3 &a1, const double f);
	friend void operator+=(array3 &a1, const array3 &a2);
	friend void operator-=(array3 &a1, const array3 &a2);
	friend void operator*=(array3 &a1, const double f);
	friend void operator/=(array3 &a1, const double f);
	
	// Stream operator friend
	friend std::ostream& operator<<(std::ostream &out, const array3 &a1);
	
private:
	double x;
	double y;
	double z;
};

// Arithmetic operators
array3 operator-(const array3 &a1);
array3 operator+(const array3 &a1, const array3 &a2);
array3 operator-(const array3 &a1, const array3 &a2);
array3 operator*(const array3 &a1, const double f);
array3 operator*(const double f, const array3 &a2);
array3 operator/(const array3 &a1, const double f);
void operator+=(array3 &a1, const array3 &a2);
void operator-=(array3 &a1, const array3 &a2);
void operator*=(array3 &a1, const double f);
void operator/=(array3 &a1, const double f);

// Stream operators
std::ostream& operator<<(std::ostream &out, const array3 &a1);

#endif