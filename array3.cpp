// File: array.cpp
// Purpose: Implementation file for the array class
// Author: Daniel Bujno

#include "array3.h"
#include <stdexcept>
#include <iostream>
#include <iomanip>

// Constructor
array3::array3(double xp, double yp, double zp) {
	x = xp;
	y = yp;
	z = zp;
}

// Getter
double array3::get(int n) const {
	if (n == 0) return x;
	if (n == 1) return y;
	if (n == 2) return z;
	throw std::runtime_error("invalid index");
}

// Setter
void array3::set(int n, double kp) {
	     if (n == 0) x = kp;
	else if (n == 1) y = kp;
	else if (n == 2) z = kp;
	else throw std::runtime_error("invalid index");
}

// Unzip
void array3::unzip(double &vx, double &vy, double &vz) {
	vx = x;
	vy = y;
	vz = z;
}


// Operator overloading

// Negation of array
array3 operator-(const array3 &a1) {
	return array3(-a1.x, -a1.y, -a1.z);
}

// Addition of two arrays
array3 operator+(const array3 &a1, const array3 &a2) {
	double vx = a1.x + a2.x;
	double vy = a1.y + a2.y;
	double vz = a1.z + a2.z;
	return array3(vx, vy, vz);
}

// Subtraction of two arrays
array3 operator-(const array3 &a1, const array3 &a2) {
	double vx = a1.x - a2.x;
	double vy = a1.y - a2.y;
	double vz = a1.z - a2.z;
	return array3(vx, vy, vz);
}

// Scaling of array by double
array3 operator*(const array3 &a1, const double f) {
	double vx = a1.x * f;
	double vy = a1.y * f;
	double vz = a1.z * f;
	return array3(vx, vy, vz);
}

// Scaling of array by double
array3 operator*(const double f, const array3 &a2) {
	double vx = a2.x * f;
	double vy = a2.y * f;
	double vz = a2.z * f;
	return array3(vx, vy, vz);
}

// Scaling of array by double
array3 operator/(const array3 &a1, const double f) {
	double vx = a1.get(0) / f;
	double vy = a1.get(1) / f;
	double vz = a1.get(2) / f;
	return array3(vx, vy, vz);
}

// Addition assignment of array and array
void operator+=(array3 &a1, const array3 &a2) {
	a1.x += a2.x;
	a1.y += a2.y;
	a1.z += a2.z;
}

// Subtraction assignment of array and array
void operator-=(array3 &a1, const array3 &a2) {
	a1.x -= a2.x;
	a1.y -= a2.y;
	a1.z -= a2.z;
}

// Multiplication assignment of array and double
void operator*=(array3 &a1, const double f) {
	a1.x *= f;
	a1.y *= f;
	a1.z *= f;
}
// Division assignment of array and double
void operator/=(array3 &a1, const double f) {
	a1.x /= f;
	a1.y /= f;
	a1.z /= f;
}


// Stream output
std::ostream& operator<<(std::ostream &out, const array3 &a1) {
	return out << std::setprecision(18) << "[" << a1.x << ", " << a1.y << ", " << a1.z << "]";
}
