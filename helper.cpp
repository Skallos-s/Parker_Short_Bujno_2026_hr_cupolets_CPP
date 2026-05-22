// File: helper.cpp
// Purpose: Numpy and Scipy methods implementation file
// Author: Daniel Bujno

#include "helper.h"
#include "array3.h"

#include <sstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <vector>
#include <fstream>
#include <iostream>
#include <experimental/filesystem>


// Load data from text file into 4 vectors
void loadtxt_4(std::vector<double> &vec_t, std::vector<double> &vec_x, std::vector<double> &vec_y, std::vector<double> &vec_z, const std::string file_path) {
	
	std::ifstream file;
	std::string line;
	std::istringstream i_line;
	
	double vt, vx, vy, vz;
	
	// Load file
	file.open(file_path);
	
	if (file.is_open()) {
		// Read file
		while (getline(file, line)) {
			if (line[0] == '#') continue;
			
			i_line.clear();
			i_line.str(line);
			
			if (i_line >> vt >> vx >> vy >> vz) {
				vec_t.push_back(vt);
				vec_x.push_back(vx);
				vec_y.push_back(vy);
				vec_z.push_back(vz);
			} else {
				std::cout << "Error Reading Line " << line << std::endl; 
			}
		}
		file.close();
	}
	
	return;
}

// Load data from text file into 3 vectors
void loadtxt_3(std::vector<double> &psx, std::vector<double> &psy, std::vector<double> &psz, const std::string file_path) {
	
	std::ifstream file;
	std::string line;
	std::istringstream i_line;
	
	double vx, vy, vz;
	
	// Load file
	file.open(file_path);
	
	if (file.is_open()) {
		// Read file
		while (getline(file, line)) {
			if (line[0] == '#') continue;
			
			i_line.clear();
			i_line.str(line);
			
			if (i_line >> vx >> vy >> vz) {
				psx.push_back(vx);
				psy.push_back(vy);
				psz.push_back(vz);
			} else {
				std::cout << "Error Reading Line " << line << std::endl; 
			}
		}
		file.close();
	}
	
	return;
}

// Load data from text file into vector of vector<doubles>
void loadtxt_1(std::vector<std::vector<double>> &pts, const std::string file_path) {
	
	std::ifstream file;
	std::string line;
	std::istringstream i_line;
	
	double vt, vx, vy, vz;
	
	// Load file
	file.open(file_path);
	
	if (file.is_open()) {
		// Read file
		while (getline(file, line)) {
			if (line[0] == '#') continue;
			
			i_line.clear();
			i_line.str(line);
			
			std::vector<double> vec;
			
			if (i_line >> vt >> vx >> vy >> vz) {
				vec.push_back(vt);
				vec.push_back(vx);
				vec.push_back(vy);
				vec.push_back(vz);
				pts.push_back(vec);
			} else {
				std::cout << "Error Reading Line " << line << std::endl; 
			}
		}
		file.close();
	}
	
	return;
}

// Load polynomial coordinate data from text file into vector of doubles
void loadtxt_poly(std::vector<double> &poly_coords, const std::string file_path) {
	
	std::ifstream file;
	std::string line;
	std::istringstream i_line;
	
	double pc;
	
	// Load file
	file.open(file_path);
	
	if (file.is_open()) {
		// Read file
		while (getline(file, line)) {
			if (line[0] == '#') continue;
			
			i_line.clear();
			i_line.str(line);
			
			if (i_line >> pc) {
				poly_coords.push_back(pc);
			} else {
				std::cout << "Error Reading Line " << line << std::endl; 
			}
		}
		file.close();
	}
	
	return;
}

// Find local maxima of vector
// Returns vector of indices
std::vector<unsigned int> find_peaks(std::vector<double> &vec) {
	std::vector<unsigned int> vect_of_indices;
	unsigned int run_start, run_end;
	bool in_run = false;
	
	for (unsigned int i = 1; i < vec.size() - 1; i++) {
		// When current value decreases from previous, then run ends
		if (vec[i] < vec[i-1]) in_run = false;
		
		// When current value increases from previous, a new run starts
		if (vec[i] > vec[i-1]) {
			run_start = i;
			in_run = true;
		}
		
		// When next value decreases, a local maximum is found
		if ((vec[i] > vec[i+1]) && in_run) {
			run_end = i;
			
			// Store midpoint
			vect_of_indices.push_back((run_start + run_end)/2);
			in_run = false;
		}
	}
	
	return vect_of_indices;
}

// Like find_peaks, but negates input and finds peaks within [0,1]
std::vector<unsigned int> find_mid_peaks(std::vector<double> &vec) {
	std::vector<unsigned int> vect_of_indices;
	unsigned int run_start, run_end;
	bool in_run = false;
	
	for (unsigned int i = 1; i < vec.size() - 1; i++) {
		// When current value decreases from previous, then run ends
		if (vec[i] > vec[i-1]) in_run = false;
		
		// When current value increases from previous, a new run starts
		if (vec[i] < vec[i-1]) {
			run_start = i;
			in_run = true;
		}
		
		// When next value decreases, a local maximum is found
		if ((vec[i] < vec[i+1]) && in_run) {
			run_end = i;
			
			// Store midpoint if peak is within acceptable range
			if ((0 >= vec[i]) && (vec[i] >= -1)) vect_of_indices.push_back((run_start + run_end)/2);
			in_run = false;
		}
	}
	
	return vect_of_indices;
}

// Write 2D array to file
void save_array3_vector(std::vector<array3> &data, const std::string header, const std::string location) {
	std::ofstream file;
	file.open(location);
	file << "# " << header << std::endl;
	file << std::scientific;
	file << std::setprecision(18); // double precision
	for (unsigned int i = 0; i < data.size(); i++) {
		file << data[i].get(0) << "\t" << data[i].get(1) << "\t" << data[i].get(2) << std::endl;
	}
	file.close();
}

// Write 2D array to file
// data contains vector<double> of length 5 (t x y z ps)
void save_vector_vector(std::vector<std::vector<double>> &data, const std::string header, const std::string location) {
	std::ofstream file;
	file.open(location);
	file << "# " << header << std::endl;
	file << std::scientific;
	file << std::setprecision(18); // double precision
	for (unsigned int i = 0; i < data.size(); i++) {
		for (int j = 0; j < 4; j++) {
			file << data[i][j] << "\t";
		}
		file << data[i][4] << std::endl;
	}
	file.close();
}

// Write 1d array to file
void save_double_vector(std::vector<double> &poly_coords, const std::string header, const std::string location) {
	std::ofstream file;
	file.open(location);
	file << "# " << header << std::endl;
	file << std::scientific;
	file << std::setprecision(18); // double precision
	for (unsigned int i = 0; i < poly_coords.size(); i++) {
		file << poly_coords[i] << std::endl;
	}
	file.close();
}

// Write vector of binary sequence to file
void save_binary(std::vector<std::vector<unsigned int>> &map, const std::string header, const std::string location) {
	std::ofstream file;
	file.open(location);
	file << "# " << header << std::endl;
	for (unsigned int i = 0; i < map.size(); i++) {
		for (unsigned int j = 0; j < map[i].size(); j++) {
			file << map[i][j];
		}
		file << std::endl;
	}
	file.close();
}

// Evaluate a polynomial with ax^n+...+b where poly_coords=<a,...,b>
double eval_poly(std::vector<double> &poly_coords, double x) {
	double y = 0;
	
	for (unsigned int i = 0; i < poly_coords.size(); i++) {
		y *= x;
		y += poly_coords[i];
	}
	
	return y;
}