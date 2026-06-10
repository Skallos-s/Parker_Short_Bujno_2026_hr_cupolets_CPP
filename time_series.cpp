// File: time_series.cpp
// Purpose: Generate time series data of Hindmarsh Rose Attractor
// Author: Daniel Bujno


#include "hindmarsh_rose.h"
#include "time_series.h"
#include "keep_data.h"
#include "helper.h"
#include "rk4.h"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>

// Generate time series data of hr with step size dt on the interval [0,t] with initial conditions IC.
// Save to direc/hr_system/hr_time_series.txt
void generate_time_series(hindmarsh_rose hr, array3 IC, const std::string direc, double dt, double t_max) {
	
	std::cout << "Generating time series data..." << std::endl;
	
	// Directory to store time series data
	const std::string store_direc = direc + "/hr_system";
	
	// Create hr system directory if it does not exist
	check_direc(store_direc);
	
	// Initial conditions
	array3 curr = IC;
	
	// Number of iterations
	unsigned int N = static_cast<unsigned int>(std::round(t_max / dt));
	
	// Time series data
	std::vector<std::vector<double>> time_series;
	
	// Generate time series data
	for (unsigned int i = 0; i < N; i++) {
		// Save current position
		std::vector<double> row{i * dt, curr.get(0), curr.get(1), curr.get(2)};
		time_series.push_back(row);
		
		// Integrate one step forward
		array3 next = rk4(curr, dt, &hindmarsh_rose::hr_dynamics, hr);
		
		// Set new values to old values
		curr = next;
	}
	
	// Save time series data to file
	save_vector_vector(time_series, "t x y z data from HR integration", store_direc + "/hr_time_series.txt");
	
	return;
}

// Generate time series data of hr with step size dt for N iterations with initial conditions IC.
// Save to direc/hr_system/hr_time_series.txt
void generate_time_series(hindmarsh_rose hr, array3 IC, const std::string direc, double dt, unsigned int N) {
	
	std::cout << "Generating time series data..." << std::endl;
	
	// Directory to store time series data
	const std::string store_direc = direc + "/hr_system";
	
	// Create hr system directory if it does not exist
	check_direc(store_direc);
	
	// Initial conditions
	array3 curr = IC;
	
	// Time series data
	std::vector<std::vector<double>> time_series;
	
	// Generate time series data
	for (unsigned int i = 0; i < N; i++) {
		// Save current position
		std::vector<double> row{i * dt, curr.get(0), curr.get(1), curr.get(2)};
		time_series.push_back(row);
		
		// Integrate one step forward
		array3 next = rk4(curr, dt, &hindmarsh_rose::hr_dynamics, hr);
		
		// Set new values to old values
		curr = next;
	}
	
	// Save time series data to file
	save_vector_vector(time_series, "t x y z data from HR integration", store_direc + "/hr_time_series.txt");
	
	return;
}