// File: keep_data.cpp
// Purpose: Implementation file for Keep Data methods
// Author: Daniel Bujno

#include "keep_data.h"
#include <string>
#include <experimental/filesystem>

namespace fs = std::experimental::filesystem;

// Create directory if it doesn't exist.
void check_direc(const std::string file_path) {
	fs::create_directories(file_path);
	return;
}