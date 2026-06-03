// File: macrocontrol.h
// Purpose: macrocontrol header file
// Author: Daniel Bujno

#ifndef __macrocontrol_h_
#define __macrocontrol_h_

#include <string>
#include <vector>

// Reads in the top level directory name. Generates the macrocontrol sequence that corresponds
// to the number of bins and crossings and saves in directory /direc/macrocontrol.
void establish_macrocontrol(const std::string bin_rn_direc);

// Reads in a coding sequence rn. Write to an array ctrl that maps bin B to the nearest
// smallest difference of the coding sequence and an array diff that says the size of that difference.
void macro_control_map(std::vector<unsigned int> &ctrl, std::vector<double> &diff, std::vector<double> rn);

#endif