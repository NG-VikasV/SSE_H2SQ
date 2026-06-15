#ifndef READ_INPUT_HPP
#define READ_INPUT_HPP

#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <unordered_map>

#include "variables.hpp"

// --------------------------------------------------
// Read Parameters from input file
// --------------------------------------------------
inline void read_parameters(const std::string& filename,
                            Parameters& prm)
{
    std::ifstream fin(filename);
    if (!fin) {
        std::cerr << "ERROR: Cannot open input file: "
                  << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }

    std::unordered_map<std::string, double> values;
    std::string line;

    while (std::getline(fin, line)) {
        // Skip comments and empty lines
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        std::string key;
        double val;

        if (!(iss >> key >> val)) continue;
        values[key] = val;
    }

    // Mandatory parameters
    prm.Lx = int(values.at("Lx"));
    prm.Ly = int(values.at("Ly"));
    prm.Lz = int(values.at("Lz"));

    prm.beta = values.at("beta");

    prm.J0 = values.at("J0");
    prm.J1 = values.at("J1");
    prm.J2 = values.at("J2");
    prm.J3 = values.at("J3");
    prm.hx = values.at("hx");

    prm.n_therm    = int(values.at("n_therm"));
    prm.n_measure = int(values.at("n_measure"));

    // Derived quantities
    prm.Ns = prm.Lx * prm.Ly * prm.Lz;
    prm.T  = 1.0 / prm.beta;

    // Plaquettes (for square lattice, now 3D stack)
    // Note: This is total number of plaquettes in the system
    prm.Np = prm.Lx * prm.Ly * prm.Lz / 2;
}

#endif
