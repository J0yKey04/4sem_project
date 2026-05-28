#include "MieTable.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <limits>
#include <cmath>

bool MieTable::loadCSV(const std::string& filename) {
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Warning: cannot open Mie table: " << filename << "\n";
        loaded = false;
        return false;
    }

    samples.clear();

    std::string line;
    std::getline(file, line);

    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }

        std::stringstream ss(line);
        std::string cell;

        Sample s{};

        if (!std::getline(ss, cell, ',')) continue;
        s.radius_um = std::stod(cell);

        if (!std::getline(ss, cell, ',')) continue;
        s.wavelength_nm = std::stod(cell);

        if (!std::getline(ss, cell, ',')) continue;
        s.theta_deg = std::stod(cell);

        if (!std::getline(ss, cell, ',')) continue;
        s.intensity = std::stod(cell);

        samples.push_back(s);
    }

    loaded = !samples.empty();

    std::cout << "Loaded Mie samples: " << samples.size() << "\n";

    return loaded;
}

double MieTable::lookupNearest(
    double radius_um,
    double wavelength_nm,
    double theta_deg
) const {
    if (!loaded || samples.empty()) {
        return 1.0;
    }

    double bestDistance = std::numeric_limits<double>::max();
    double bestIntensity = 1.0;

    for (const auto& s : samples) {
        double dr = (s.radius_um - radius_um) / 100.0;
        double dw = (s.wavelength_nm - wavelength_nm) / 300.0;
        double dt = (s.theta_deg - theta_deg) / 180.0;

        double distance = dr * dr + dw * dw + dt * dt;

        if (distance < bestDistance) {
            bestDistance = distance;
            bestIntensity = s.intensity;
        }
    }

    return bestIntensity;
}


