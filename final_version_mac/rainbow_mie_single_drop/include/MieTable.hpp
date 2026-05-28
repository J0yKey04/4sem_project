#pragma once

#include <string>
#include <vector>

class MieTable {
public:
    bool loadCSV(const std::string& filename);

    double lookupNearest(
        double radius_um,
        double wavelength_nm,
        double theta_deg
    ) const;

    bool isLoaded() const {
        return loaded;
    }

private:
    struct Sample {
        double radius_um;
        double wavelength_nm;
        double theta_deg;
        double intensity;
    };

    std::vector<Sample> samples;
    bool loaded = false;
};



