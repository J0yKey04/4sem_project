#pragma once
#include <string>
#include <vector>

struct RGB {
    double r = 0.0;
    double g = 0.0;
    double b = 0.0;
};

void write_csv(const std::string& path,
               const std::vector<double>& alpha_deg,
               const std::vector<double>& red,
               const std::vector<double>& green,
               const std::vector<double>& blue,
               const std::vector<double>& white);

void write_ppm(const std::string& path,
               int width,
               int height,
               const std::vector<RGB>& pixels);
