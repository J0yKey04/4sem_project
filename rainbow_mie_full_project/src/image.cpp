#include "image.hpp"
#include <algorithm>
#include <fstream>
#include <stdexcept>

namespace {
int to_byte(double v) {
    v = std::max(0.0, std::min(1.0, v));
    return static_cast<int>(255.0 * v + 0.5);
}
}

void write_csv(const std::string& path,
               const std::vector<double>& alpha_deg,
               const std::vector<double>& red,
               const std::vector<double>& green,
               const std::vector<double>& blue,
               const std::vector<double>& white) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot open CSV file: " + path);
    f << "alpha_deg,red_650nm,green_550nm,blue_450nm,white_sum\n";
    for (size_t i = 0; i < alpha_deg.size(); ++i) {
        f << alpha_deg[i] << ',' << red[i] << ',' << green[i] << ',' << blue[i] << ',' << white[i] << '\n';
    }
}

void write_ppm(const std::string& path,
               int width,
               int height,
               const std::vector<RGB>& pixels) {
    if (static_cast<int>(pixels.size()) != width * height) {
        throw std::runtime_error("Pixel buffer size does not match image dimensions");
    }
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot open PPM file: " + path);
    f << "P6\n" << width << ' ' << height << "\n255\n";
    for (const RGB& p : pixels) {
        unsigned char rgb[3] = {
            static_cast<unsigned char>(to_byte(p.r)),
            static_cast<unsigned char>(to_byte(p.g)),
            static_cast<unsigned char>(to_byte(p.b))
        };
        f.write(reinterpret_cast<const char*>(rgb), 3);
    }
}
