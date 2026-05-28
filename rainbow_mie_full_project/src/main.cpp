#include "mie.hpp"
#include "dispersion.hpp"
#include "image.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;
constexpr double PI = 3.141592653589793238462643383279502884;

struct ChannelResult {
    double wavelength_nm;
    std::vector<MiePoint> phase;
};

static double get_intensity_by_alpha(const std::vector<MiePoint>& phase, double alpha_deg) {
    if (alpha_deg <= 0.0) return phase.back().i_unpolarized;
    if (alpha_deg >= 180.0) return phase.front().i_unpolarized;

    // phase is stored by theta increasing, alpha = 180 - theta decreasing.
    const double theta_deg = 180.0 - alpha_deg;
    const double pos = theta_deg / 180.0 * static_cast<double>(phase.size() - 1);
    const int i0 = static_cast<int>(std::floor(pos));
    const int i1 = std::min(i0 + 1, static_cast<int>(phase.size() - 1));
    const double t = pos - i0;
    return (1.0 - t) * phase[i0].i_unpolarized + t * phase[i1].i_unpolarized;
}

static std::vector<double> log_normalize(const std::vector<double>& v) {
    std::vector<double> out(v.size());
    double maxv = 0.0;
    for (double x : v) maxv = std::max(maxv, std::log10(1.0 + x));
    if (maxv <= 0.0) return out;
    for (size_t i = 0; i < v.size(); ++i) out[i] = std::log10(1.0 + v[i]) / maxv;
    return out;
}

int main(int argc, char** argv) {
    try {
        double radius_um = 50.0;
        int angle_samples = 3601;
        int width = 1200;
        int height = 1200;

        if (argc >= 2) radius_um = std::stod(argv[1]);
        if (argc >= 3) angle_samples = std::stoi(argv[2]);
        if (argc >= 4) width = std::stoi(argv[3]);
        if (argc >= 5) height = std::stoi(argv[4]);

        fs::create_directories("results");

        const std::vector<double> wavelengths = {650.0, 550.0, 450.0};
        std::vector<ChannelResult> channels;

        std::cout << "Mie rainbow simulation\n";
        std::cout << "radius_um = " << radius_um << "\n";
        std::cout << "angle_samples = " << angle_samples << "\n";

        for (double wl : wavelengths) {
            MieParameters p;
            p.radius_um = radius_um;
            p.wavelength_nm = wl;
            p.refractive_index = {water_refractive_index_visible(wl), 0.0};
            p.angle_samples = angle_samples;
            std::cout << "Computing lambda = " << wl << " nm, n = " << p.refractive_index.real() << "\n";
            channels.push_back({wl, compute_mie_phase_function(p)});
        }

        std::vector<double> alpha(angle_samples);
        std::vector<double> red(angle_samples), green(angle_samples), blue(angle_samples), white(angle_samples);

        // Reorder by alpha increasing from 0 to 180 degrees.
        for (int i = 0; i < angle_samples; ++i) {
            const double a = 180.0 * static_cast<double>(i) / static_cast<double>(angle_samples - 1);
            alpha[i] = a;
            red[i] = get_intensity_by_alpha(channels[0].phase, a);
            green[i] = get_intensity_by_alpha(channels[1].phase, a);
            blue[i] = get_intensity_by_alpha(channels[2].phase, a);
            white[i] = red[i] + green[i] + blue[i];
        }

        write_csv("results/angular_intensity.csv", alpha, red, green, blue, white);

        const auto rnorm = log_normalize(red);
        const auto gnorm = log_normalize(green);
        const auto bnorm = log_normalize(blue);

        std::vector<RGB> pixels(static_cast<size_t>(width * height));
        const double cx = 0.5 * (width - 1);
        const double cy = 0.5 * (height - 1);
        const double max_radius_px = 0.48 * std::min(width, height);

        // Full 360-degree azimuthal rainbow around the antisolar point.
        // Image radius maps alpha from 0 to 90 degrees: enough to show primary ~42 deg and secondary ~51 deg.
        const double alpha_max_deg = 90.0;

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                const double dx = x - cx;
                const double dy = y - cy;
                const double rho = std::sqrt(dx * dx + dy * dy);
                const double alpha_deg = alpha_max_deg * rho / max_radius_px;

                RGB color;
                if (alpha_deg <= alpha_max_deg) {
                    const double pos = alpha_deg / 180.0 * static_cast<double>(angle_samples - 1);
                    const int i0 = static_cast<int>(std::floor(pos));
                    const int i1 = std::min(i0 + 1, angle_samples - 1);
                    const double t = pos - i0;
                    color.r = (1.0 - t) * rnorm[i0] + t * rnorm[i1];
                    color.g = (1.0 - t) * gnorm[i0] + t * gnorm[i1];
                    color.b = (1.0 - t) * bnorm[i0] + t * bnorm[i1];

                    // Dark blue background so weak scattering is visible.
                    color.r = 0.015 + 0.95 * color.r;
                    color.g = 0.020 + 0.95 * color.g;
                    color.b = 0.035 + 0.95 * color.b;
                }
                pixels[static_cast<size_t>(y * width + x)] = color;
            }
        }

        write_ppm("results/full_360_rainbow.ppm", width, height, pixels);

        std::cout << "Done. Files written:\n";
        std::cout << "  results/angular_intensity.csv\n";
        std::cout << "  results/full_360_rainbow.ppm\n";
        std::cout << "\nInterpretation:\n";
        std::cout << "  alpha = angular distance from the antisolar point.\n";
        std::cout << "  Primary rainbow is near alpha ~ 40-43 deg.\n";
        std::cout << "  Secondary rainbow is near alpha ~ 50-54 deg, weaker and color-reversed.\n";
        std::cout << "  Mie oscillations are physical interference/supernumerary structure.\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
