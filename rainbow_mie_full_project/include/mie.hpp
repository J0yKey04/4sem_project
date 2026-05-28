#pragma once
#include <complex>
#include <vector>

struct MiePoint {
    double theta_deg;        // scattering angle: 0 forward, 180 backward
    double alpha_deg;        // rainbow angle from antisolar point: alpha = 180 - theta
    double i_parallel;       // |S2|^2
    double i_perpendicular;  // |S1|^2
    double i_unpolarized;    // 0.5 * (|S1|^2 + |S2|^2)
};

struct MieParameters {
    double radius_um = 50.0;
    double wavelength_nm = 550.0;
    std::complex<double> refractive_index = {1.333, 0.0};
    int angle_samples = 1801;
};

std::vector<MiePoint> compute_mie_phase_function(const MieParameters& p);
