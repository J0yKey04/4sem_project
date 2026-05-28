#include "dispersion.hpp"
#include <cmath>

// Empirical visible-light approximation for liquid water near room temperature.
// lambda is in micrometers. Accuracy is sufficient for educational rainbow modelling.
double water_refractive_index_visible(double wavelength_nm) {
    const double lambda_um = wavelength_nm * 1e-3;
    const double l2 = lambda_um * lambda_um;

    // Cauchy-type fit around visible wavelengths.
    // Gives approximately n(400 nm) ~= 1.344, n(550 nm) ~= 1.335, n(700 nm) ~= 1.331.
    const double A = 1.32292;
    const double B = 0.003145;
    const double C = 0.000273;
    return A + B / l2 + C / (l2 * l2);
}
