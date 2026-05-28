#include "mie.hpp"
#include <cmath>
#include <stdexcept>

namespace {
constexpr double PI = 3.141592653589793238462643383279502884;
using cd = std::complex<double>;

int estimate_series_terms(double x) {
    return static_cast<int>(std::ceil(x + 4.0 * std::cbrt(x) + 2.0));
}

std::vector<cd> logarithmic_derivative_D(cd z, int nmax) {
    // Downward recurrence for D_n(z) = psi'_n(z) / psi_n(z).
    // Start sufficiently far above nmax for stability.
    const int nstop = nmax + static_cast<int>(std::ceil(15.0 + std::abs(z))) + 50;
    std::vector<cd> D(nstop + 2, cd(0.0, 0.0));
    for (int n = nstop; n >= 1; --n) {
        D[n - 1] = cd(static_cast<double>(n), 0.0) / z - cd(1.0, 0.0) / (D[n] + cd(static_cast<double>(n), 0.0) / z);
    }
    D.resize(nmax + 2);
    return D;
}

void compute_mie_coefficients(double x, cd m, std::vector<cd>& a, std::vector<cd>& b) {
    if (x <= 0.0) throw std::runtime_error("Size parameter x must be positive");

    const int nmax = estimate_series_terms(x);
    a.assign(nmax + 1, cd(0.0, 0.0));
    b.assign(nmax + 1, cd(0.0, 0.0));

    const cd mx = m * x;
    const auto D = logarithmic_derivative_D(mx, nmax + 1);

    // Riccati-Bessel functions psi_n(x), chi_n(x).
    // psi_0 = sin x; psi_1 = sin x / x - cos x.
    // chi_0 = cos x; chi_1 = cos x / x + sin x.
    double psi_nm1 = std::sin(x);
    double psi_n = std::sin(x) / x - std::cos(x);
    double chi_nm1 = std::cos(x);
    double chi_n = std::cos(x) / x + std::sin(x);

    for (int n = 1; n <= nmax; ++n) {
        const cd xi_n(psi_n, -chi_n);
        const cd xi_nm1(psi_nm1, -chi_nm1);
        const cd Dn = D[n];
        const double nn = static_cast<double>(n);

        const cd numerator_a = (Dn / m + nn / x) * psi_n - psi_nm1;
        const cd denominator_a = (Dn / m + nn / x) * xi_n - xi_nm1;

        const cd numerator_b = (m * Dn + nn / x) * psi_n - psi_nm1;
        const cd denominator_b = (m * Dn + nn / x) * xi_n - xi_nm1;

        a[n] = numerator_a / denominator_a;
        b[n] = numerator_b / denominator_b;

        const double psi_np1 = ((2.0 * nn + 1.0) / x) * psi_n - psi_nm1;
        const double chi_np1 = ((2.0 * nn + 1.0) / x) * chi_n - chi_nm1;
        psi_nm1 = psi_n;
        psi_n = psi_np1;
        chi_nm1 = chi_n;
        chi_n = chi_np1;
    }
}

} // namespace

std::vector<MiePoint> compute_mie_phase_function(const MieParameters& p) {
    if (p.angle_samples < 2) throw std::runtime_error("angle_samples must be >= 2");

    const double wavelength_um = p.wavelength_nm * 1e-3;
    const double x = 2.0 * PI * p.radius_um / wavelength_um;

    std::vector<cd> a, b;
    compute_mie_coefficients(x, p.refractive_index, a, b);
    const int nmax = static_cast<int>(a.size()) - 1;

    std::vector<MiePoint> out;
    out.reserve(static_cast<size_t>(p.angle_samples));

    for (int i = 0; i < p.angle_samples; ++i) {
        const double theta = PI * static_cast<double>(i) / static_cast<double>(p.angle_samples - 1);
        const double mu = std::cos(theta);

        double pi_nm1 = 0.0; // pi_0
        double pi_n = 1.0;   // pi_1
        cd S1(0.0, 0.0);
        cd S2(0.0, 0.0);

        for (int n = 1; n <= nmax; ++n) {
            const double nn = static_cast<double>(n);
            const double tau_n = nn * mu * pi_n - (nn + 1.0) * pi_nm1;
            const double factor = (2.0 * nn + 1.0) / (nn * (nn + 1.0));

            S1 += factor * (a[n] * pi_n + b[n] * tau_n);
            S2 += factor * (a[n] * tau_n + b[n] * pi_n);

            const double pi_np1 = ((2.0 * nn + 1.0) / (nn + 1.0)) * mu * pi_n
                                - (nn / (nn + 1.0)) * pi_nm1;
            pi_nm1 = pi_n;
            pi_n = pi_np1;
        }

        MiePoint mp;
        mp.theta_deg = theta * 180.0 / PI;
        mp.alpha_deg = 180.0 - mp.theta_deg;
        mp.i_perpendicular = std::norm(S1);
        mp.i_parallel = std::norm(S2);
        mp.i_unpolarized = 0.5 * (mp.i_perpendicular + mp.i_parallel);
        out.push_back(mp);
    }
    return out;
}
