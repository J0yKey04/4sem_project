import os
import numpy as np
import pandas as pd
from tqdm import tqdm
import miepython


def water_refractive_index(wavelength_nm: float) -> float:
    lam_um = wavelength_nm / 1000.0
    return 1.322 + 0.003 / (lam_um * lam_um)


def compute_s1_s2(m: complex, x: float, mu: np.ndarray):
    if hasattr(miepython, "S1_S2"):
        return miepython.S1_S2(m, x, mu)

    if hasattr(miepython, "mie_S1_S2"):
        return miepython.mie_S1_S2(m, x, mu)

    raise RuntimeError("Cannot find S1_S2 function in miepython.")


def main():
    os.makedirs("data", exist_ok=True)

    radii_um = np.array([5.0, 10.0, 20.0, 50.0, 100.0])
    wavelengths_nm = np.arange(400.0, 701.0, 10.0)
    angles_deg = np.arange(0.0, 181.0, 1.0)

    rows = []

    for radius_um in tqdm(radii_um, desc="Generating Mie table"):
        radius_nm = radius_um * 1000.0

        for wavelength_nm in wavelengths_nm:
            n_water = water_refractive_index(wavelength_nm)
            m = complex(n_water, 0.0)

            x = 2.0 * np.pi * radius_nm / wavelength_nm
            mu = np.cos(np.deg2rad(angles_deg))

            s1, s2 = compute_s1_s2(m, x, mu)

            intensity = 0.5 * (np.abs(s1) ** 2 + np.abs(s2) ** 2)

            max_i = np.max(intensity)
            if max_i > 0.0:
                intensity = intensity / max_i

            for theta_deg, value in zip(angles_deg, intensity):
                rows.append(
                    {
                        "radius_um": radius_um,
                        "wavelength_nm": wavelength_nm,
                        "theta_deg": theta_deg,
                        "intensity": float(value),
                    }
                )

    df = pd.DataFrame(rows)
    df.to_csv("data/mie_table.csv", index=False)

    print("Saved: data/mie_table.csv")
    print(f"Rows: {len(df)}")


if __name__ == "__main__":
    main()


    