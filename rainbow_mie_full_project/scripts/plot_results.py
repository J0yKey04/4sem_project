#!/usr/bin/env python3
import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "results"
CSV = RESULTS / "angular_intensity.csv"

alpha = []
red = []
green = []
blue = []
white = []

with CSV.open() as f:
    reader = csv.DictReader(f)
    for row in reader:
        alpha.append(float(row["alpha_deg"]))
        red.append(float(row["red_650nm"]))
        green.append(float(row["green_550nm"]))
        blue.append(float(row["blue_450nm"]))
        white.append(float(row["white_sum"]))

def log_norm(values):
    return [math.log10(1.0 + x) for x in values]

plt.figure(figsize=(10, 6))
plt.plot(alpha, log_norm(red), label="650 nm red")
plt.plot(alpha, log_norm(green), label="550 nm green")
plt.plot(alpha, log_norm(blue), label="450 nm blue")
plt.xlabel("alpha, degrees from antisolar point")
plt.ylabel("log10(1 + Mie intensity), arbitrary units")
plt.title("Mie angular scattering: rainbow region")
plt.xlim(30, 60)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(RESULTS / "rainbow_region_plot.png", dpi=180)

plt.figure(figsize=(10, 6))
plt.plot(alpha, log_norm(white), label="white = RGB sum")
plt.xlabel("alpha, degrees from antisolar point")
plt.ylabel("log10(1 + intensity), arbitrary units")
plt.title("Full Mie phase function, backscattering hemisphere")
plt.xlim(0, 90)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(RESULTS / "full_backscatter_plot.png", dpi=180)

print("Written:")
print(RESULTS / "rainbow_region_plot.png")
print(RESULTS / "full_backscatter_plot.png")
