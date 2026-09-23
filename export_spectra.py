"""Export filtered SpekPy spectra to CSV for the MATLAB engine."""
import numpy as np
import spekpy as sp

for kvp in (50, 80, 120):
    s = sp.Spek(kvp=kvp, th=12)
    s.filter('Al', 2.5)
    energies, fluence = s.get_spectrum()
    keep = fluence > 1e-6 * fluence.max()
    data = np.column_stack([energies[keep], fluence[keep]])
    fname = f"spectrum_{kvp}kvp.csv"
    np.savetxt(fname, data, delimiter=",",
               header="energy_keV,fluence", comments="")
    print(f"Saved {fname}  ({data.shape[0]} energy bins)")
