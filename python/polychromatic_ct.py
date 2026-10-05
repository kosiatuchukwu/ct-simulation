import numpy as np
import matplotlib.pyplot as plt
from skimage.transform import radon, iradon

KVP          = 120
FILTER_MM_AL = 2.5
OBJ_DIAM_CM  = 10.0
GRID         = 180
CIRC_PIX     = 160
N_ANGLES     = 180
FILTER_NAME  = "hamming"

PIXEL_SIZE_CM = OBJ_DIAM_CM / CIRC_PIX

E_tab, mu_rho_tab = [], []
with open("water_mu.txt") as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 2:
            try:
                e, m = float(parts[0]), float(parts[1])
            except ValueError:
                continue
            E_tab.append(e * 1000.0)
            mu_rho_tab.append(m)
E_tab  = np.array(E_tab)
mu_tab = np.array(mu_rho_tab) * 1.0
print(f"Loaded {len(E_tab)} XCOM points, {E_tab.min():.0f}-{E_tab.max():.0f} keV")

def mu_water(E_keV):
    return np.exp(np.interp(np.log(E_keV), np.log(E_tab), np.log(mu_tab)))

try:
    import spekpy as sp
    s = sp.Spek(kvp=KVP, th=12)
    s.filter('Al', FILTER_MM_AL)
    energies, fluence = s.get_spectrum()
    print(f"Spectrum: SpekPy, W anode {KVP} kVp + {FILTER_MM_AL} mm Al")
except ImportError:
    energies = np.arange(16.0, KVP + 0.5, 1.0)
    fluence  = (energies - 15.0)**1.5 * (KVP - energies)**2
    print("Spectrum: SYNTHETIC fallback (install spekpy for the real one)")

keep     = fluence > 1e-6 * fluence.max()
energies = energies[keep]
fluence  = fluence[keep]
mean_E   = np.sum(energies * fluence) / np.sum(fluence)
print(f"{len(energies)} energy bins, mean energy = {mean_E:.1f} keV")

yy, xx = np.mgrid[:GRID, :GRID]
c      = (GRID - 1) / 2
mask   = (((xx - c)**2 + (yy - c)**2) <= (CIRC_PIX / 2)**2).astype(float)

theta   = np.linspace(0.0, 180.0, N_ANGLES, endpoint=False)
path_cm = radon(mask, theta=theta) * PIXEL_SIZE_CM

N0_tot = np.zeros_like(path_cm)
N_tot  = np.zeros_like(path_cm)
for E, N0_E in zip(energies, fluence):
    sino_E  = mu_water(E) * path_cm
    N_tot  += N0_E * np.exp(-sino_E)
    N0_tot += N0_E

sino_poly = -np.log(N_tot / N0_tot)
sino_mono = mu_water(mean_E) * path_cm

recon_poly = iradon(sino_poly, theta=theta, filter_name=FILTER_NAME) / PIXEL_SIZE_CM
recon_mono = iradon(sino_mono, theta=theta, filter_name=FILTER_NAME) / PIXEL_SIZE_CM

r = np.sqrt((xx - c)**2 + (yy - c)**2)
centre_zone = r < 0.15 * (CIRC_PIX / 2)
edge_zone   = (r > 0.80 * (CIRC_PIX / 2)) & (r < 0.92 * (CIRC_PIX / 2))

def cupping(recon):
    mu_c = recon[centre_zone].mean()
    mu_e = recon[edge_zone].mean()
    return mu_c, mu_e, (mu_e - mu_c) / mu_e

for name, rec in [("poly", recon_poly), ("mono", recon_mono)]:
    mu_c, mu_e, C = cupping(rec)
    print(f"{name}:  mu_centre={mu_c:.4f}  mu_edge={mu_e:.4f}  cupping={C*100:.2f}%")

mid = GRID // 2
fig, axes = plt.subplots(2, 2, figsize=(11, 9))
im0 = axes[0,0].imshow(recon_poly, cmap="gray")
axes[0,0].set_title(f"Polychromatic ({KVP} kVp spectrum)")
axes[0,0].axis("off"); fig.colorbar(im0, ax=axes[0,0], fraction=0.046)
im1 = axes[0,1].imshow(recon_mono, cmap="gray")
axes[0,1].set_title(f"Monochromatic ({mean_E:.0f} keV)")
axes[0,1].axis("off"); fig.colorbar(im1, ax=axes[0,1], fraction=0.046)
axes[1,0].plot(recon_poly[mid, :], label="polychromatic", color="#c0392b")
axes[1,0].plot(recon_mono[mid, :], label="monochromatic", color="#2c3e50")
axes[1,0].set_title("Profile through centre")
axes[1,0].set_xlabel("pixel"); axes[1,0].set_ylabel("mu (1/cm)")
axes[1,0].legend(); axes[1,0].grid(alpha=0.3)
axes[1,1].plot(energies, fluence, color="#2c3e50")
axes[1,1].set_title("Spectrum used N0(E)")
axes[1,1].set_xlabel("Energy (keV)"); axes[1,1].set_ylabel("Fluence")
axes[1,1].grid(alpha=0.3)
plt.tight_layout()
plt.savefig("polychromatic_ct.png", dpi=120, bbox_inches="tight")
plt.show()
print("Saved polychromatic_ct.png")
