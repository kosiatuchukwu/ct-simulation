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

E_tab, mu_rho = [], []
with open("water_mu.txt") as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 2:
            try:
                e, m = float(parts[0]), float(parts[1])
            except ValueError:
                continue
            E_tab.append(e * 1000.0); mu_rho.append(m)
E_tab, mu_tab = np.array(E_tab), np.array(mu_rho) * 1.0

def mu_water(E):
    return np.exp(np.interp(np.log(E), np.log(E_tab), np.log(mu_tab)))

try:
    import spekpy as sp
    s = sp.Spek(kvp=KVP, th=12); s.filter('Al', FILTER_MM_AL)
    energies, fluence = s.get_spectrum()
    print(f"Spectrum: SpekPy, W anode {KVP} kVp + {FILTER_MM_AL} mm Al")
except ImportError:
    energies = np.arange(16.0, KVP + 0.5, 1.0)
    fluence  = (energies - 15.0)**1.5 * (KVP - energies)**2
    print("Spectrum: SYNTHETIC fallback")
keep = fluence > 1e-6 * fluence.max()
energies, fluence = energies[keep], fluence[keep]

mean_E_counting    = np.sum(energies * fluence) / np.sum(fluence)
mean_E_integrating = np.sum(energies**2 * fluence) / np.sum(energies * fluence)
print(f"Mean energy seen by counting detector:    {mean_E_counting:.1f} keV")
print(f"Mean energy seen by integrating detector: {mean_E_integrating:.1f} keV")

yy, xx = np.mgrid[:GRID, :GRID]
c = (GRID - 1) / 2
mask = (((xx - c)**2 + (yy - c)**2) <= (CIRC_PIX / 2)**2).astype(float)
theta = np.linspace(0.0, 180.0, N_ANGLES, endpoint=False)
path_cm = radon(mask, theta=theta) * PIXEL_SIZE_CM

sig  = {"counting": np.zeros_like(path_cm), "integrating": np.zeros_like(path_cm)}
sig0 = {"counting": 0.0, "integrating": 0.0}
for E, N0 in zip(energies, fluence):
    T = np.exp(-mu_water(E) * path_cm)
    sig["counting"]     += N0 * T
    sig["integrating"]  += E * N0 * T
    sig0["counting"]    += N0
    sig0["integrating"] += E * N0

r = np.sqrt((xx - c)**2 + (yy - c)**2)
centre_zone = r < 0.15 * (CIRC_PIX / 2)
edge_zone   = (r > 0.80 * (CIRC_PIX / 2)) & (r < 0.92 * (CIRC_PIX / 2))

recon, cup = {}, {}
for mode in ("counting", "integrating"):
    sino = -np.log(sig[mode] / sig0[mode])
    recon[mode] = iradon(sino, theta=theta, filter_name=FILTER_NAME) / PIXEL_SIZE_CM
    mu_c = recon[mode][centre_zone].mean()
    mu_e = recon[mode][edge_zone].mean()
    cup[mode] = (mu_e - mu_c) / mu_e
    print(f"{mode:12s}: mu_centre={mu_c:.4f}  mu_edge={mu_e:.4f}  cupping={cup[mode]*100:.2f}%")

# ---- standalone step-5 image: energy-integrating reconstruction --------
plt.figure(figsize=(5.5, 5))
plt.imshow(recon["integrating"], cmap="gray")
plt.title(f"Energy-integrating reconstruction, {KVP} kVp")
plt.axis("off")
plt.colorbar(fraction=0.046)
plt.savefig(f"image_integrating_{KVP}kvp.png", dpi=120, bbox_inches="tight")
print(f"Saved image_integrating_{KVP}kvp.png")

# ---- four-panel overview ------------------------------------------------
mid = GRID // 2
fig, axes = plt.subplots(2, 2, figsize=(11.5, 9.5))
im0 = axes[0, 0].imshow(recon["counting"], cmap="gray")
axes[0, 0].set_title(f"Image: photon-counting ({KVP} kVp)")
axes[0, 0].axis("off"); fig.colorbar(im0, ax=axes[0, 0], fraction=0.046)
im1 = axes[0, 1].imshow(recon["integrating"], cmap="gray")
axes[0, 1].set_title(f"Image: energy-integrating ({KVP} kVp)")
axes[0, 1].axis("off"); fig.colorbar(im1, ax=axes[0, 1], fraction=0.046)
axes[1, 0].plot(recon["counting"][mid, :],    color="#c0392b",
                label=f"photon-counting ({cup['counting']*100:.2f}%)")
axes[1, 0].plot(recon["integrating"][mid, :], color="#2c3e50",
                label=f"energy-integrating ({cup['integrating']*100:.2f}%)")
axes[1, 0].set_title("Centre profiles")
axes[1, 0].set_xlabel("pixel"); axes[1, 0].set_ylabel("mu (1/cm)")
axes[1, 0].legend(); axes[1, 0].grid(alpha=0.3)
axes[1, 1].plot(energies, fluence / fluence.max(), color="#c0392b",
                label="what counting sees: N0(E)")
w = energies * fluence
axes[1, 1].plot(energies, w / w.max(), color="#2c3e50",
                label="what integrating sees: E*N0(E)")
axes[1, 1].set_title("Effective spectrum per detector model")
axes[1, 1].set_xlabel("Energy (keV)"); axes[1, 1].set_ylabel("relative weight")
axes[1, 1].legend(); axes[1, 1].grid(alpha=0.3)
plt.tight_layout()
plt.savefig(f"detector_models_{KVP}kvp.png", dpi=120, bbox_inches="tight")
print(f"Saved detector_models_{KVP}kvp.png")
