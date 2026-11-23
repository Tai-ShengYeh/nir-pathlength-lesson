"""Write exercises/orange/ch08_references.csv: the chapter-8 reference-method table (deterministic, no random seed).

    python exercises/orange/tools/make_ch08_references.py

One row per sample (water + 3 bottles).  Same numbers as tools/nir-pathlength/make_data_ch08.py
(-> site/nir-pathlength/lib/ch08data.js, CH8):
  titr_mL    NaOH volume of the titration (0 for water); K = f * 0.006 * (100/25) * (100/10), f = 1.07875
  titr_pct   = titr_mL * K  (acidity, w/v %)
  d2_2258    1 mm, replicate 1, Savitzky-Golay (9 pt, quadratic) 2nd derivative x 1e3 at 2258 nm (data/spectra_long.csv)
  pal30s_pct PAL-30S refractometer reading of 11/17 (acidity %, 0 for water): 2.8 / 4.8 / 10.4 (ATAGO PAL-30S, not conductivity)
  titr_pct_f097 = titr_pct * 0.97  (illustrative -3 % titre-factor correction, NOT a measured value)
Header uses Orange's flag#name format: c = class (d2_2258), m = meta, C = numeric, S = string."""
import os, csv
import numpy as np, pandas as pd
from scipy.signal import savgol_filter

here = os.path.dirname(os.path.abspath(__file__))
orange_dir = os.path.normpath(os.path.join(here, ".."))
root = os.path.normpath(os.path.join(here, "..", "..", ".."))

F = 0.4315 / 40 / 0.1 * 10                      # 1.07875 (NaOH factor, not standardised)
K = F * 0.006 * (100 / 25) * (100 / 10)         # % acidity per mL NaOH (= 0.2589)
bottles = [("water", 0.0, 0.0, 0.0, 0.0, 0.0), ("vinegar_3%", 3.0, 3.32, 3.01, 2.8, 8.8),
           ("vinegar_5%", 5.0, 4.87, 4.12, 4.8, 17.0), ("vinegar_10%", 10.0, 7.45, 6.86, 10.4, 44.0)]

d = pd.read_csv(os.path.join(root, "data", "spectra_long.csv"))
d = d[(d.replicate == 1) & (d.path_mm == 1.0)]
g = d.pivot_table(index="sample", columns="wavelength_nm", values="absorbance")
wl = g.columns.values; j = int(np.argmin(abs(wl - 2258)))
D = savgol_filter(g.values, 9, 2, deriv=2, axis=1) * 1e3
sig = {s: float(D[i, j]) for i, s in enumerate(g.index)}

out = os.path.join(orange_dir, "ch08_references.csv")
with open(out, "w", encoding="utf-8-sig", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["mS#bottle", "mC#label_pct", "mC#atago_1102", "mC#atago_1103", "mC#pal30s_pct", "mC#titr_mL", "mC#titr_pct",
                "cC#d2_2258", "mC#titr_pct_f097"])
    for name, lab, a1, a2, pal, v in bottles:
        t = v * K
        w.writerow([name, lab, a1, a2, pal, v, "%.6f" % t, "%.6f" % sig[name], "%.6f" % (t * 0.97)])
print("wrote", os.path.basename(out), len(bottles), "rows; K = %.6f" % K)
