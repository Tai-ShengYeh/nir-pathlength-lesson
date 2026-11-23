"""Write exercises/orange/ch09_validation.csv: the chapter-9 validation table (deterministic, no random seed).

    python exercises/orange/tools/make_ch09_validation.py

One row per sample (water + 3 bottles).  Same numbers as tools/nir-pathlength/make_data_ch09.py
(-> site/nir-pathlength/lib/ch09data.js, CH9):
  d2_2258 / d2_1670   1 mm, replicate 1, Savitzky-Golay (9 pt, quadratic) 2nd derivative x 1e3 (data/spectra_long.csv)
  label_pct           bottle label 0 / 3 / 5 / 10 %
  atago_1102          ATAGO PAL-Easy ACID 181 (conductivity meter) reading of 11/02 (= acid_pct in spectra_long.csv)
  titr_pct            titration acidity = V * K (V = 0 / 8.8 / 17 / 44 mL NaOH, K = f*0.006*(100/25)*(100/10), f = 1.07875)
  pal30s_pct          PAL-30S refractometer (measured 11/17, shows acidity % directly)
Header uses Orange's flag#name format: c = class (d2_2258), m = meta, C = numeric, S = string."""
import os, csv
import numpy as np, pandas as pd
from scipy.signal import savgol_filter

here = os.path.dirname(os.path.abspath(__file__))
orange_dir = os.path.normpath(os.path.join(here, ".."))
root = os.path.normpath(os.path.join(here, "..", "..", ".."))

K = 0.4315 / 40 / 0.1 * 10 * 0.006 * (100 / 25) * (100 / 10)    # % acidity per mL NaOH (= 0.2589)
bottles = [("water", 0.0, 0.0, 0.0, 0.0), ("vinegar_3%", 3.0, 3.32, 8.8, 2.8),
           ("vinegar_5%", 5.0, 4.87, 17.0, 4.8), ("vinegar_10%", 10.0, 7.45, 44.0, 10.4)]   # name, label, atago 11/02, V (mL), PAL-30S

d = pd.read_csv(os.path.join(root, "data", "spectra_long.csv"))
d = d[(d.replicate == 1) & (d.path_mm == 1.0)]
g = d.pivot_table(index="acid_pct", columns="wavelength_nm", values="absorbance").sort_index()
wl = g.columns.values
D = savgol_filter(g.values, 9, 2, deriv=2, axis=1) * 1e3
j58, j70 = int(np.argmin(abs(wl - 2258))), int(np.argmin(abs(wl - 1670)))
assert np.allclose(g.index.values, [b[2] for b in bottles])      # rows of the pivot = Atago 11/02 values, ascending

out = os.path.join(orange_dir, "ch09_validation.csv")
with open(out, "w", encoding="utf-8-sig", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["mS#bottle", "cC#d2_2258", "mC#d2_1670", "mC#label_pct", "mC#atago_1102", "mC#titr_pct", "mC#pal30s_pct"])
    for i, (name, lab, a1, v, pal) in enumerate(bottles):
        w.writerow([name, "%.6f" % D[i, j58], "%.6f" % D[i, j70], lab, a1, "%.6f" % (v * K), pal])
print("wrote", os.path.basename(out), len(bottles), "rows; K = %.6f" % K)
