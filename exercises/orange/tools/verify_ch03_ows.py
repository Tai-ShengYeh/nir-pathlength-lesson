"""Load ch03_wavelength_pitfall.ows headless, let the signals propagate, and print what the
students should see. Orange's numbers are compared with a numpy replica of the same DIRECT
regression (c = a + b*d2) and with the lesson's INVERSE calibration (d2 = a + b*c).

    python exercises/orange/tools/verify_ch03_ows.py [path/to/ch03_wavelength_pitfall.ows] [--img DIR]
"""
import os, sys, argparse
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import owsutil, owsbuild
import numpy as np
from scipy.signal import savgol_filter

ap = argparse.ArgumentParser()
ap.add_argument("ows", nargs="?", default=os.path.join(here, "..", "ch03_wavelength_pitfall.ows"))
ap.add_argument("--img", default=None)
args = ap.parse_args()
model = owsbuild.load_ows(args.ows)
W = lambda t: owsutil.by_title(model, t)[1]

wf = W("File")
print("File path loaded:", wf.recent_paths[0].abspath if wf.recent_paths else None,
      "| rows:", None if wf.data is None else len(wf.data))
pp = W("Preprocess Spectra (d2, 88 wl)")
rank = W("Rank (outside CV)")
d = rank.data
names = [a.name for a in d.domain.attributes]
print("Preprocess list:", pp.storedsettings)
print("after Preprocess Spectra: rows", len(d), "| conc (file order):", list(d.get_column("acid_pct")),
      "| n wavelengths:", len(names), "|", names[0], "...", names[-1])

# exact mask check: the 88 names must be the wavelengths where all four raw 1 mm spectra have A < 1.5
raw = W("Select Rows (1 mm)")
import Orange
sel_raw = raw.data[(raw.data.get_column("path_mm") == 1) & (raw.data.get_column("replicate") == 1)]
usable = (sel_raw.X < 1.5).all(axis=0)
mask_names = [a.name for a, u in zip(sel_raw.domain.attributes, usable) if u]
print("88-wavelength mask identical to A<1.5 mask:", names == mask_names, "(", len(mask_names), ")")


def rank_top(w, k=4):
    m = w.ranksModel
    hdr = [m.headerData(c, 1) for c in range(m.columnCount())]
    rows = []
    for r in range(min(k, m.rowCount())):
        rows.append([m.data(m.index(r, c)) for c in range(m.columnCount())])
    return hdr, rows


for t in ["Rank (outside CV)", "Rank (shuffled target)"]:
    w = W(t)
    hdr, rows = rank_top(w)
    out = w.Outputs.reduced_data
    print(t, "| methods:", sorted(w.selected_methods), "| sorting:", w.sorting, "| nSelected:", w.nSelected)
    print("   header:", hdr)
    for r in rows:
        F = float(r[2]) if r[2] is not None else float("nan")
        print("   row:", r[0].name if hasattr(r[0], "name") else r[0], "F = %.3f" % F, "-> R2 = F/(F+2) = %.4f" % (F / (F + 2)))

tsA = W("Test and Score (A: outside)")
print("Test and Score (A) input features:", [a.name for a in tsA.data.domain.attributes], "| resampling", tsA.resampling)
print("   row:", owsbuild.score_row(model, "Test and Score (A: outside)"))
resA = list(tsA.learners.values())[0].results.value
print("   LOO (actual, predicted, error):", [(round(float(a), 2), round(float(p), 3), round(float(p - a), 3))
                                           for a, p in zip(resA.actual, resA.predicted[0])])
tsB = W("Test and Score (B: inside)")
print("Test and Score (B) input features:", len(tsB.data.domain.attributes), "| preprocessor:", tsB.preprocessor,
      "| resampling", tsB.resampling)
print("   row:", owsbuild.score_row(model, "Test and Score (B: inside)"))
resB = list(tsB.learners.values())[0].results.value
print("   LOO (actual, predicted, error):", [(round(float(a), 2), round(float(p), 3), round(float(p - a), 3))
                                           for a, p in zip(resB.actual, resB.predicted[0])])

rz = W("Randomize (shuffle target)")
dS = W("Rank (shuffled target)").data
print("Randomize: replicable =", rz.random_seed, "| scope", rz.scope_prop, "% | shuffled acid_pct (file order):",
      list(dS.get_column("acid_pct")), "| original:", list(d.get_column("acid_pct")))
print("messages:", owsutil.messages(model))

# Which wavelength does each fold pick? (Orange API, the same preprocessor the widget sends)
from Orange.preprocess.fss import SelectBestFeatures
from Orange.preprocess.score import UnivariateLinearRegression
picks = []
for k in range(len(d)):
    tr = d[[i for i in range(len(d)) if i != k]]
    picks.append((float(d[k]["acid_pct"]), SelectBestFeatures(UnivariateLinearRegression(), k=1)(tr).domain.attributes[0].name))
print("inside-CV picks (left-out conc, wavelength):", picks)

# ---- numpy replica -------------------------------------------------------------------
import pandas as pd
csv = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(args.ows)), "../../data/spectra_wide_orange.csv"))
df = pd.read_csv(csv); wc = df.columns[4:]
s = df[(df["mC#path_mm"] == 1) & (df["mC#replicate"] == 1)]
X = s[wc].values; y = s["cC#acid_pct"].values
D = savgol_filter(X, 9, 2, deriv=2, axis=1, mode="nearest"); us = (X < 1.5).all(axis=0)
r2 = lambda a, b: float(np.corrcoef(a, b)[0, 1] ** 2)
R = np.where(us, [r2(D[:, j], y) for j in range(len(wc))], -1); ju = int(np.argmax(R))
eo_d, ei_d, eo_i, ei_i = [], [], [], []
for k in range(4):
    tr = [i for i in range(4) if i != k]
    jr = int(np.argmax(np.where(us, [r2(D[tr, j], y[tr]) for j in range(len(wc))], -1)))
    for j, ed, ei in [(ju, eo_d, eo_i), (jr, ei_d, ei_i)]:
        b, a = np.polyfit(D[tr, j], y[tr], 1); ed.append(a + b * D[k, j] - y[k])
        b, a = np.polyfit(y[tr], D[tr, j], 1); ei.append((D[k, j] - a) / b - y[k])
rm = lambda e: float(np.sqrt(np.mean(np.square(e))))
print("\nnumpy replica: best usable %s nm R2 %.4f" % (wc[ju], R[ju]))
print("               | Orange LOO RMSE | numpy direct | numpy inverse (lesson)")
print("A outside CV   | %15.3f | %12.3f | %22.3f" % (float(owsbuild.score_row(model, "Test and Score (A: outside)")["RMSE"]), rm(eo_d), rm(eo_i)))
print("B inside CV    | %15.3f | %12.3f | %22.3f" % (float(owsbuild.score_row(model, "Test and Score (B: inside)")["RMSE"]), rm(ei_d), rm(ei_i)))
yS = dS.get_column("acid_pct")
RS = np.where(us, [r2(D[:, j], yS) for j in range(len(wc))], -1)
print("shuffled target: best usable %s nm R2 %.4f | fixed 2258.067 nm R2 %.4f" % (
    wc[int(np.argmax(RS))], RS.max(), r2(D[:, list(wc).index("2258.067")], yS)))

if args.img:
    os.makedirs(args.img, exist_ok=True)
    owsbuild.shot(model, "Preprocess Spectra (d2, 88 wl)", os.path.join(args.img, "ch03_orange_preprocess.png"), (1100, 900))
    owsbuild.shot(model, "Rank (outside CV)", os.path.join(args.img, "ch03_orange_rank.png"), (900, 560))
    owsbuild.shot(model, "Preprocess (select 1 inside CV)", os.path.join(args.img, "ch03_orange_fss.png"), (900, 560))
    owsbuild.shot(model, "Test and Score (A: outside)", os.path.join(args.img, "ch03_orange_testscore_outside.png"), (1000, 560))
    owsbuild.shot(model, "Test and Score (B: inside)", os.path.join(args.img, "ch03_orange_testscore_inside.png"), (1000, 560))
    owsbuild.shot(model, "Rank (shuffled target)", os.path.join(args.img, "ch03_orange_rank_shuffled.png"), (900, 560))
