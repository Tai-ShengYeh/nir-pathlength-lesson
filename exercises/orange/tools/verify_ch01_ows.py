"""Load ch01_pathlength.ows headless, let the signals propagate, and print what the students
should see; then switch path (Select Rows) and wavelength (Select Columns) on the live widgets
and read Test and Score again. Orange's numbers are compared with a numpy replica of the same
DIRECT regression (c = a + b*d2) and with the lesson's INVERSE calibration (d2 = a + b*c).

    python exercises/orange/tools/verify_ch01_ows.py [path/to/ch01_pathlength.ows] [--img DIR]
"""
import os, sys, argparse
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import owsutil, owsbuild
import numpy as np
from scipy.signal import savgol_filter

ap = argparse.ArgumentParser()
ap.add_argument("ows", nargs="?", default=os.path.join(here, "..", "ch01_pathlength.ows"))
ap.add_argument("--img", default=None)
args = ap.parse_args()
model = owsbuild.load_ows(args.ows)
W = lambda t: owsutil.by_title(model, t)[1]

wf = W("File")
print("File path loaded:", wf.recent_paths[0].abspath if wf.recent_paths else None,
      "| rows:", None if wf.data is None else len(wf.data))
sp = W("Spectra (all 21)")
X = sp.curveplot.data.X
print("Spectra: rows", len(X), "| colour:", getattr(sp.curveplot.feature_color, "name", None),
      "| max A per path:", {p: round(float(X[sp.curveplot.data.get_column("path_mm") == p].max()), 3)
                            for p in [0.25, 0.5, 1, 2, 5]})
sel = W("Select Rows (1 mm)")
print("Select Rows conditions:", sel.conditions)
cols = W("Select Columns (2258 nm)")
print("Select Columns features:", [v.name for v in cols.used_attrs], "target:", [v.name for v in cols.class_attrs])
ts = W("Test and Score")
print("Test and Score resampling:", ts.resampling, "(3 = Leave one out)")
print("  LOO row:", owsbuild.score_row(model))
sc = W("Scatter Plot")
print("Scatter:", sc.attr_x.name, "vs", sc.attr_y.name, "| regression line:", sc.graph.show_reg_line)
res = list(ts.learners.values())[0].results.value
print("  LOO per row (actual, predicted):",
      [(round(float(a), 2), round(float(p), 3)) for a, p in zip(res.actual, res.predicted[0])])
print("messages:", owsutil.messages(model))

if args.img:
    os.makedirs(args.img, exist_ok=True)
    owsbuild.shot(model, "Spectra (all 21)", os.path.join(args.img, "ch01_orange_spectra.png"), (1100, 700))
    owsbuild.shot(model, "Select Columns (2258 nm)", os.path.join(args.img, "ch01_orange_selectcolumns.png"), (900, 620))
    owsbuild.shot(model, "Test and Score", os.path.join(args.img, "ch01_orange_testscore.png"), (1000, 560))
    owsbuild.shot(model, "Scatter Plot", os.path.join(args.img, "ch01_orange_scatter.png"), (1000, 700), hide_legend=True)

# ---- numpy replica (same rows, same SG, mode="nearest" as Spectroscopy uses) ------------
import pandas as pd
csv = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(args.ows)), "../../data/spectra_wide_orange.csv"))
df = pd.read_csv(csv); wc = df.columns[4:]


def replica(path, col):
    s = df[(df["mC#path_mm"] == path) & (df["mC#replicate"] == 1)]
    D = savgol_filter(s[wc].values, 9, 2, deriv=2, axis=1, mode="nearest")
    d = D[:, list(wc).index(col)]; c = s["cC#acid_pct"].values
    e_dir, e_inv = [], []
    for k in range(len(c)):
        tr = [i for i in range(len(c)) if i != k]
        b, a = np.polyfit(d[tr], c[tr], 1); e_dir.append(a + b * d[k] - c[k])
        b, a = np.polyfit(c[tr], d[tr], 1); e_inv.append((d[k] - a) / b - c[k])
    rm = lambda e: float(np.sqrt(np.mean(np.square(e))))
    return float(np.corrcoef(d, c)[0, 1] ** 2), rm(e_dir), rm(e_inv)


print("\npath  band      | Orange LOO RMSE | Orange R2 (train) | numpy direct LOO | numpy inverse LOO (lesson) | R2")
for col in ["2258.067", "1669.901"]:
    owsbuild_cols = cols
    for path in [0.25, 0.5, 1, 2, 5]:
        sel.remove_all_rows()
        dom = sel.data.domain
        sel.add_row(dom["path_mm"], 0, [("%g" % path)]); sel.add_row(dom["replicate"], 0, ["1"])
        sel.conditions_changed(); sel.commit.now(); owsutil.settle(model, timeout=60)
        d = cols.data.domain
        if [v.name for v in cols.used_attrs] != [col]:
            allv = list(d.attributes) + list(d.class_vars)
            cols.used_attrs[:] = [d[col]]
            cols.available_attrs[:] = [v for v in d.attributes if v.name != col]
            cols.commit.now(); owsutil.settle(model, timeout=60)
        owsbuild.set_resampling(model, "Test and Score", 3)
        loo = owsbuild.score_row(model)
        owsbuild.set_resampling(model, "Test and Score", 4)   # 4 == Test on train data
        trn = owsbuild.score_row(model)
        r2, rd, ri = replica(path, col)
        print("%-5g %-9s | %15.3f | %17.4f | %16.3f | %26.3f | %.4f" % (
            path, col, float(loo["RMSE"]), float(trn["R2"]), rd, ri, r2))
