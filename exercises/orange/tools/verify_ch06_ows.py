"""Load ch06_data_quality.ows headless, let the signals propagate, and print what the students
should see, next to the lesson's numbers (site/nir-pathlength/lib/ch06data.js):

  1. water-band height ratios (Integrate Spectra -> Formula -> Group by median -> Merge Data -> Formula)
     vs the lesson's exact definition (window mean 1440-1460 minus window mean 1580-1650, numpy);
  2. the 0.25 mm spectra in the Spectra widget;
  3. what the Outliers widget flags with its defaults and with the other three methods, Isolation
     Forest repeated, the Mahalanobis column of the covariance method, and the PCA scores / a numpy
     PCA-Q ranking (Orange has no Q-residual widget).

    python exercises/orange/tools/verify_ch06_ows.py [path/to/ch06_data_quality.ows] [--img DIR]
"""
import os, sys, json, argparse
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import owsutil, owsbuild
import numpy as np
from AnyQt.QtTest import QTest

ap = argparse.ArgumentParser()
ap.add_argument("ows", nargs="?", default=os.path.join(here, "..", "ch06_data_quality.ows"))
ap.add_argument("--img", default=None)
ap.add_argument("--lesson", default=os.path.join(here, "..", "..", "..", "site", "nir-pathlength", "lib", "ch06data.js"))
args = ap.parse_args()
model = owsbuild.load_ows(args.ows)
W = lambda t: owsutil.by_title(model, t)[1]
js = open(args.lesson, encoding="utf-8").read()
L = json.loads(js[js.index("{"): js.rindex("}") + 1])

wf = W("File")
print("File loaded:", wf.recent_paths[0].abspath if wf.recent_paths else None, "| rows:", len(wf.data))
print("Select Rows (replicate = 1) conditions:", W("Select Rows (replicate = 1)").conditions)
integ = W("Integrate Spectra")
print("  integrals:", integ.preprocessormodel.rowCount(), "(expected 3: 1440-1460 area, 1580-1650 area, 1440-1460 peak)")
g = W("Group by path")
print("Group by:", [v.name for v in g.gb_attrs], {k.name: sorted(v) for k, v in g.aggregations.items() if v})
mgw = W("Merge Data")
print("Merge Data row matching:", [(a.name, b.name) for a, b in mgw.attr_boxes.current_state()])
print("messages:", owsutil.messages(model))

# ---- 1. water-band ratios --------------------------------------------------------------------
t = W("Data Table (ratios)").input.table
path = t.get_column("path_mm").astype(float); acid = t.Y.astype(float)
smp = [str(s) for s in t.get_column("sample")]
ratio = t.get_column("ratio").astype(float); flag = t.get_column("flag").astype(float)
height = t.get_column("height").astype(float); peak = t.get_column("peak").astype(float)
med = t.get_column("height - Median").astype(float)

# exact definition, straight from the CSV (window means, numpy)
import pandas as pd
d = pd.read_csv(os.path.join(os.path.dirname(args.ows), "..", "..", "data", "spectra_wide_orange.csv"))
wl = np.array([float(c) for c in d.columns[4:]]); X = d.iloc[:, 4:].values
r1 = (d["mC#replicate"] == 1).values
WB = (wl >= 1440) & (wl <= 1460); FL = (wl >= 1580) & (wl <= 1650)
exact_h = X[:, WB].mean(1) - X[:, FL].mean(1)
Pp = d["mC#path_mm"].values; Cc = d["cC#acid_pct"].values
print("\nwindow points: 1440-1460 ->", wl[WB], "| 1580-1650 -> %d points %.3f..%.3f" % (FL.sum(), wl[FL][0], wl[FL][-1]))
print("\npath | acid  | Orange height | exact height | Orange ratio | exact ratio | lesson ratio | peak (1440-1460) | flag")
worst = 0.0; worst_usable = 0.0; flagged = []
for p in [0.25, 0.5, 1.0, 2.0, 5.0]:
    sel = (path == p)
    ex = {c: h for c, h in zip(Cc[(Pp == p) & r1], exact_h[(Pp == p) & r1])}
    exmed = np.median(list(ex.values()))
    for i in np.where(sel)[0]:
        c = round(acid[i], 2)
        les = L["wb"][str(p)]["ratio"][[0.0, 3.32, 4.87, 7.45].index(c)]
        er = ex[c] / exmed
        dd = max(abs(ratio[i] - er), abs(ratio[i] - les))
        worst = max(worst, dd)
        if p <= 0.5:
            worst_usable = max(worst_usable, dd)
        print("%4g | %5.2f | %13.4f | %12.4f | %12.4f | %11.4f | %12.4g | %16.3f | %d"
              % (p, c, height[i], ex[c], ratio[i], er, les, peak[i], flag[i]))
        if flag[i] == 1:
            flagged.append((p, c, ratio[i]))
print("max |Orange ratio - exact or lesson ratio|: 0.25 and 0.5 mm = %.5f, all 20 (incl. saturated 2 and 5 mm) = %.5f" % (worst_usable, worst))
assert worst_usable < 6e-4
print("flagged by Orange (usable and |ratio-1| > 0.10):", flagged)
assert [(p, c) for p, c, _ in flagged] == [(0.25, 4.87), (0.5, 3.32)]
print("lesson flags: 0.25 mm index 2 (4.87 %) ->", L["wb"]["0.25"]["flag"], "| 0.5 mm ->", L["wb"]["0.5"]["flag"],
      "| usable:", {k: v["usable"] for k, v in L["wb"].items()}, "| peaks:", {k: v["peak"] for k, v in L["wb"].items()})
print("without the saturation rule (|ratio-1| > 0.10 only) Orange would also flag:",
      [(float(path[i]), round(float(acid[i]), 2), round(float(ratio[i]), 4)) for i in range(20)
       if abs(ratio[i] - 1) > 0.10 and flag[i] == 0])
tf = W("Data Table (flagged)").input.table
print("Data Table (flagged) rows:", [(float(r["path_mm"]), str(r["sample"]), round(float(r["ratio"]), 4)) for r in tf])

# ---- 2. the 0.25 mm spectra -------------------------------------------------------------------
sp = W("Spectra (0.25 mm)")
sd = sp.curveplot.data
print("\nSpectra (0.25 mm): rows", len(sd), "| colour:", getattr(sp.curveplot.feature_color, "name", None),
      "| samples", [str(s) for s in sd.get_column("sample")])
xs = np.array([float(a.name) for a in sd.domain.attributes]); Xs = sd.X
wbm = (xs >= 1440) & (xs <= 1460)
print("  mean A in 1440-1460 nm:", {str(s): round(float(Xs[i, wbm].mean()), 4) for i, s in enumerate(sd.get_column("sample"))})

# ---- 3. Outliers widget -----------------------------------------------------------------------
o = W("Outliers")
names = ["One Class SVM", "Covariance Estimator", "Local Outlier Factor", "Isolation Forest"]
dt = W("Data Table (Outliers)")


def outlier_state(label):
    owsutil.settle(model, timeout=60)
    tb = dt.input.table
    col = tb.get_column("Outlier")
    dom = [v for v in tb.domain.metas + tb.domain.attributes if v.name == "Outlier"][0]
    vals = [dom.str_val(v) if hasattr(dom, "str_val") else str(v) for v in col]
    labs = {v for v in vals}
    yes = [i for i, v in enumerate(vals) if v in ("0", "0.0", "No", "no")]   # filled below once the coding is known
    return tb, vals


tb, vals = outlier_state("default")
print("\nOutliers default method:", names[o.outlier_method], "| editor", {k: getattr(o.current_editor, k, None)
      for k in ("cont", "n_neighbors", "metric_index")}, "| Outlier column values:", sorted(set(vals)))


def flagged_of(tb, vals):
    return [(float(tb.get_column("path_mm")[i]), str(tb.get_column("sample")[i])) for i, v in enumerate(vals) if v == "Yes"]


print("  flagged with defaults (LOF, 20 neighbours, 10 %%): %d -> %s" % (len(flagged_of(tb, vals)), flagged_of(tb, vals)))
rows = {}
for k in (0, 1, 2, 3):
    o.outlier_method = k
    o._OWOutliers__method_changed()
    tb, vals = outlier_state(names[k])
    fl_ = flagged_of(tb, vals)
    ed = o.current_editor
    prm = {a: getattr(ed, a, None) for a in ("cont", "nu", "gamma", "n_neighbors", "support_fraction", "empirical_covariance", "replicable")
           if getattr(ed, a, None) is not None}
    print("  %-22s params %-60s flagged %2d: %s" % (names[k], prm, len(fl_), fl_))
    rows[k] = (tb, fl_)
tb = rows[1][0]
mh = tb.get_column("Mahalanobis").astype(float)
print("  covariance method, Mahalanobis column: min %.4f max %.4f (n-1 = 19, all distances equal when features >> samples)"
      % (mh.min(), mh.max()))
o.outlier_method = 3; o._OWOutliers__method_changed()
runs = []
for i in range(8):
    o.outlier_method = 2; o._OWOutliers__method_changed(); owsutil.settle(model, timeout=60)   # switch away and back = refit
    o.outlier_method = 3; o._OWOutliers__method_changed(); owsutil.settle(model, timeout=60)
    tb, vals = outlier_state("isf"); runs.append(tuple(flagged_of(tb, vals)))
print("  Isolation Forest, 8 refits, 'Replicable training' off:", sorted(set(runs)))
o.isf_editor.replicable = True; o.isf_editor.param_changed.emit(); owsutil.settle(model, timeout=60)
print("  with Replicable training on:", flagged_of(*outlier_state("isf")))
o.isf_editor.replicable = False
o.outlier_method = 2; o._OWOutliers__method_changed(); outlier_state("restore")

# ---- PCA ----------------------------------------------------------------------------------------
pw = W("PCA (20 raw spectra)")
print("\nPCA widget: explained variance ratio PC1 %.4f PC2 %.4f" % tuple(pw._variance_ratio[:2]))
ps = W("Scatter Plot (PCA)")
pt = ps.data
print("  PC1 / PC2 per spectrum (path, acid):")
for i in np.argsort(pt.get_column("PC1").astype(float)):
    print("   %4g mm %5.2f %%  PC1 %8.3f  PC2 %8.3f" % (float(pt.get_column("path_mm")[i]), float(pt.Y[i]),
                                                        float(pt.get_column("PC1")[i]), float(pt.get_column("PC2")[i])))
Xr = X[r1]; Pr = Pp[r1]; Cr = Cc[r1]
mu = Xr.mean(0); U, S, Vt = np.linalg.svd(Xr - mu, full_matrices=False)
Sc = (Xr - mu) @ Vt[:2].T; Q = ((Xr - mu - Sc @ Vt[:2]) ** 2).sum(1)
order = np.argsort(-Q)
print("numpy PCA (2 comps) Q-residual ranking (top 8):", [(float(Pr[i]), float(Cr[i]), round(float(Q[i]), 3)) for i in order[:8]])
rk = lambda p, c: int(np.where((Pr[order] == p) & (Cr[order] == c))[0][0]) + 1
print("  rank of 0.25 mm / 4.87 %%: %d | of 0.5 mm / 3.32 %%: %d | lesson: %s"
      % (rk(0.25, 4.87), rk(0.5, 3.32), L["pcaQ"]["rank_of_flagged"]))
print("  lesson top 4 by Q:", [(L["pcaQ"]["path"][i], L["pcaQ"]["conc"][i]) for i in L["pcaQ"]["rank"][:4]])
assert rk(0.25, 4.87) == L["pcaQ"]["rank_of_flagged"]["0.25-4.87"] and rk(0.5, 3.32) == L["pcaQ"]["rank_of_flagged"]["0.5-3.32"]
print("messages:", owsutil.messages(model))

if args.img:
    os.makedirs(args.img, exist_ok=True)
    P = lambda n: os.path.join(args.img, n)
    owsbuild.shot(model, "Integrate Spectra", P("ch06_orange_integrate.png"), (1100, 640))
    owsbuild.shot(model, "Data Table (ratios)", P("ch06_orange_ratios.png"), (1280, 640))
    owsbuild.shot(model, "Spectra (0.25 mm)", P("ch06_orange_spectra.png"), (1100, 640))
    owsbuild.shot(model, "Outliers", P("ch06_orange_outliers.png"), (420, 420))
    owsbuild.shot(model, "Scatter Plot (PCA)", P("ch06_orange_pca.png"), (1000, 700), hide_legend=True)
