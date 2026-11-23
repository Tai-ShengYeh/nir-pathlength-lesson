"""Load ch05_prediction.ows headless, let the signals propagate, and print what the students
should see in the two Predictions widgets. Orange's DIRECT regression (c = a + b*signal) is
compared with a numpy replica of the same regression and with the lesson's INVERSE calibration
(signal = a + b*c, c0 = (y0 - a)/b; site/nir-pathlength/lib/ch05data.js). Also checks that
data/unknown_2258_orange.csv holds exactly Orange's own preprocessed 2258.067 values (x 10^3),
and shows the same-name-column trap (a plain column called 2258.067 -> prediction = intercept).

    python exercises/orange/tools/verify_ch05_ows.py [path/to/ch05_prediction.ows] [--img DIR]
"""
import os, sys, re, json, argparse
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import owsutil, owsbuild
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("ows", nargs="?", default=os.path.join(here, "..", "ch05_prediction.ows"))
ap.add_argument("--img", default=None)
ap.add_argument("--lesson", default=os.path.join(here, "..", "..", "..", "site", "nir-pathlength", "lib", "ch05data.js"))
args = ap.parse_args()
model = owsbuild.load_ows(args.ows)
W = lambda t: owsutil.by_title(model, t)[1]

for t in ["File", "File (unknown 2258)"]:
    wf = W(t)
    print(t, "loaded:", wf.recent_paths[0].abspath if wf.recent_paths else None,
          "| rows:", None if wf.data is None else len(wf.data))
for t in ["Select Rows (1 mm)", "Select Rows (hide 4.87 %)", "Select Rows (set = unknown)"]:
    print(t, "conditions:", W(t).conditions)
cols = W("Select Columns (2258 nm)")
print("Select Columns features:", [v.name for v in cols.used_attrs], "target:", [v.name for v in cols.class_attrs])


def pred_widget(title):
    w = W(title)
    slot = w.predictors[0]
    res = slot.results
    m = w.predictionsview.model()
    shown = [m.data(m.index(r, 0)) for r in range(m.rowCount())]
    data = w.data
    return dict(name=slot.name, predicted=[float(v) for v in res.predicted[0]],
                actual=[float(v) for v in data.Y], shown=shown,
                samples=list(data.get_column("sample")))


out = {}
for t in ["Predictions (4.87 % bottle)", "Predictions (y0 = 4.20)"]:
    out[t] = p = pred_widget(t)
    print(t, "->", p)
for t in ["Linear Regression (3 bottles)", "Linear Regression (4 bottles)"]:
    mdl = W(t).model
    print(t, "| trained on", len(W(t).data), "rows | direct c = a + b*x: a = %.4f, b = %.4f"
          % (mdl.intercept, mdl.coefficients[0]), "| x =", mdl.domain.attributes[0].name)
ct = W("Data Table (coefficients)")
print("Data Table (coefficients):", [(str(r["name"]), round(float(r[0]), 4)) for r in ct.input.table]
      if hasattr(ct, "input") else None)
print("messages:", owsutil.messages(model))

# ---- the small CSV holds Orange's own preprocessed values -----------------------------------
from Orange.data import Domain
dpre = cols.data.transform(Domain(list(cols.used_attrs), list(cols.class_attrs), list(cols.meta_attrs)))  # = its output
pre = {float(c): float(x) * 1e3 for c, x in zip(dpre.Y, dpre.get_column("2258.067"))}
du = W("File (unknown 2258)").data
cal = {float(c): float(x) for c, x, s in zip(du.Y, du.X[:, 0], du.get_column("set")) if s == "calibration"}
diff = max(abs(pre[c] - cal[c]) for c in pre)
print("unknown_2258_orange.csv calibration vs Orange Preprocess Spectra x1e3 (matched by acid_pct): max |diff| = %.2e" % diff)
assert sorted(pre) == sorted(cal) and diff < 1e-5

# ---- numpy replicas (direct) and the lesson (inverse) --------------------------------------
js = open(args.lesson, encoding="utf-8").read()
L = json.loads(js[js.index("{"): js.rindex("}") + 1])
c = np.array(sorted(pre)); y = np.array([pre[k] for k in c])          # x1e3, by concentration
tr = [i for i in range(4) if c[i] != 4.87]; k = list(c).index(4.87)
b, a = np.polyfit(y[tr], c[tr], 1); dirA = a + b * y[k]
b2, a2 = np.polyfit(c[tr], y[tr], 1); invA = (y[k] - a2) / b2
b, a = np.polyfit(y, c, 1); dirB = a + b * 4.20
b2, a2 = np.polyfit(c, y, 1); invB = (4.20 - a2) / b2
oA = out["Predictions (4.87 % bottle)"]["predicted"][0]; oB = out["Predictions (y0 = 4.20)"]["predicted"][0]
print("\ncase                         | Orange (direct) | numpy direct | numpy inverse | lesson inverse (ch05data.js)")
print("A: hide 4.87 %%, 3 bottles     | %15.3f | %12.3f | %13.3f | %.3f  (actual 4.87, signal %.4f)"
      % (oA, dirA, invA, L["hide3"]["pred"]["x0"], y[k]))
print("B: y0 = 4.20, 4 bottles       | %15.3f | %12.3f | %13.3f | %.3f" % (oB, dirB, invB, L["four"]["pred"]["x0"]))
assert abs(oA - dirA) < 1e-9 and abs(oB - dirB) < 1e-5   # B: CSV stores 6 decimals
assert abs(invA - L["hide3"]["pred"]["x0"]) < 5e-4 and abs(invB - L["four"]["pred"]["x0"]) < 5e-4
print("lesson Fieller 95 %% interval: A %s..%s (unbounded, g = %.2f); B %s..%s (g = %.3f); B with m = 3: %s..%s"
      % (L["hide3"]["pred"]["lo"], L["hide3"]["pred"]["hi"], L["hide3"]["calib"]["g"],
         L["four"]["pred"]["lo"], L["four"]["pred"]["hi"], L["four"]["calib"]["g"],
         L["four"]["pred_m3"]["lo"], L["four"]["pred_m3"]["hi"]))

# ---- the trap: a plain column named like the preprocessed one -------------------------------
import Orange
from Orange.data import Table, Domain, ContinuousVariable
m4 = Orange.regression.LinearRegressionLearner()(dpre)                 # 4 bottles, spectral domain
plain = Table.from_numpy(Domain([ContinuousVariable("2258.067")]), np.array([[4.20e-3]]))
print("\ntrap: 4-bottle model on Preprocess Spectra output, unknown = plain column '2258.067' = 0.0042 ->",
      "%.3f (intercept %.3f; correct direct prediction %.3f)" % (float(m4(plain)[0]), m4.intercept, dirB),
      "| value Orange actually used:", float(plain.transform(m4.domain).X[0, 0]))

if args.img:
    os.makedirs(args.img, exist_ok=True)
    owsbuild.shot(model, "Select Rows (hide 4.87 %)", os.path.join(args.img, "ch05_orange_selectrows.png"), (760, 420))
    owsbuild.shot(model, "Predictions (4.87 % bottle)", os.path.join(args.img, "ch05_orange_pred_hide.png"), (1000, 360))
    owsbuild.shot(model, "Predictions (y0 = 4.20)", os.path.join(args.img, "ch05_orange_pred_y0.png"), (1000, 360))
    owsbuild.shot(model, "Data Table (coefficients)", os.path.join(args.img, "ch05_orange_coefficients.png"), (760, 300))
