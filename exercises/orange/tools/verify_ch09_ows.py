"""Load ch09_validation.ows headless and compare what Orange shows with the lesson's numbers
(site/nir-pathlength/lib/ch09data.js, CH9.spec.r2 and CH9.acc) and with numpy.

    python exercises/orange/tools/verify_ch09_ows.py [path/to/ch09_validation.ows] [--img DIR]

Orange regresses reference = a + b * d2_2258 (direct): R2 (Test on train) equals the lesson's R2 exactly; the leave-one-out
RMSE (in % acid) is a DIRECT-regression number and differs from the lesson's inverse-calibration RMSECV (checked below)."""
import os, sys, json, argparse
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import owsutil, owsbuild
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("ows", nargs="?", default=os.path.join(here, "..", "ch09_validation.ows"))
ap.add_argument("--img", default=None)
ap.add_argument("--lesson", default=os.path.join(here, "..", "..", "..", "site", "nir-pathlength", "lib", "ch09data.js"))
args = ap.parse_args()
model = owsbuild.load_ows(args.ows)
W = lambda t: owsutil.by_title(model, t)[1]
js = open(args.lesson, encoding="utf-8").read()
L = json.loads(js[js.index("{"): js.rindex("}") + 1])

wf = W("File")
print("File loaded:", wf.recent_paths[0].relpath if wf.recent_paths else None, "| rows:", len(wf.data))
assert len(wf.data) == 4
print("messages:", owsutil.messages(model))
assert owsutil.messages(model) == []

d = wf.data
sig = np.array([float(v) for v in d.Y])
col = lambda n: np.array([float(v) for v in d.get_column(n)])
refs = {"label": col("label_pct"), "Atago": col("atago_1102"), "titration": col("titr_pct"), "PAL-30S": col("pal30s_pct")}
key2258 = {"label": L["spec"]["r2"]["2258"]["label"], "Atago": L["spec"]["r2"]["2258"]["atago"],
           "titration": L["spec"]["r2"]["2258"]["titr"], "PAL-30S": L["spec"]["r2"]["2258"]["pal30s"]}
key_rmsecv = {"Atago": L["acc"]["rmsecv_atago"], "titration": L["acc"]["rmsecv_titr"], "PAL-30S": L["acc"]["rmsecv_pal30s"], "label": 0.6769}  # label: ch08 nirfit


def loo(x, y, direct):
    e = []
    for i in range(len(x)):
        m = np.arange(len(x)) != i
        if direct:
            b, a = np.polyfit(y[m], x[m], 1); e.append(a + b * y[i] - x[i])      # reference = a + b * signal  (Orange)
        else:
            b, a = np.polyfit(x[m], y[m], 1); e.append((y[i] - a) / b - x[i])    # signal = a + b * reference, back-calculated (lesson)
    return float(np.sqrt(np.mean(np.square(e))))


print("\n%-10s | Orange R2 train | lesson R2 (CH9) || Orange LOO RMSE | numpy direct LOO | numpy inverse LOO | lesson RMSECV" % "reference")
for nm, x in refs.items():
    hl, rl = owsutil.test_score_table(W("Test and Score LOO (%s)" % nm)); ht, rt = owsutil.test_score_table(W("Test and Score train (%s)" % nm))
    rmse = dict(zip(hl, rl[0]))["RMSE"]; r2 = dict(zip(ht, rt[0]))["R2"]
    r2n = float(np.corrcoef(x, sig)[0, 1] ** 2); dn, inv = loo(x, sig, True), loo(x, sig, False)
    print("%-10s | %15.5f | %16.4f || %15.4f | %16.4f | %17.4f | %13.4f" % (nm, r2, key2258[nm], rmse, dn, inv, key_rmsecv[nm]))
    assert abs(r2 - r2n) < 1e-6 and abs(r2 - key2258[nm]) <= 1e-3 * key2258[nm]      # R2 identical to the lesson
    assert abs(rmse - dn) < 1e-5                                                      # Orange LOO RMSE = numpy direct regression
    assert abs(inv - key_rmsecv[nm]) <= 1e-3 * key_rmsecv[nm]                         # lesson (inverse) RMSECV reproduced by numpy, not by Orange
# 1670 nm: R2 is also checkable from the same CSV (swap Features to d2_1670 in Select Columns)
s70 = col("d2_1670")
for nm, k in (("Atago", "atago"), ("titration", "titr"), ("label", "label"), ("PAL-30S", "pal30s")):
    r = float(np.corrcoef(refs[nm], s70)[0, 1] ** 2); assert abs(r - L["spec"]["r2"]["1670"][k]) <= 1e-3 * L["spec"]["r2"]["1670"][k]
print("1670 nm R2 (numpy from the CSV) matches CH9.spec.r2['1670'] for all four references")
print("OK: R2 identical to the lesson; Orange LOO RMSE (direct) equals numpy direct, differs from the lesson's inverse RMSECV")

if args.img:
    os.makedirs(args.img, exist_ok=True)
    P = lambda n: os.path.join(args.img, n)
    owsbuild.shot(model, "Data Table (validation)", P("ch09_orange_table.png"), (1250, 260))
    owsbuild.shot(model, "Select Columns (titration)", P("ch09_orange_selectcolumns.png"), (900, 620))
    owsbuild.shot(model, "Test and Score LOO (titration)", P("ch09_orange_loo_titr.png"), (1300, 640))
    owsbuild.shot(model, "Test and Score LOO (PAL-30S)", P("ch09_orange_loo_pal30s.png"), (1300, 640))
    owsbuild.shot(model, "Test and Score train (titration)", P("ch09_orange_train_titr.png"), (1300, 640))
    owsbuild.shot(model, "Test and Score train (PAL-30S)", P("ch09_orange_train_pal30s.png"), (1300, 640))
