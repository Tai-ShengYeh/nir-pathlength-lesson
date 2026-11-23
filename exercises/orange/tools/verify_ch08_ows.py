"""Load ch08_references.ows headless and compare what Orange shows with the lesson's numbers
(site/nir-pathlength/lib/ch08data.js, CH8.factor / CH8.agree) and with numpy.

    python exercises/orange/tools/verify_ch08_ows.py [path/to/ch08_references.ows] [--img DIR]
"""
import os, sys, json, argparse
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import owsutil, owsbuild
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("ows", nargs="?", default=os.path.join(here, "..", "ch08_references.ows"))
ap.add_argument("--img", default=None)
ap.add_argument("--lesson", default=os.path.join(here, "..", "..", "..", "site", "nir-pathlength", "lib", "ch08data.js"))
args = ap.parse_args()
model = owsbuild.load_ows(args.ows)
W = lambda t: owsutil.by_title(model, t)[1]
js = open(args.lesson, encoding="utf-8").read()
L = json.loads(js[js.index("{"): js.rindex("}") + 1])

wf = W("File")
print("File loaded:", wf.recent_paths[0].relpath if wf.recent_paths else None, "| rows:", len(wf.data))
assert len(wf.data) == 4
print("messages:", owsutil.messages(model))

# --- Formula: titr_pct_calc = 0.2589 * V ; gap in mL vs the lesson
fd = W("Data Table (gap in mL)").input.table
col = lambda n: np.array([float(v) for v in fd.get_column(n)])
calc, need, gap = col("titr_pct_calc"), col("need_mL_label"), col("gap_mL")
print("titr_pct_calc", np.round(calc, 4), "| CSV titr_pct", np.round(col("titr_pct"), 4))
assert np.allclose(calc, col("titr_pct"), atol=1e-4)           # 0.2589 vs unrounded K
assert np.allclose(calc[1:], L["refs"]["titr"], rtol=2e-3)
print("need_mL_label", np.round(need, 3), "lesson", L["agree"]["label"]["need_mL"])
print("gap_mL        ", np.round(gap, 3), "lesson", L["agree"]["label"]["gap_mL"])
assert np.allclose(need[1:], L["agree"]["label"]["need_mL"], atol=0.01)
assert np.allclose(gap[1:], L["agree"]["label"]["gap_mL"], atol=0.01)

# --- regressions: slope (Orange) vs numpy vs lesson
fa = L["factor"]
x = col("titr_pct"); y = np.array([float(v) for v in fd.Y]); xs = x * 0.97
for key, xv, exp in (("titration", x, fa["slope"]), ("titration x 0.97", xs, fa["slope_scaled"])):
    t = W("Data Table (slope, %s)" % key).input.table
    coef = {str(r["name"]): float(r[0]) for r in t}
    b = [v for k, v in coef.items() if k != "intercept"][0]
    bn, an = np.polyfit(xv, y, 1)
    hdr, rows = owsutil.test_score_table(W("Test and Score (%s)" % key))
    r2 = dict(zip(hdr, rows[0]))["R2"]
    r2n = float(np.corrcoef(xv, y)[0, 1] ** 2)
    print("%-17s Orange slope %.5f intercept %.5f R2 %.5f | numpy slope %.5f R2 %.5f | lesson slope %.4f R2 %.4f"
          % (key, b, coef["intercept"], r2, bn, r2n, exp, fa["r2"]))
    assert abs(b - bn) < 1e-5 and abs(r2 - r2n) < 1e-6 and abs(b - exp) < 5e-4 and abs(r2 - fa["r2"]) < 5e-4
print("slope ratio scaled/unscaled: %.5f (1/0.97 = %.5f)" % (fa["slope_scaled"] / fa["slope"], 1 / 0.97))
for t in ("Scatter Plot (d2 vs label)", "Scatter Plot (d2 vs Atago 11-02)", "Scatter Plot (d2 vs Atago 11-03)",
          "Scatter Plot (d2 vs titration)"):
    s = W(t); print(t, "| x", s.attr_x.name, "y", s.attr_y.name)

if args.img:
    os.makedirs(args.img, exist_ok=True)
    P = lambda n: os.path.join(args.img, n)
    owsbuild.shot(model, "Data Table (references)", P("ch08_orange_table.png"), (1000, 260))
    owsbuild.shot(model, "Data Table (gap in mL)", P("ch08_orange_gap.png"), (1750, 260))
    owsbuild.shot(model, "Scatter Plot (d2 vs titration)", P("ch08_orange_scatter_titr.png"), (800, 560), hide_legend=True)
    owsbuild.shot(model, "Scatter Plot (d2 vs Atago 11-03)", P("ch08_orange_scatter_atago.png"), (800, 560), hide_legend=True)
    owsbuild.shot(model, "Data Table (slope, titration)", P("ch08_orange_slope_titr.png"), (760, 200))
    owsbuild.shot(model, "Data Table (slope, titration x 0.97)", P("ch08_orange_slope_f097.png"), (760, 200))
    owsbuild.shot(model, "Test and Score (titration x 0.97)", P("ch08_orange_testscore.png"), (1400, 380))
