"""Load ch07_design.ows headless and compare what Orange shows with the lesson's numbers
(site/nir-pathlength/lib/ch07data.js, CH7.order) and with numpy.

    python exercises/orange/tools/verify_ch07_ows.py [path/to/ch07_design.ows] [--img DIR]
"""
import os, sys, json, argparse
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import owsutil, owsbuild
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("ows", nargs="?", default=os.path.join(here, "..", "ch07_design.ows"))
ap.add_argument("--img", default=None)
ap.add_argument("--lesson", default=os.path.join(here, "..", "..", "..", "site", "nir-pathlength", "lib", "ch07data.js"))
args = ap.parse_args()
model = owsbuild.load_ows(args.ows)
W = lambda t: owsutil.by_title(model, t)[1]
js = open(args.lesson, encoding="utf-8").read()
L = json.loads(js[js.index("{"): js.rindex("}") + 1])["order"]

wf = W("File")
print("File loaded:", wf.recent_paths[0].abspath if wf.recent_paths else None, "| rows:", len(wf.data))
assert len(wf.data) == 48
print("messages:", owsutil.messages(model))

exp = {"ascending": 1 + L["asc_bias"] / 100, "blocked": 1 + L["example_bias"] / 100}
print("lesson: asc_bias %.4f %% -> slope %.5f | example_bias %.5f %% -> slope %.5f"
      % (L["asc_bias"], exp["ascending"], L["example_bias"], exp["blocked"]))
for key in ("ascending", "blocked"):
    sel = W("Select Rows (%s)" % key)
    t = W("Data Table (slope, %s)" % key).input.table
    coef = {str(r["name"]): float(r[0]) for r in t}
    d = sel.data  # all 48; the matching 24 are what the regression sees
    lr = W("Linear Regression (%s)" % key).data
    x = lr.get_column("level_pct").astype(float); y = lr.Y.astype(float)
    order = {str(v) for v in lr.get_column("order_type")}
    b_np, a_np = np.polyfit(x, y, 1)
    print("%-9s order_type %s | rows %d | Orange slope %.5f intercept %.5f | numpy slope %.5f intercept %.5f | lesson slope %.5f"
          % (key, order, len(lr), coef["level_pct"], coef["intercept"], b_np, a_np, exp[key]))
    assert len(lr) == 24
    assert abs(coef["level_pct"] - b_np) < 1e-9 and abs(coef["level_pct"] - exp[key]) < 5e-5
print("run vs level correlation in the ascending design: %.4f ; blocked: %.4f" % tuple(
    np.corrcoef(np.arange(1, 25), W("Linear Regression (%s)" % k).data.get_column("level_pct").astype(float))[0, 1]
    for k in ("ascending", "blocked")))
for t in ("Scatter Plot (run vs level)", "Scatter Plot (run vs error)"):
    s = W(t); print(t, "| x", s.attr_x.name, "y", s.attr_y.name, "color", getattr(s.attr_color, "name", None))

if args.img:
    os.makedirs(args.img, exist_ok=True)
    P = lambda n: os.path.join(args.img, n)
    owsbuild.shot(model, "Data Table (design)", P("ch07_orange_table.png"), (1180, 640))
    owsbuild.shot(model, "Scatter Plot (run vs level)", P("ch07_orange_run_level.png"), (900, 640), hide_legend=False)
    owsbuild.shot(model, "Scatter Plot (run vs error)", P("ch07_orange_run_error.png"), (900, 640))
    owsbuild.shot(model, "Data Table (slope, ascending)", P("ch07_orange_slope_ascending.png"), (760, 200))
    owsbuild.shot(model, "Data Table (slope, blocked)", P("ch07_orange_slope_blocked.png"), (760, 200))
    owsbuild.shot(model, "Linear Regression (ascending)", P("ch07_orange_lr.png"), None)
