"""Load ch02_references.ows headless, let the signals propagate, and print what the students
should see; then switch the target (4 references) and the wavelength (2258 / 1670 nm) on the live
Select Columns widget and read Test and Score ("Test on train data" R2, and Leave one out RMSE).

    python exercises/orange/tools/verify_ch02_ows.py [path/to/ch02_references.ows] [--img DIR]
"""
import os, sys, argparse
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import owsutil, owsbuild
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("ows", nargs="?", default=os.path.join(here, "..", "ch02_references.ows"))
ap.add_argument("--img", default=None)
args = ap.parse_args()
model = owsbuild.load_ows(args.ows)
W = lambda t: owsutil.by_title(model, t)[1]

wf = W("File")
print("File path loaded:", wf.recent_paths[0].abspath if wf.recent_paths else None,
      "| rows:", None if wf.data is None else len(wf.data))
cols = W("Select Columns (2258 nm)")
d = cols.data
print("Formula output metas:", [m.name for m in d.domain.metas])
for nm in ["sample", "label_pct", "atago_1102", "atago_1103", "naoh_ml", "titration", "ratio"]:
    v = d.get_column(nm)
    print("  %-10s" % nm, [x if isinstance(x, str) else round(float(x), 4) for x in v])
print("Select Columns features:", [v.name for v in cols.used_attrs], "target:", [v.name for v in cols.class_attrs])
ts = W("Test and Score")
print("Test and Score resampling:", ts.resampling, "(4 = Test on train data)")
print("  row:", owsbuild.score_row(model))
sc = W("Scatter Plot")
print("Scatter:", sc.attr_x.name, "vs", sc.attr_y.name, "| label:", getattr(sc.attr_label, "name", None),
      "| regression line:", sc.graph.show_reg_line)
print("messages:", owsutil.messages(model))

if args.img:
    os.makedirs(args.img, exist_ok=True)
    owsbuild.shot(model, "Formula (titration, ratio)", os.path.join(args.img, "ch02_orange_formula.png"), (900, 520))
    owsbuild.shot(model, "Data Table", os.path.join(args.img, "ch02_orange_table.png"), (1400, 300))
    owsbuild.shot(model, "Select Columns (2258 nm)", os.path.join(args.img, "ch02_orange_selectcolumns.png"), (900, 620))
    owsbuild.shot(model, "Test and Score", os.path.join(args.img, "ch02_orange_testscore.png"), (1000, 560))
    owsbuild.shot(model, "Scatter Plot", os.path.join(args.img, "ch02_orange_scatter.png"), (1000, 700))


def inverse_direct(dv, c):
    e_dir, e_inv = [], []
    for k in range(len(c)):
        tr = [i for i in range(len(c)) if i != k]
        b, a = np.polyfit(dv[tr], c[tr], 1); e_dir.append(a + b * dv[k] - c[k])
        b, a = np.polyfit(c[tr], dv[tr], 1); e_inv.append((dv[k] - a) / b - c[k])
    rm = lambda e: float(np.sqrt(np.mean(np.square(e))))
    return float(np.corrcoef(dv, c)[0, 1] ** 2), rm(e_dir), rm(e_inv)


KEY = {"label_pct": (0.9955, 0.9983), "titration": (0.9898, 0.9866), "atago_1103": (0.9555, 0.9633),
       "atago_1102": (0.9342, 0.9448)}
print("\ntarget      band      | Orange R2 (train) | lesson R2 | Orange LOO RMSE | numpy direct | numpy inverse")
for band, kb in [("2258.067", 0), ("1669.901", 1)]:
    for tgt in ["label_pct", "titration", "atago_1103", "atago_1102"]:
        dom = cols.data.domain
        feats = [dom[band]]; t = [dom[tgt]]
        mets = [m for m in dom.metas if m.name != tgt]
        cols.used_attrs[:] = feats; cols.class_attrs[:] = t; cols.meta_attrs[:] = mets
        cols.available_attrs[:] = [a for a in dom.attributes if a.name != band]
        cols.commit.now(); owsutil.settle(model, timeout=60)
        owsbuild.set_resampling(model, "Test and Score", 4)
        trn = owsbuild.score_row(model)
        owsbuild.set_resampling(model, "Test and Score", 3)
        loo = owsbuild.score_row(model)
        out = W("Linear Regression").data
        r2, rd, ri = inverse_direct(out.X[:, 0], out.Y)
        ok = abs(float(trn["R2"]) - KEY[tgt][kb]) < 5e-5
        print("%-11s %-9s | %17.4f | %9.4f %s | %15.3f | %12.3f | %13.3f" % (
            tgt, band, float(trn["R2"]), KEY[tgt][kb], "OK " if ok else "BAD", float(loo["RMSE"]), rd, ri))

# What-if on the live widgets (not saved): NaOH factor f off by 5 % -> titration x 1.05, R2 unchanged.
from Orange.widgets.data.owfeatureconstructor import ContinuousDescriptor
dom = cols.data.domain
cols.used_attrs[:] = [dom["2258.067"]]; cols.class_attrs[:] = [dom["titration"]]
cols.meta_attrs[:] = [m for m in dom.metas if m.name != "titration"]
cols.available_attrs[:] = [a for a in dom.attributes if a.name != "2258.067"]
cols.commit.now(); owsutil.settle(model, timeout=60)
fw = W("Formula (titration, ratio)")
fw.featuremodel[0] = ContinuousDescriptor("titration", "naoh_ml * 1.07875 * 1.05 * 0.24", 3, True)
fw.apply(); owsutil.settle(model, timeout=60, quiet=2)
owsbuild.set_resampling(model, "Test and Score", 4)
print("what-if f x 1.05: titration =", [round(float(x), 4) for x in W("Linear Regression").data.Y],
      "| R2 (train) =", round(float(owsbuild.score_row(model)["R2"]), 4))
