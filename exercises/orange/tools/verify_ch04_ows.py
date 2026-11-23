"""Load ch04_pca_pls.ows headless the way Orange.canvas.run does, let the signals propagate,
and print the numbers the students should see.

    python exercises/orange/tools/verify_ch04_ows.py [path/to/ch04_pca_pls.ows] [--img DIR]
"""
import os, sys, argparse
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import owsutil
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("ows", nargs="?", default=os.path.join(here, "..", "ch04_pca_pls.ows"))
ap.add_argument("--img", default=None)
args = ap.parse_args()
ows = os.path.abspath(args.ows)

reg, model = owsutil.registry_and_model(os.path.dirname(ows))
from orangecanvas.scheme import signalmanager
sm = model.findChild(signalmanager.SignalManager)
sm.pause()
with open(ows, "rb") as fh:
    model.load_from(fh, registry=reg)
for n in model.nodes:
    model.widget_for_node(n)
sm.resume()
ok = owsutil.settle(model, timeout=120)
print("settled:", ok)

W = lambda t: owsutil.by_title(model, t)[1]
wf = W("File")
print("File path loaded:", wf.recent_paths[0].abspath if wf.recent_paths else None,
      "| rows:", None if wf.data is None else len(wf.data))

pA = W("PCA (20 raw spectra)")
print("Branch A  PCA input rows:", len(pA.data), "| normalize:", pA.normalize,
      "| explained variance:", np.round(pA._variance_ratio[:2], 3))
sc = W("Scatter Plot")
d = sc.data
pc1 = d.get_column("PC1"); path = d.get_column("path_mm"); conc = d.get_column("acid_pct")
print("Scatter:", sc.attr_x.name, "vs", sc.attr_y.name, "| colour:", getattr(sc.attr_color, "name", None),
      "| |r(PC1,path)| =", round(abs(np.corrcoef(pc1, path)[0, 1]), 3),
      "| |r(PC1,conc)| =", round(abs(np.corrcoef(pc1, conc)[0, 1]), 3))

pp = W("Preprocess Spectra")
out = pp.data  # input
pB = W("PCA (1 mm, 2nd derivative)")
dB = pB.data
wl = np.array([float(a.name) for a in dB.domain.attributes])
print("Branch B  rows:", len(dB), "| conc:", list(dB.get_column("acid_pct")),
      "| n wavelengths:", len(wl), "| ranges:", round(wl.min(), 1), "-", round(wl[wl < 2000].max(), 1),
      "and", round(wl[wl > 2000].min(), 1), "-", round(wl.max(), 1))
print("Preprocess list:", pp.storedsettings)
print("Branch B  PCA normalize:", pB.normalize, "| explained variance:", np.round(pB._variance_ratio[:2], 3))
pl = W("PLS")
print("PLS n_components:", pl.n_components, "| scale:", pl.scale)
ts = W("Test and Score")
print("Test and Score resampling:", ts.resampling, "(3 = Leave one out)")
hdr, rows = owsutil.test_score_table(ts)
print("  header:", hdr)
for r in rows:
    print("  row:", r)
res = list(ts.learners.values())[0].results.value
print("  LOO per row (actual, predicted, error):",
      [(round(float(a), 2), round(float(p), 3), round(float(p - a), 3)) for a, p in zip(res.actual, res.predicted[0])],
      "| row_indices:", list(res.row_indices))
print("messages:", owsutil.messages(model))

if args.img:
    os.makedirs(args.img, exist_ok=True)
    for title, fn in [("Preprocess Spectra", "ch04_orange_preprocess.png"), ("PLS", "ch04_orange_pls.png"),
                      ("Test and Score", "ch04_orange_testscore.png"), ("Scatter Plot", "ch04_orange_scatter.png"),
                      ("Select Rows (1 mm)", "ch04_orange_selectrows.png"), ("PCA (20 raw spectra)", "ch04_orange_pca.png")]:
        w = W(title)
        w.resize(w.sizeHint().expandedTo(w.minimumSizeHint()))
        if title == "Scatter Plot":  # legend box sits on the 5 mm cluster; hide it for the picture only
            from AnyQt.QtWidgets import QCheckBox
            for cb in w.findChildren(QCheckBox):
                if cb.text() == "Show legend" and cb.isChecked():
                    cb.click()
        w.show()
        big = {"Scatter Plot": (1000, 820), "Preprocess Spectra": (1100, 900), "PCA (20 raw spectra)": (1000, 560)}
        if title in big:
            w.resize(*big[title])
        owsutil.settle(model, timeout=10, quiet=0.5)
        w.grab().save(os.path.join(args.img, fn)); w.hide()
        print("saved", fn)

# What-if checks on the live widgets (not saved): PLS widget default "Scale" on, and 2 components.
for label, nc, scale in [("Scale features and target CHECKED, 1 comp", 1, True), ("unscaled, 2 comps", 2, False),
                         ("back to workflow setting", 1, False)]:
    pl.controls.n_components.setValue(nc)
    if pl.controls.scale.isChecked() != scale:
        pl.controls.scale.click()
    pl.apply()
    owsutil.settle(model, timeout=60)
    hdr, rows = owsutil.test_score_table(ts)
    print("variant:", label, "-> RMSE", round(float(rows[0][hdr.index("RMSE")]), 3))
