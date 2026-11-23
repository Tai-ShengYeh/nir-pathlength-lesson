"""Build exercises/orange/ch07_design.ows with Orange's own scheme API (Orange 3.40.0), same approach as
build_ch06_ows.py.   python exercises/orange/tools/make_ch07_design.py   (once, writes ch07_design.csv)
                     python exercises/orange/tools/build_ch07_ows.py

File (ch07_design.csv, 48 rows) -> Data Table (design), Scatter Plot (run vs level_pct, colour order_type),
  Formula (error = signal - level_pct) -> Scatter Plot (run vs error)
  Select Rows (order_type = 由低到高)  -> Linear Regression -> Coefficients -> Data Table (slope, ascending)
  Select Rows (order_type = 分3塊打亂) -> Linear Regression -> Coefficients -> Data Table (slope, blocked)
Linear Regression has no regularisation (Ordinary), signal = a + b * level_pct, so b is the lesson's slope."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import owsutil
from owsbuild import Builder
from Orange.widgets.data.owfeatureconstructor import ContinuousDescriptor

ORANGE_DIR = os.path.normpath(os.path.join(here, ".."))
OWS = os.path.join(ORANGE_DIR, "ch07_design.ows")
B = Builder(ORANGE_DIR, "ch07_design.csv")

SR = "Orange.widgets.data.owselectrows.OWSelectRows"
FC = "Orange.widgets.data.owfeatureconstructor.OWFeatureConstructor"
TAB = "Orange.widgets.data.owtable.OWTable"
SCAT = "Orange.widgets.visualize.owscatterplot.OWScatterPlot"
LR = "Orange.widgets.model.owlinearregression.OWLinearRegression"

f = B.file_node((60, 220))
tab = B.node(TAB, "Data Table (design)", (260, 40))
sc1 = B.node(SCAT, "Scatter Plot (run vs level)", (260, 130))
err = B.node(FC, "Formula (error)", (260, 220), {"currentIndex": 0, "descriptors": [
    ContinuousDescriptor("error", "signal - level_pct", 4, True)]})
sc2 = B.node(SCAT, "Scatter Plot (run vs error)", (460, 220))
sa = B.node(SR, "Select Rows (ascending)", (260, 340))
sb = B.node(SR, "Select Rows (blocked)", (260, 460))
la = B.node(LR, "Linear Regression (ascending)", (470, 340), {"learner_name": "LR ascending", "reg_type": 0, "fit_intercept": True})
lb = B.node(LR, "Linear Regression (blocked)", (470, 460), {"learner_name": "LR blocked", "reg_type": 0, "fit_intercept": True})
ca = B.node(TAB, "Data Table (slope, ascending)", (690, 340))
cb = B.node(TAB, "Data Table (slope, blocked)", (690, 460))

B.link(f, "Data", tab, "Data"); B.link(f, "Data", sc1, "Data"); B.link(f, "Data", err, "Data")
B.link(err, "Data", sc2, "Data")
B.link(f, "Data", sa, "Data"); B.link(f, "Data", sb, "Data")
B.link(sa, "Matching Data", la, "Data"); B.link(sb, "Matching Data", lb, "Data")
B.link(la, "Coefficients", ca, "Data"); B.link(lb, "Coefficients", cb, "Data")
B.create_widgets()

# discrete value in the Select Rows combo = 1-based index into the sorted values ("" is item 0)
vals = list(B.w("File").data.domain["order_type"].values)
B.set_conditions("Select Rows (ascending)", [("order_type", 0, [str(vals.index("由低到高") + 1)])])
B.set_conditions("Select Rows (blocked)", [("order_type", 0, [str(vals.index("分3塊打亂") + 1)])])

for title, ykey, ckey in (("Scatter Plot (run vs level)", "level_pct", "order_type"),
                          ("Scatter Plot (run vs error)", "error", None)):   # the drift is the same in both orders
    s = B.w(title)
    d = s.data.domain
    s.attr_x, s.attr_y, s.attr_color = d["run"], d[ykey], (d[ckey] if ckey else None)
    s.attr_changed() if hasattr(s, "attr_changed") else None
    s.colors_changed()
    owsutil.settle(B.model)

for t in ("Linear Regression (ascending)", "Linear Regression (blocked)"):
    print(t, "rows in:", len(B.w(t).data))
    assert len(B.w(t).data) == 24
for t in ("Data Table (slope, ascending)", "Data Table (slope, blocked)"):
    print(t, [(str(r["name"]), round(float(r[0]), 5)) for r in B.w(t).input.table])

B.save(OWS, "Ch.7 measurement order (acetic acid NIR)",
       "Illustrative design table ch07_design.csv: 24 runs, 8 levels x 3, two orders (ascending vs 3 shuffled blocks), "
       "signal = level + linear drift (0.5 over the 24 runs, illustrative). Select Rows per order -> Linear Regression "
       "-> Coefficients: slope 1.0652 vs 1.0007. Scatter plots of run vs level and run vs error show the confounding.")
