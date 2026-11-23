"""Build exercises/orange/ch08_references.ows with Orange's own scheme API (Orange 3.40.0), same approach as
build_ch07_ows.py.   python exercises/orange/tools/make_ch08_references.py   (once, writes ch08_references.csv)
                     python exercises/orange/tools/build_ch08_ows.py

File (ch08_references.csv, 4 rows) -> Data Table (references)
  File -> Formula (titr_pct_calc = 0.2589 * titr_mL; need_mL_label = label_pct / 0.2589; gap_mL = need - titrated)
       -> Data Table (gaps in mL), Scatter Plot x4 (d2_2258 vs label / Atago 11-02 / Atago 11-03 / titration)
  File -> Select Columns (titr_pct)      -> Linear Regression -> Coefficients -> Data Table, + Test and Score (train)
  File -> Select Columns (titr_pct_f097) -> Linear Regression -> Coefficients -> Data Table, + Test and Score (train)
Linear Regression has no regularisation: d2_2258 = a + b * reference (direct regression, 4 points incl. water)."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import owsutil
from owsbuild import Builder
from Orange.widgets.data.owfeatureconstructor import ContinuousDescriptor

ORANGE_DIR = os.path.normpath(os.path.join(here, ".."))
OWS = os.path.join(ORANGE_DIR, "ch08_references.ows")
B = Builder(ORANGE_DIR, "ch08_references.csv")

FC = "Orange.widgets.data.owfeatureconstructor.OWFeatureConstructor"
TAB = "Orange.widgets.data.owtable.OWTable"
SCAT = "Orange.widgets.visualize.owscatterplot.OWScatterPlot"
LR = "Orange.widgets.model.owlinearregression.OWLinearRegression"
SC = "Orange.widgets.data.owselectcolumns.OWSelectAttributes"
TAS = "Orange.widgets.evaluate.owtestandscore.OWTestAndScore"

f = B.file_node((50, 260))
tab = B.node(TAB, "Data Table (references)", (260, 40))
form = B.node(FC, "Formula (titration, gap)", (260, 180), {"currentIndex": 0, "descriptors": [
    ContinuousDescriptor("titr_pct_calc", "0.2589 * titr_mL", 4, True),
    ContinuousDescriptor("need_mL_label", "label_pct / 0.2589", 3, True),
    ContinuousDescriptor("gap_mL", "label_pct / 0.2589 - titr_mL", 3, True)]})
gap = B.node(TAB, "Data Table (gap in mL)", (480, 180))
scs = [("label", "label_pct", 20), ("Atago 11-02", "atago_1102", 140), ("Atago 11-03", "atago_1103", 260),
       ("titration", "titr_pct", 380)]
sc_nodes = [B.node(SCAT, "Scatter Plot (d2 vs %s)" % n, (480, 270 + 0 * y + i * 90)) for i, (n, _, y) in enumerate(scs)]
colA = B.node(SC, "Select Columns (titration)", (260, 660))
colB = B.node(SC, "Select Columns (titration x 0.97)", (260, 800))
lrA = B.node(LR, "Linear Regression (titration)", (480, 660), {"learner_name": "LR titration", "reg_type": 0, "fit_intercept": True})
lrB = B.node(LR, "Linear Regression (titration x 0.97)", (480, 800), {"learner_name": "LR titration x 0.97", "reg_type": 0, "fit_intercept": True})
cfA = B.node(TAB, "Data Table (slope, titration)", (720, 620))
cfB = B.node(TAB, "Data Table (slope, titration x 0.97)", (720, 760))
tsA = B.node(TAS, "Test and Score (titration)", (720, 700), {"resampling": 4})   # 4 == Test on train data
tsB = B.node(TAS, "Test and Score (titration x 0.97)", (720, 840), {"resampling": 4})

B.link(f, "Data", tab, "Data"); B.link(f, "Data", form, "Data"); B.link(form, "Data", gap, "Data")
for s in sc_nodes: B.link(form, "Data", s, "Data")
B.link(f, "Data", colA, "Data"); B.link(f, "Data", colB, "Data")
B.link(colA, "Data", lrA, "Data"); B.link(colB, "Data", lrB, "Data")
B.link(lrA, "Coefficients", cfA, "Data"); B.link(lrB, "Coefficients", cfB, "Data")
B.link(colA, "Data", tsA, "Data"); B.link(lrA, "Learner", tsA, "Learner")
B.link(colB, "Data", tsB, "Data"); B.link(lrB, "Learner", tsB, "Learner")
B.create_widgets()
owsutil.settle(B.model, timeout=30, quiet=3)   # Formula runs in a thread

B.set_columns("Select Columns (titration)", ["titr_pct"], "d2_2258", metas=[])
B.set_columns("Select Columns (titration x 0.97)", ["titr_pct_f097"], "d2_2258", metas=[])

for (n, col, _), node in zip(scs, sc_nodes):
    s = B.w(node.title)
    d = s.data.domain
    s.attr_x, s.attr_y = d[col], d["d2_2258"]
    s.attr_label = d["bottle"]
    s.attr_changed() if hasattr(s, "attr_changed") else None
    s.graph.show_reg_line = True
    s.graph.update_regression_line() if hasattr(s.graph, "update_regression_line") else None
    owsutil.settle(B.model)

for t in ("Data Table (slope, titration)", "Data Table (slope, titration x 0.97)"):
    print(t, [(str(r["name"]), round(float(r[0]), 5)) for r in B.w(t).input.table])
for t in ("Test and Score (titration)", "Test and Score (titration x 0.97)"):
    print(t, owsutil.test_score_table(B.w(t)))

B.save(OWS, "Ch.8 reference-method errors (acetic acid NIR)",
       "ch08_references.csv: water + 3 vinegars with label, two Atago readings, NaOH titration volume and the 1 mm "
       "2258 nm second-derivative signal. Formula: titration acidity 0.2589 x V, mL needed for the label value and the gap "
       "in mL. Scatter plots of the signal vs each reference. Linear Regression (no regularisation) on titration vs "
       "titration x 0.97: slope changes by 1/0.97, R2 identical.")
