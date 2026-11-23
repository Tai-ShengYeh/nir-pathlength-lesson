"""Build exercises/orange/ch09_validation.ows with Orange's own scheme API (Orange 3.40.0), same approach as build_ch08_ows.py.
    python exercises/orange/tools/make_ch09_validation.py   (once, writes ch09_validation.csv)
    python exercises/orange/tools/build_ch09_ows.py

File (ch09_validation.csv, 4 rows) -> Data Table
  File -> Select Columns (feature d2_2258, target = one reference: label_pct / atago_1102 / titr_pct / pal30s_pct) x 4
       -> Test and Score (Leave one out)  x 4   [RMSE in % acid, R2 of the leave-one-out predictions]
       -> Test and Score (Test on train)  x 4   [R2 = calibration R2 = squared correlation, same as the lesson]
  one Linear Regression learner (No regularization, fit intercept) feeds all eight Test and Score widgets.
Orange regresses DIRECTLY: reference (%) = a + b * d2_2258, so its RMSE is already in % acid. The lesson (CH9) uses INVERSE
calibration (d2 = a + b * reference, then back-calculate), so the leave-one-out RMSE differs slightly; R2 (train) is identical."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import owsutil
from owsbuild import Builder

ORANGE_DIR = os.path.normpath(os.path.join(here, ".."))
OWS = os.path.join(ORANGE_DIR, "ch09_validation.ows")
B = Builder(ORANGE_DIR, "ch09_validation.csv")

TAB = "Orange.widgets.data.owtable.OWTable"
LR = "Orange.widgets.model.owlinearregression.OWLinearRegression"
SC = "Orange.widgets.data.owselectcolumns.OWSelectAttributes"
TAS = "Orange.widgets.evaluate.owtestandscore.OWTestAndScore"

refs = [("label", "label_pct"), ("Atago", "atago_1102"), ("titration", "titr_pct"), ("PAL-30S", "pal30s_pct")]
f = B.file_node((50, 300))
tab = B.node(TAB, "Data Table (validation)", (260, 40))
lr = B.node(LR, "Linear Regression", (480, 40), {"learner_name": "Linear Regression", "reg_type": 0, "fit_intercept": True})
B.link(f, "Data", tab, "Data")
cols, loo, trn = [], [], []
for i, (nm, col) in enumerate(refs):
    y = 160 + i * 150
    c = B.node(SC, "Select Columns (%s)" % nm, (260, y))
    a = B.node(TAS, "Test and Score LOO (%s)" % nm, (500, y - 30), {"resampling": 3})     # 3 == Leave one out
    b = B.node(TAS, "Test and Score train (%s)" % nm, (500, y + 50), {"resampling": 4})   # 4 == Test on train data
    cols.append(c); loo.append(a); trn.append(b)
    B.link(f, "Data", c, "Data")
    for t in (a, b):
        B.link(c, "Data", t, "Data"); B.link(lr, "Learner", t, "Learner")
B.create_widgets()
owsutil.settle(B.model, timeout=30, quiet=3)

for (nm, col) in refs:
    B.set_columns("Select Columns (%s)" % nm, ["d2_2258"], col, metas=[])

for (nm, col), a, b in zip(refs, loo, trn):
    for t in (a, b):
        w = B.w(t.title)
        print(t.title, owsutil.test_score_table(w))

B.save(OWS, "Ch.9 method validation evidence (acetic acid NIR)",
       "ch09_validation.csv: water + 3 vinegars with the 1 mm 2nd-derivative signal at 2258 and 1670 nm and four references "
       "(label, ATAGO PAL-Easy ACID 181, titration, PAL-30S). Select Columns puts d2_2258 as the feature and one reference as "
       "target; Linear Regression (no regularization) -> Test and Score with Leave one out (RMSE in % acid) and Test on train "
       "(R2). Orange regresses reference on signal (direct); the lesson regresses signal on reference (inverse), so only R2 "
       "matches exactly. Fieller intervals, LOD/LOQ and the scorecard are not possible in Orange.")
