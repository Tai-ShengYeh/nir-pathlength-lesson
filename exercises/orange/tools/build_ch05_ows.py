"""Build exercises/orange/ch05_prediction.ows with Orange's own scheme API (Orange 3.40.0 +
Spectroscopy 0.9.2), same approach as build_ch01_ows.py.

    python tools/nir-pathlength/make_orange_unknown.py      # data/unknown_2258_orange.csv (branch B)
    python exercises/orange/tools/build_ch05_ows.py

Branch A (pretend unknown = the 4.87 % bottle):
  File -> Select Rows (1 mm, replicate 1) -> Preprocess Spectra (SG 9/2/d2) -> Select Columns (2258.067)
       -> Select Rows (sample = vinegar_5%):  Unmatched Data (3 rows) -> Linear Regression (3 bottles)
                                              Matching Data (1 row)  -> Predictions (4.87 % bottle)
Branch B (4-point calibration, hypothetical unknown signal 4.20 x 10^-3):
  File (unknown_2258) -> Select Rows (set = unknown): Unmatched Data (4 rows) -> Linear Regression (4 bottles)
                                                      Matching Data (1 row)  -> Predictions (y0 = 4.20)
  Linear Regression (4 bottles) -> Coefficients -> Data Table (coefficients)
Orange's Linear Regression is DIRECT (conc = a + b*signal) and gives a point prediction only."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import owsutil
from owsbuild import Builder
from orangewidget.utils.filedialogs import RecentPath

ORANGE_DIR = os.path.normpath(os.path.join(here, ".."))
OWS = os.path.join(ORANGE_DIR, "ch05_prediction.ows")
REL_U = "../../data/unknown_2258_orange.csv"
B = Builder(ORANGE_DIR, "../../data/spectra_wide_orange.csv")


def rel_u():
    # abspath 不寫入本機路徑（公開 repo）；Orange 以 basedir + relpath 解析
    return RecentPath(REL_U, "basedir", REL_U, file_format="Orange.data.io.CSVReader")


SR = "Orange.widgets.data.owselectrows.OWSelectRows"
LR = "Orange.widgets.model.owlinearregression.OWLinearRegression"
PR = "Orange.widgets.evaluate.owpredictions.OWPredictions"

# ---- branch A ----------------------------------------------------------------------------
f = B.file_node((60, 160))
sel = B.node(SR, "Select Rows (1 mm)", (220, 160))
prep = B.node("orangecontrib.spectroscopy.widgets.owpreprocess.OWPreprocess", "Preprocess Spectra", (380, 160),
              {"autocommit": True, "storedsettings": {"name": "", "preprocessors": [
                  ("orangecontrib.spectroscopy.savitzkygolay", {"window": 9, "polyorder": 2, "deriv": 2})]}})
cols = B.node("Orange.widgets.data.owselectcolumns.OWSelectAttributes", "Select Columns (2258 nm)", (540, 160))
hide = B.node(SR, "Select Rows (hide 4.87 %)", (700, 160))
lr3 = B.node(LR, "Linear Regression (3 bottles)", (870, 80),
             {"learner_name": "LR (3 bottles)", "reg_type": 0, "fit_intercept": True})
pr3 = B.node(PR, "Predictions (4.87 % bottle)", (1040, 160))

# ---- branch B ----------------------------------------------------------------------------
fu = B.node("Orange.widgets.data.owfile.OWFile", "File (unknown 2258)", (60, 420),
            {"recent_paths": [rel_u()], "source": 0})
selu = B.node(SR, "Select Rows (set = unknown)", (300, 420))
lr4 = B.node(LR, "Linear Regression (4 bottles)", (560, 340),
             {"learner_name": "LR (4 bottles)", "reg_type": 0, "fit_intercept": True})
pr4 = B.node(PR, "Predictions (y0 = 4.20)", (800, 420))
coef = B.node("Orange.widgets.data.owtable.OWTable", "Data Table (coefficients)", (800, 300))

B.link(f, "Data", sel, "Data"); B.link(sel, "Matching Data", prep, "Data")
B.link(prep, "Preprocessed Data", cols, "Data"); B.link(cols, "Data", hide, "Data")
B.link(hide, "Unmatched Data", lr3, "Data"); B.link(hide, "Matching Data", pr3, "Data")
B.link(lr3, "Model", pr3, "Predictors")
B.link(fu, "Data", selu, "Data")
B.link(selu, "Unmatched Data", lr4, "Data"); B.link(selu, "Matching Data", pr4, "Data")
B.link(lr4, "Model", pr4, "Predictors"); B.link(lr4, "Coefficients", coef, "Data")
B.create_widgets()

B.set_conditions("Select Rows (1 mm)", [("path_mm", 0, ["1"]), ("replicate", 0, ["1"])])   # 0 == "equals"
B.set_columns("Select Columns (2258 nm)", ["2258.067"], "acid_pct")
B.set_conditions("Select Rows (hide 4.87 %)", [("sample", 0, ["vinegar_5%"])])              # string 0 == "equals"
B.set_conditions("Select Rows (set = unknown)", [("set", 0, ["unknown"])])


# set_conditions prints the INPUT size; check the real output sizes at the downstream widgets
n3 = len(B.w("Linear Regression (3 bottles)").data); n1 = len(B.w("Predictions (4.87 % bottle)").data)
n4 = len(B.w("Linear Regression (4 bottles)").data); nu = len(B.w("Predictions (y0 = 4.20)").data)
print("rows reaching LR(3) / Pred(4.87) / LR(4) / Pred(y0):", n3, n1, n4, nu)
assert (n3, n1, n4, nu) == (3, 1, 4, 1)

# Orange's File widget re-writes recent_paths as absolute on init; put the relative entry back
B.w("File (unknown 2258)").recent_paths[:] = [rel_u()]
B.save(OWS, "Ch.5 predicting an unknown (acetic acid NIR)",
       "A: 1 mm -> SG 9/2/d2 -> 2258.067 nm only; hide the 4.87 % bottle, Linear Regression on the other 3, "
       "Predictions for the hidden bottle. B: 4 calibration signals (x1e3) + one unknown y0 = 4.20 from "
       "unknown_2258_orange.csv, Linear Regression on 4, Predictions for the unknown. Orange's regression is "
       "direct (conc = a + b*signal) and gives a point prediction only; the lesson's inverse calibration "
       "and Fieller interval are in R / Python. Data: ../../data/spectra_wide_orange.csv, ../../data/unknown_2258_orange.csv")
