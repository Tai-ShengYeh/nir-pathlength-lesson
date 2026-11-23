"""Build exercises/orange/ch01_pathlength.ows with Orange's own scheme API (Orange 3.40.0 +
Spectroscopy 0.9.2), same approach as build_ch04_ows.py.

    python exercises/orange/tools/build_ch01_ows.py

File -> Formula (categorical copy "path" of path_mm) -> Spectra (all 21 spectra, coloured by path)
File -> Select Rows (path_mm = 1, replicate = 1) -> Preprocess Spectra (SG 9/2/d2)
     -> Select Columns (feature 2258.067 nm, target acid_pct) -> Linear Regression
     -> Test and Score (Leave one out);  Select Columns -> Scatter Plot (d2 vs acid_pct)"""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import owsutil
from owsbuild import Builder

ORANGE_DIR = os.path.normpath(os.path.join(here, ".."))
OWS = os.path.join(ORANGE_DIR, "ch01_pathlength.ows")
B = Builder(ORANGE_DIR, "../../data/spectra_wide_orange.csv")

f = B.file_node((60, 260))
from Orange.widgets.data.owfeatureconstructor import DiscreteDescriptor
# Spectra can only colour by a categorical variable; path_mm is numeric -> make a categorical copy
form = B.node("Orange.widgets.data.owfeatureconstructor.OWFeatureConstructor", "Formula (path label)", (200, 90),
              {"descriptors": [DiscreteDescriptor("path", "str(path_mm) + ' mm'", (), False, True)],
               "currentIndex": 0})
spec = B.node("orangecontrib.spectroscopy.widgets.owspectra.OWSpectra", "Spectra (all 21)", (370, 90))
sel = B.node("Orange.widgets.data.owselectrows.OWSelectRows", "Select Rows (1 mm)", (240, 300))
prep = B.node("orangecontrib.spectroscopy.widgets.owpreprocess.OWPreprocess", "Preprocess Spectra", (410, 300),
              {"autocommit": True, "storedsettings": {"name": "", "preprocessors": [
                  ("orangecontrib.spectroscopy.savitzkygolay", {"window": 9, "polyorder": 2, "deriv": 2})]}})
cols = B.node("Orange.widgets.data.owselectcolumns.OWSelectAttributes", "Select Columns (2258 nm)", (580, 300))
lr = B.node("Orange.widgets.model.owlinearregression.OWLinearRegression", "Linear Regression", (750, 200),
            {"reg_type": 0, "fit_intercept": True, "auto_apply": True})
tas = B.node("Orange.widgets.evaluate.owtestandscore.OWTestAndScore", "Test and Score", (920, 300),
             {"resampling": 3})  # 3 == OWTestAndScore.LeaveOneOut
scat = B.node("Orange.widgets.visualize.owscatterplot.OWScatterPlot", "Scatter Plot", (750, 420))

B.link(f, "Data", form, "Data"); B.link(form, "Data", spec, "Data")
B.link(f, "Data", sel, "Data"); B.link(sel, "Matching Data", prep, "Data")
B.link(prep, "Preprocessed Data", cols, "Data")
B.link(cols, "Data", lr, "Data"); B.link(cols, "Data", tas, "Data"); B.link(lr, "Learner", tas, "Learner")
B.link(cols, "Data", scat, "Data")
B.create_widgets()

B.set_conditions("Select Rows (1 mm)", [("path_mm", 0, ["1"]), ("replicate", 0, ["1"])])   # 0 == "equals"
B.set_columns("Select Columns (2258 nm)", ["2258.067"], "acid_pct")

ws = B.w("Spectra (all 21)")
ws.curveplot.feature_color = ws.curveplot.data.domain["path"]
ws.curveplot.update_view() if hasattr(ws.curveplot, "update_view") else None
owsutil.settle(B.model)

sc = B.w("Scatter Plot")
dom = sc.data.domain
sc.attr_x, sc.attr_y = dom["acid_pct"], dom["2258.067"]
sc.attr_changed() if hasattr(sc, "attr_changed") else None
sc.graph.show_reg_line = True
sc.graph.update_regression_line() if hasattr(sc.graph, "update_regression_line") else None
owsutil.settle(B.model)

B.save(OWS, "Ch.1 path length (acetic acid NIR)",
       "Spectra of all 21 rows (saturation plateaus); 1 mm rows -> SG 9/2/d2 -> 2258.067 nm only -> "
       "Linear Regression, Test and Score Leave one out. Change Select Rows to another path, or "
       "Select Columns to 1669.901 nm. Data: ../../data/spectra_wide_orange.csv")
