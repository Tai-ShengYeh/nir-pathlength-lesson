"""Build exercises/orange/ch02_references.ows with Orange's own scheme API (Orange 3.40.0 +
Spectroscopy 0.9.2), same approach as build_ch04_ows.py.

    python tools/nir-pathlength/make_orange_refs.py      # data/spectra_1mm_refs_orange.csv
    python exercises/orange/tools/build_ch02_ows.py

File (1 mm, 4 rows, references as metas) -> Preprocess Spectra (SG 9/2/d2)
  -> Formula (titration = naoh_ml * 1.07875 * 0.24; ratio = titration / atago_1103; both metas)
  -> Select Columns (feature 2258.067 nm, target titration) -> Linear Regression
  -> Test and Score (Test on train data -> R2 of the calibration line)
  Formula -> Scatter Plot (atago_1103 vs titration);  Formula -> Data Table (see the ratio)"""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import owsutil
from owsbuild import Builder
from Orange.widgets.data.owfeatureconstructor import ContinuousDescriptor

ORANGE_DIR = os.path.normpath(os.path.join(here, ".."))
OWS = os.path.join(ORANGE_DIR, "ch02_references.ows")
B = Builder(ORANGE_DIR, "../../data/spectra_1mm_refs_orange.csv")

f = B.file_node((40, 260))
prep = B.node("orangecontrib.spectroscopy.widgets.owpreprocess.OWPreprocess", "Preprocess Spectra", (190, 260),
              {"autocommit": True, "storedsettings": {"name": "", "preprocessors": [
                  ("orangecontrib.spectroscopy.savitzkygolay", {"window": 9, "polyorder": 2, "deriv": 2})]}})
form = B.node("Orange.widgets.data.owfeatureconstructor.OWFeatureConstructor", "Formula (titration, ratio)",
              (340, 260),
              {"descriptors": [ContinuousDescriptor("titration", "naoh_ml * 1.07875 * 0.24", 3, True),
                               ContinuousDescriptor("ratio", "naoh_ml * 1.07875 * 0.24 / atago_1103 if atago_1103 > 0 else nan", 3, True)],
               "currentIndex": 0})
cols = B.node("Orange.widgets.data.owselectcolumns.OWSelectAttributes", "Select Columns (2258 nm)", (500, 260))
lr = B.node("Orange.widgets.model.owlinearregression.OWLinearRegression", "Linear Regression", (660, 160),
            {"reg_type": 0, "fit_intercept": True, "auto_apply": True})
tas = B.node("Orange.widgets.evaluate.owtestandscore.OWTestAndScore", "Test and Score", (820, 260),
             {"resampling": 4})  # 4 == OWTestAndScore.TestOnTrain ("Test on train data")
scat = B.node("Orange.widgets.visualize.owscatterplot.OWScatterPlot", "Scatter Plot", (500, 420))
table = B.node("Orange.widgets.data.owtable.OWTable", "Data Table", (500, 110))

B.link(f, "Data", prep, "Data"); B.link(prep, "Preprocessed Data", form, "Data")
B.link(form, "Data", cols, "Data")
B.link(cols, "Data", lr, "Data"); B.link(cols, "Data", tas, "Data"); B.link(lr, "Learner", tas, "Learner")
B.link(form, "Data", scat, "Data"); B.link(form, "Data", table, "Data")
B.create_widgets()

owsutil.settle(B.model, timeout=30, quiet=3)   # Formula runs in a thread; let it finish
B.set_columns("Select Columns (2258 nm)", ["2258.067"], "titration")

sc = B.w("Scatter Plot")
dom = sc.data.domain
sc.attr_x, sc.attr_y = dom["atago_1103"], dom["titration"]
sc.attr_label = dom["sample"]
sc.attr_changed() if hasattr(sc, "attr_changed") else None
sc.graph.show_reg_line = True
sc.graph.update_regression_line() if hasattr(sc.graph, "update_regression_line") else None
owsutil.settle(B.model)

B.save(OWS, "Ch.2 reference values (acetic acid NIR)",
       "1 mm spectra of water + 3 vinegars with four reference values. Formula computes the titration "
       "acidity (V x 1.07875 x 0.24) and its ratio to Atago 11/03; Select Columns picks 2258.067 nm and the "
       "target; Test and Score (Test on train data) gives the R2 of each reference. "
       "Data: ../../data/spectra_1mm_refs_orange.csv")
