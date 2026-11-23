"""Build exercises/orange/ch03_wavelength_pitfall.ows with Orange's own scheme API (Orange 3.40.0 +
Spectroscopy 0.9.2), same approach as build_ch04_ows.py.

    python exercises/orange/tools/build_ch03_ows.py

File -> Select Rows (1 mm) -> Preprocess Spectra (SG 9/2/d2; Cut keep 1340-2290; Cut remove
1415-1477; Cut remove 1875-2110 = the 88 wavelengths where all four 1 mm spectra have A < 1.5)
  (A) selection OUTSIDE the CV loop: Rank (Univariate Regression, best 1) -> Test and Score (A)
  (B) selection INSIDE the CV loop:  Preprocess (Select Relevant Features, Univariate Linear
      Regression, fixed 1) -> Test and Score (B) "Preprocessor" input; data = all 88 wavelengths
  Linear Regression is the learner of both Test and Score widgets (Leave one out).
  Optional: Randomize (shuffle the target, replicable) -> Rank (shuffled target)."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import owsutil
from owsbuild import Builder

ORANGE_DIR = os.path.normpath(os.path.join(here, ".."))
OWS = os.path.join(ORANGE_DIR, "ch03_wavelength_pitfall.ows")
B = Builder(ORANGE_DIR, "../../data/spectra_wide_orange.csv")
RANK = "Orange.widgets.data.owrank.OWRank"
# Rank keeps one set of method names for classification and regression; this is the default set
# with RReliefF unticked (what a user gets after unchecking it in the GUI).
RANK_PROPS = {"selected_methods": {"Information Gain Ratio", "Gini Decrease", "Univariate Regression"},
              "sorting": (0, 1), "auto_apply": True}   # column 0 = Univar. reg., 1 = descending

f = B.file_node((40, 300))
sel = B.node("Orange.widgets.data.owselectrows.OWSelectRows", "Select Rows (1 mm)", (180, 300))
prep = B.node("orangecontrib.spectroscopy.widgets.owpreprocess.OWPreprocess", "Preprocess Spectra (d2, 88 wl)",
              (330, 300),
              {"autocommit": True, "storedsettings": {"name": "", "preprocessors": [
                  ("orangecontrib.spectroscopy.savitzkygolay", {"window": 9, "polyorder": 2, "deriv": 2}),
                  ("orangecontrib.spectroscopy.cut", {"lowlim": 1340.0, "highlim": 2290.0, "inverse": False}),
                  ("orangecontrib.spectroscopy.cut", {"lowlim": 1415.0, "highlim": 1477.0, "inverse": True}),
                  ("orangecontrib.spectroscopy.cut", {"lowlim": 1875.0, "highlim": 2110.0, "inverse": True}),
              ]}})
rank = B.node(RANK, "Rank (outside CV)", (500, 160), RANK_PROPS)
lr = B.node("Orange.widgets.model.owlinearregression.OWLinearRegression", "Linear Regression", (500, 300),
            {"reg_type": 0, "fit_intercept": True, "auto_apply": True})
tsA = B.node("Orange.widgets.evaluate.owtestandscore.OWTestAndScore", "Test and Score (A: outside)", (680, 160),
             {"resampling": 3})  # 3 == OWTestAndScore.LeaveOneOut
fss = B.node("Orange.widgets.data.owpreprocess.OWPreprocess", "Preprocess (select 1 inside CV)", (500, 440),
             {"autocommit": True, "storedsettings": {"name": "", "preprocessors": [
                 ("orange.preprocess.fss", {"score": 8, "strategy": 1, "k": 1, "p": 75.0})]}})
tsB = B.node("Orange.widgets.evaluate.owtestandscore.OWTestAndScore", "Test and Score (B: inside)", (680, 440),
             {"resampling": 3})
rnd = B.node("Orange.widgets.data.owrandomize.OWRandomize", "Randomize (shuffle target)", (500, 590),
             {"shuffle_class": True, "shuffle_attrs": False, "shuffle_metas": False, "scope_prop": 100,
              "random_seed": True, "auto_apply": True})
rankS = B.node(RANK, "Rank (shuffled target)", (680, 590), RANK_PROPS)

B.link(f, "Data", sel, "Data"); B.link(sel, "Matching Data", prep, "Data")
B.link(prep, "Preprocessed Data", rank, "Data")
B.link(rank, "Reduced Data", tsA, "Data"); B.link(lr, "Learner", tsA, "Learner")
B.link(prep, "Preprocessed Data", tsB, "Data"); B.link(lr, "Learner", tsB, "Learner")
B.link(fss, "Preprocessor", tsB, "Preprocessor")
B.link(prep, "Preprocessed Data", rnd, "Data"); B.link(rnd, "Data", rankS, "Data")
B.create_widgets()

B.set_conditions("Select Rows (1 mm)", [("path_mm", 0, ["1"]), ("replicate", 0, ["1"])])   # 0 == "equals"
for t in ["Rank (outside CV)", "Rank (shuffled target)"]:
    w = B.w(t)
    w.nSelected = 1
    w.setSelectionMethod(w.SelectNBest) if hasattr(w, "setSelectionMethod") else None
    w.autoSelection() if hasattr(w, "autoSelection") else None
    w.commit.now() if hasattr(w.commit, "now") else w.commit()
owsutil.settle(B.model)
for t in ["Rank (outside CV)", "Rank (shuffled target)"]:
    w = B.w(t)
    print(t, "methods", w.selected_methods, "sorting", w.sorting, "nSelected", w.nSelected,
          "selectionMethod", w.selectionMethod)

B.save(OWS, "Ch.3 wavelength selection pitfall (acetic acid NIR)",
       "1 mm, SG 9/2/d2, 88 wavelengths with A < 1.5. (A) Rank picks the best wavelength with all 4 "
       "samples, then Leave one out; (B) Preprocess 'Select Relevant Features' on the Preprocessor input "
       "of Test and Score re-selects inside every fold. Data: ../../data/spectra_wide_orange.csv")
