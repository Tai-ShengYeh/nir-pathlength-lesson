"""Build exercises/orange/ch06_data_quality.ows with Orange's own scheme API (Orange 3.40.0 +
Spectroscopy 0.9.2), same approach as build_ch01_ows.py.

    python exercises/orange/tools/build_ch06_ows.py

All branches start from File -> Select Rows (replicate = 1)  (20 spectra; the dataset has one extra
replicate: water at 0.25 mm).

Branch W (water-band check, lesson point 1):
  Integrate Spectra (Integral from 0 over 1440-1460 and 1580-1650, Peak from 0 over 1440-1460)
  -> Formula (height = I1/15.637 - I2/64.717, peak) -> Select Columns (path_mm, height)
  -> Group by path_mm (Median of height) -> Merge Data (path_mm = path_mm, with the Formula output)
  -> Formula (ratio = height / median, flag = usable AND |ratio - 1| > 0.10) -> Select Columns (summary) -> Data Table, Select Rows (flag = 1), Scatter Plot
Branch V (visual check, lesson point 2): Select Rows (0.25 mm) -> Formula (bottle) -> Spectra
Branch O (automatic-detection trap, lesson point 3): Outliers (defaults) -> Data Table;  PCA -> Scatter Plot

Orange's Integral is a trapezoid sum over the measured points inside the window, so dividing by the
spanned width gives the window mean up to a ~0.0002 difference in the ratios (see verify_ch06_ows.py)."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import owsutil
from owsbuild import Builder
from Orange.widgets.data.owfeatureconstructor import ContinuousDescriptor, DiscreteDescriptor

ORANGE_DIR = os.path.normpath(os.path.join(here, ".."))
OWS = os.path.join(ORANGE_DIR, "ch06_data_quality.ows")
B = Builder(ORANGE_DIR, "../../data/spectra_wide_orange.csv")

SR = "Orange.widgets.data.owselectrows.OWSelectRows"
FC = "Orange.widgets.data.owfeatureconstructor.OWFeatureConstructor"
TAB = "Orange.widgets.data.owtable.OWTable"
SCAT = "Orange.widgets.visualize.owscatterplot.OWScatterPlot"

# widths actually spanned by the measured points inside the two windows (nm)
W_WATER, W_FLAT = 15.637, 64.717     # 1443.165..1458.802 and 1584.725..1649.442

f = B.file_node((60, 300))
r1 = B.node(SR, "Select Rows (replicate = 1)", (220, 300))

# ---- branch W ----------------------------------------------------------------------------
integ = B.node("orangecontrib.spectroscopy.widgets.owintegrate.OWIntegrate", "Integrate Spectra", (400, 60),
               {"autocommit": True, "output_metas": True, "storedsettings": {"name": "", "preprocessors": [
                   ("orangecontrib.infrared.integrate.simple", {"Low limit": 1440.0, "High limit": 1460.0}),
                   ("orangecontrib.infrared.integrate.simple", {"Low limit": 1580.0, "High limit": 1650.0}),
                   ("orangecontrib.infrared.integrate.peak_max", {"Low limit": 1440.0, "High limit": 1460.0})]}})
hgt = B.node(FC, "Formula (water-band height)", (580, 60), {"currentIndex": 0, "descriptors": [
    ContinuousDescriptor("height", "_1440_0___1460_0 / %s - _1580_0___1650_0 / %s" % (W_WATER, W_FLAT), 4, True),
    ContinuousDescriptor("peak", "_1440_0___1460_0__1_", 4, True)]})
selc = B.node("Orange.widgets.data.owselectcolumns.OWSelectAttributes", "Select Columns (path, height)", (760, 140))
gb = B.node("Orange.widgets.data.owgroupby.OWGroupBy", "Group by path", (940, 140))
mg = B.node("Orange.widgets.data.owmergedata.OWMergeData", "Merge Data", (1120, 60))
rat = B.node(FC, "Formula (ratio, flag)", (1300, 60), {"currentIndex": 0, "descriptors": [
    ContinuousDescriptor("ratio", "height / height___Median", 4, True),
    ContinuousDescriptor("flag", "1 if (peak < 1.5 and abs(height / height___Median - 1) > 0.10) else 0", 0, True)]})
summ = B.node("Orange.widgets.data.owselectcolumns.OWSelectAttributes", "Select Columns (summary)", (1480, 60))
tab = B.node(TAB, "Data Table (ratios)", (1660, 0))
fl = B.node(SR, "Select Rows (flag = 1)", (1660, 120))
tabf = B.node(TAB, "Data Table (flagged)", (1840, 120))
scr = B.node(SCAT, "Scatter Plot (ratio vs path)", (1660, 240))

# ---- branch V ----------------------------------------------------------------------------
s25 = B.node(SR, "Select Rows (0.25 mm)", (400, 300))
lab = B.node(FC, "Formula (bottle label)", (580, 300), {"currentIndex": 0, "descriptors": [
    DiscreteDescriptor("bottle", "sample", (), False, True)]})
spec = B.node("orangecontrib.spectroscopy.widgets.owspectra.OWSpectra", "Spectra (0.25 mm)", (760, 300))

# ---- branch O ----------------------------------------------------------------------------
out = B.node("Orange.widgets.data.owoutliers.OWOutliers", "Outliers", (400, 480))
tabo = B.node(TAB, "Data Table (Outliers)", (580, 480))
pca = B.node("Orange.widgets.unsupervised.owpca.OWPCA", "PCA (20 raw spectra)", (400, 620),
             {"ncomponents": 2, "normalize": False})
scp = B.node(SCAT, "Scatter Plot (PCA)", (580, 620))

B.link(f, "Data", r1, "Data")
B.link(r1, "Matching Data", integ, "Data"); B.link(integ, "Integrated Data", hgt, "Data")
B.link(hgt, "Data", selc, "Data"); B.link(selc, "Data", gb, "Data")
B.link(hgt, "Data", mg, "Data"); B.link(gb, "Data", mg, "Extra Data")
B.link(mg, "Data", rat, "Data"); B.link(rat, "Data", summ, "Data"); B.link(summ, "Data", tab, "Data")
B.link(summ, "Data", fl, "Data"); B.link(fl, "Matching Data", tabf, "Data")
B.link(summ, "Data", scr, "Data")
B.link(r1, "Matching Data", s25, "Data"); B.link(s25, "Matching Data", lab, "Data"); B.link(lab, "Data", spec, "Data")
B.link(r1, "Matching Data", out, "Data"); B.link(out, "Data", tabo, "Data")
B.link(r1, "Matching Data", pca, "Data"); B.link(pca, "Transformed Data", scp, "Data")
B.create_widgets()

# ---- settings a user would click, applied to the live widgets -------------------------------
B.set_conditions("Select Rows (replicate = 1)", [("replicate", 0, ["1"])])        # 0 == "equals"
B.set_conditions("Select Rows (0.25 mm)", [("path_mm", 0, ["0.25"])])
B.set_columns("Select Columns (path, height)", [], None, metas=["path_mm", "height"])

g = B.w("Group by path")
dom = g.data.domain
assert [v.name for v in g.gb_attrs] == ["path_mm"], g.gb_attrs
g.aggregations[dom["height"]] = {"Median"}
g.aggregations[dom["path_mm"]] = set()               # path_mm is the grouping column itself
for v in (dom["height"], dom["path_mm"]):
    g.agg_table_model.update_aggregation(v)
g.commit.now()
owsutil.settle(B.model)

m = B.w("Merge Data")
m.attr_boxes.set_state([(m.data.domain["path_mm"], m.extra_data.domain["path_mm"])])
m.store_combo_state()
m.commit.now()
owsutil.settle(B.model)

B.set_columns("Select Columns (summary)", [], "acid_pct",
              metas=["path_mm", "sample", "height", "height - Median", "ratio", "flag", "peak"])
B.set_conditions("Select Rows (flag = 1)", [("flag", 0, ["1"])])

sc = B.w("Scatter Plot (ratio vs path)")
d = sc.data.domain
sc.attr_x, sc.attr_y, sc.attr_color = d["path_mm"], d["ratio"], d["flag"]
sc.attr_changed() if hasattr(sc, "attr_changed") else None
sc.colors_changed()
owsutil.settle(B.model)

sp = B.w("Spectra (0.25 mm)")
sp.curveplot.feature_color = sp.curveplot.data.domain["bottle"]
sp.curveplot.update_view() if hasattr(sp.curveplot, "update_view") else None
owsutil.settle(B.model)

sp2 = B.w("Scatter Plot (PCA)")
d2 = sp2.data.domain
sp2.attr_x, sp2.attr_y, sp2.attr_color = d2["PC1"], d2["PC2"], d2["path_mm"]
sp2.attr_changed() if hasattr(sp2, "attr_changed") else None
sp2.colors_changed()
owsutil.settle(B.model)

# report the numbers that arrive at the tables, as a sanity check before saving
t = B.w("Data Table (ratios)").input.table
print("ratios table rows:", len(t))
assert len(t) == 20 and len(B.w("Data Table (flagged)").input.table) == 2
print("flagged:", [(float(r["path_mm"]), str(r["sample"]), round(float(r["ratio"]), 4))
                   for r in B.w("Data Table (flagged)").input.table])

B.save(OWS, "Ch.6 data quality and outliers (acetic acid NIR)",
       "20 spectra (replicate 1). W: Integrate Spectra -> Formula water-band height -> Group by path (median) "
       "-> Merge Data -> Formula ratio and flag (usable and |ratio-1| > 0.10) -> Select Columns (summary) -> tables and scatter. "
       "V: 0.25 mm spectra. O: Outliers (defaults) and PCA on all 20 raw spectra (path lengths mixed). "
       "Data: ../../data/spectra_wide_orange.csv")
