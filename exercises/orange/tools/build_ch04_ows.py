"""Build exercises/orange/ch04_pca_pls.ows with Orange's own scheme API (Orange 3.40.0 +
Spectroscopy 0.9.2). Orange's scheme writer produces the file; nothing is hand-written XML.

    python exercises/orange/tools/build_ch04_ows.py

Simple settings go into node.properties before the widgets are created; context settings
(Select Rows conditions, Scatter Plot colour) are set on the live widgets after data has
arrived, exactly as a user would, and then Orange packs them when saving."""
import os, sys
here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import owsutil

ORANGE_DIR = os.path.normpath(os.path.join(here, ".."))
OWS = os.path.join(ORANGE_DIR, "ch04_pca_pls.ows")
REL = "../../data/spectra_wide_orange.csv"
CSV = os.path.normpath(os.path.join(ORANGE_DIR, REL))

reg, model = owsutil.registry_and_model(ORANGE_DIR)
from orangecanvas.scheme import SchemeNode
from orangewidget.utils.filedialogs import RecentPath


def node(qname, title, pos, props=None):
    desc = reg.widget(qname)
    props = dict(props or {})
    if props:  # stamp the current settings version, otherwise migrate_settings() treats them as v0
        import importlib
        mod, cls = qname.rsplit(".", 1)
        props["__version__"] = getattr(getattr(importlib.import_module(mod), cls), "settings_version", 1)
    n = SchemeNode(desc, title=title, position=pos, properties=props)
    model.add_node(n)
    return n


def link(a, out, b, inp):
    o = [c for c in a.output_channels() if c.name == out][0]
    i = [c for c in b.input_channels() if c.name == inp][0]
    model.new_link(a, o, b, i)


def rel_path():
    # abspath 不寫入本機路徑（公開 repo）；Orange 以 basedir + relpath 解析
    return RecentPath(REL, "basedir", REL, file_format="Orange.data.io.CSVReader")


f = node("Orange.widgets.data.owfile.OWFile", "File", (60, 300),
         {"recent_paths": [rel_path()], "source": 0})
# Branch A
selA = node("Orange.widgets.data.owselectrows.OWSelectRows", "Select Rows (replicate = 1)", (230, 150))
pcaA = node("Orange.widgets.unsupervised.owpca.OWPCA", "PCA (20 raw spectra)", (400, 150),
            {"ncomponents": 2, "normalize": False})
scat = node("Orange.widgets.visualize.owscatterplot.OWScatterPlot", "Scatter Plot", (570, 150))
# Branch B
selB = node("Orange.widgets.data.owselectrows.OWSelectRows", "Select Rows (1 mm)", (230, 450))
prep = node("orangecontrib.spectroscopy.widgets.owpreprocess.OWPreprocess", "Preprocess Spectra", (400, 450),
            {"autocommit": True, "storedsettings": {"name": "", "preprocessors": [
                ("orangecontrib.spectroscopy.savitzkygolay", {"window": 9, "polyorder": 2, "deriv": 2}),
                ("orangecontrib.spectroscopy.cut", {"lowlim": 1480.0, "highlim": 2285.0, "inverse": False}),
                ("orangecontrib.spectroscopy.cut", {"lowlim": 1875.0, "highlim": 2110.0, "inverse": True}),
            ]}})
pcaB = node("Orange.widgets.unsupervised.owpca.OWPCA", "PCA (1 mm, 2nd derivative)", (590, 360),
            {"ncomponents": 2, "normalize": False})
pls = node("Orange.widgets.model.owpls.OWPLS", "PLS", (590, 540),
           {"n_components": 1, "scale": False, "auto_apply": True})
tas = node("Orange.widgets.evaluate.owtestandscore.OWTestAndScore", "Test and Score", (770, 470),
           {"resampling": 3})  # 3 == OWTestAndScore.LeaveOneOut

link(f, "Data", selA, "Data"); link(selA, "Matching Data", pcaA, "Data"); link(pcaA, "Transformed Data", scat, "Data")
link(f, "Data", selB, "Data"); link(selB, "Matching Data", prep, "Data")
link(prep, "Preprocessed Data", pcaB, "Data"); link(prep, "Preprocessed Data", pls, "Data")
link(prep, "Preprocessed Data", tas, "Data"); link(pls, "Learner", tas, "Learner")

for n in model.nodes:
    model.widget_for_node(n)
owsutil.settle(model)


def set_conditions(title, conds):
    _, w = owsutil.by_title(model, title)
    dom = w.data.domain
    w.remove_all_rows()
    for name, op, vals in conds:
        w.add_row(dom[name], op, vals)
    w.conditions_changed()
    w.commit.now()
    owsutil.settle(model)
    print(title, "->", w.conditions, "rows out:", len(w.data) if w.data is not None else None)


set_conditions("Select Rows (replicate = 1)", [("replicate", 0, ["1"])])            # 0 == "equals"
set_conditions("Select Rows (1 mm)", [("path_mm", 0, ["1"]), ("replicate", 0, ["1"])])

_, ws = owsutil.by_title(model, "Scatter Plot")
dom = ws.data.domain
ws.attr_x, ws.attr_y = dom["PC1"], dom["PC2"]
ws.attr_color = dom["path_mm"]
ws.attr_changed() if hasattr(ws, "attr_changed") else None
ws.colors_changed()
owsutil.settle(model)

# Orange's File widget re-writes recent_paths as absolute on init (files outside the .ows
# folder are never stored relative by the GUI). Put the basedir-relative entry back.
_, wf = owsutil.by_title(model, "File")
wf.recent_paths[:] = [rel_path()]

model.title = "Ch.4 PCA and PLS (acetic acid NIR)"
model.description = ("Branch A: PCA of 20 raw spectra (replicate 1), coloured by path length. "
                     "Branch B: 1 mm, SG 9/2/d2 + Cut 1480-2285 + Cut remove 1875-2110 -> PCA and "
                     "PLS (1 component, no scaling), Leave one out. Data: " + REL)
model.save_to(OWS, pretty=True, pickle_fallback=True)  # as Orange canvas does
print("saved", OWS)
print("messages:", owsutil.messages(model))
