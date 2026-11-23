"""Helpers shared by build_ch01/02/03_ows.py (same approach as build_ch04_ows.py): create nodes
through Orange's scheme API, link them, set context settings on the live widgets, keep the
File widget path basedir-relative, save with Orange's own writer, and check the saved file
for local absolute paths."""
import os, re, base64, pickle, importlib
import owsutil


class Builder:
    def __init__(self, orange_dir, rel_csv):
        self.dir = orange_dir
        self.rel = rel_csv
        self.reg, self.model = owsutil.registry_and_model(orange_dir)

    # ---- nodes / links -------------------------------------------------------------
    def node(self, qname, title, pos, props=None):
        from orangecanvas.scheme import SchemeNode
        desc = self.reg.widget(qname)
        props = dict(props or {})
        if props:  # stamp the current settings version, otherwise migrate_settings() treats them as v0
            mod, cls = qname.rsplit(".", 1)
            props["__version__"] = getattr(getattr(importlib.import_module(mod), cls), "settings_version", 1)
        n = SchemeNode(desc, title=title, position=pos, properties=props)
        self.model.add_node(n)
        return n

    def link(self, a, out, b, inp):
        o = [c for c in a.output_channels() if c.name == out][0]
        i = [c for c in b.input_channels() if c.name == inp][0]
        self.model.new_link(a, o, b, i)

    def rel_path(self):
        # abspath 不寫入本機路徑（公開 repo）；Orange 以 basedir + relpath 解析
        from orangewidget.utils.filedialogs import RecentPath
        return RecentPath(self.rel, "basedir", self.rel, file_format="Orange.data.io.CSVReader")

    def file_node(self, pos=(60, 300)):
        return self.node("Orange.widgets.data.owfile.OWFile", "File", pos,
                         {"recent_paths": [self.rel_path()], "source": 0})

    def create_widgets(self):
        for n in self.model.nodes:
            self.model.widget_for_node(n)
        owsutil.settle(self.model)

    def w(self, title):
        return owsutil.by_title(self.model, title)[1]

    # ---- context settings set on live widgets, as a user would ---------------------
    def set_conditions(self, title, conds):
        w = self.w(title)
        dom = w.data.domain
        w.remove_all_rows()
        for name, op, vals in conds:
            w.add_row(dom[name], op, vals)
        w.conditions_changed()
        w.commit.now()
        owsutil.settle(self.model)
        print(title, "->", w.conditions, "rows out:", None if w.data is None else len(w.data))

    def set_columns(self, title, features, target, metas=None):
        """Select Columns: `features` -> Features, `target` -> Target, `metas` (None = all
        remaining metas of the input) -> Metas, everything else -> Ignored."""
        w = self.w(title)
        dom = w.data.domain
        allv = list(dom.attributes) + list(dom.class_vars) + list(dom.metas)
        feats = [dom[f] for f in features]
        tgt = [dom[target]] if target else []
        if metas is None:
            mets = [v for v in dom.metas if v not in feats + tgt]
        else:
            mets = [dom[m] for m in metas]
        used = set(feats + tgt + mets)
        w.used_attrs[:] = feats
        w.class_attrs[:] = tgt
        w.meta_attrs[:] = mets
        w.available_attrs[:] = [v for v in allv if v not in used]
        w.commit.now() if hasattr(w.commit, "now") else w.commit()
        owsutil.settle(self.model)
        out = w.data
        print(title, "-> features", [v.name for v in feats], "target", target,
              "metas", [v.name for v in mets])

    # ---- save ----------------------------------------------------------------------
    def save(self, ows, title, description):
        # Orange's File widget re-writes recent_paths as absolute on init (files outside the
        # .ows folder are never stored relative by the GUI). Put the basedir-relative entry back.
        wf = self.w("File")
        wf.recent_paths[:] = [self.rel_path()]
        self.model.title = title
        self.model.description = description
        self.model.save_to(ows, pretty=True, pickle_fallback=True)  # as Orange canvas does
        print("saved", ows)
        print("messages:", owsutil.messages(self.model))
        leak_check(ows)


BS = chr(92)
BAD = ("P:" + BS, "P:/", "p:" + BS, "acetic_acid_", "Users", "AppData", BS + BS)


def _deep(obj, depth=0):
    """Readable dump of a decoded setting, descending into Context objects (whose repr hides values)."""
    if depth > 12:
        return ""
    if hasattr(obj, "__dict__") and not isinstance(obj, type):
        return type(obj).__name__ + _deep(vars(obj), depth + 1)
    if isinstance(obj, dict):
        return "{" + ",".join(_deep(k, depth + 1) + ":" + _deep(v, depth + 1) for k, v in obj.items()) + "}"
    if isinstance(obj, (list, tuple, set)):
        return "[" + ",".join(_deep(v, depth + 1) for v in obj) + "]"
    return repr(obj)


def leak_check(ows):
    """Decode every pickled / literal setting in the saved .ows (including the values inside
    Context objects, and the raw pickle bytes) and make sure no local absolute path survived
    (the repository is public). The scheme title/description are checked as raw text."""
    txt = open(ows, encoding="utf-8").read()
    blobs = re.findall(r'<properties node_id="[^"]+" format="(\w+)">(.*?)</properties>', txt, re.S)
    found = []
    for fmt, body in blobs:
        if fmt == "pickle":
            raw = base64.b64decode(body.strip())
            dec = _deep(pickle.loads(raw)) + raw.decode("latin-1")
        else:
            dec = body
        for b in BAD:
            if b in dec:
                found.append((fmt, b, dec[max(0, dec.find(b) - 60): dec.find(b) + 60]))
    for b in BAD:
        if b in txt:
            found.append(("raw", b, ""))
    print("leak check:", "OK, %d property blobs decoded (incl. Context values), no local path" % len(blobs)
          if not found else found)
    if found:
        raise SystemExit("local path found in " + ows)


# ---- verification helpers (used by verify_ch01/02/03_ows.py) ---------------------------
def load_ows(ows, timeout=120):
    """Load a saved .ows headless the way Orange.canvas.run does and let signals propagate."""
    ows = os.path.abspath(ows)
    reg, model = owsutil.registry_and_model(os.path.dirname(ows))
    from orangecanvas.scheme import signalmanager
    sm = model.findChild(signalmanager.SignalManager)
    sm.pause()
    with open(ows, "rb") as fh:
        model.load_from(fh, registry=reg)
    for n in model.nodes:
        model.widget_for_node(n)
    sm.resume()
    print("settled:", owsutil.settle(model, timeout=timeout))
    return model


def score_row(model, title="Test and Score"):
    """{column header: value} of the first learner row of a Test and Score widget."""
    w = owsutil.by_title(model, title)[1]
    hdr, rows = owsutil.test_score_table(w)
    return {h: v for h, v in zip(hdr, rows[0])}


def set_resampling(model, title, value):
    w = owsutil.by_title(model, title)[1]
    w.resampling = value
    w._param_changed() if hasattr(w, "_param_changed") else w.kfold_changed()
    owsutil.settle(model, timeout=60)


def shot(model, title, path, size=None, hide_legend=False):
    w = owsutil.by_title(model, title)[1]
    w.resize(w.sizeHint().expandedTo(w.minimumSizeHint()))
    if hide_legend:  # legend box can sit on top of data points; hide it for the picture only
        from AnyQt.QtWidgets import QCheckBox
        for cb in w.findChildren(QCheckBox):
            if cb.text() == "Show legend" and cb.isChecked():
                cb.click()
    w.show()
    if size:
        w.resize(*size)
    owsutil.settle(model, timeout=10, quiet=0.5)
    w.grab().save(path); w.hide()
    print("saved", os.path.basename(path))
