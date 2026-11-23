"""Shared helpers: load the Orange widget registry, wait for headless signal propagation,
and read numbers out of the PCA / Test and Score widgets. Used by build_ch04_ows.py and
verify_ch04_ows.py (both run with QT_QPA_PLATFORM=offscreen)."""
import os, time

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("QT_QPA_FONTDIR", "C:/Windows/Fonts")  # offscreen Qt has no fonts otherwise

from AnyQt.QtWidgets import QApplication
from AnyQt.QtTest import QTest

_app = None


def app():
    global _app
    _app = QApplication.instance() or QApplication([])
    return _app


def registry_and_model(basedir):
    """Same set-up as Orange.canvas.run.main (Orange's own headless runner)."""
    app()
    from orangecanvas import config
    from orangecanvas.registry import WidgetRegistry, cache
    from Orange.canvas.config import Config
    cfg = Config()
    config.set_default(cfg)
    config.init()
    reg = WidgetRegistry()
    cfg.widget_discovery(reg, cached_descriptions=cache.registry_cache()).run(cfg.widgets_entry_points())
    model = cfg.workflow_constructor()
    model.set_runtime_env("basedir", os.path.abspath(basedir))
    return reg, model


def settle(model, timeout=60.0, quiet=1.0):
    """Process events until no widget is busy and the signal manager is idle for `quiet` s."""
    from orangecanvas.scheme import signalmanager
    sm = model.findChild(signalmanager.SignalManager)
    t0 = time.time(); calm_since = None
    while time.time() - t0 < timeout:
        QTest.qWait(100)
        busy = sm.has_pending() if hasattr(sm, "has_pending") else False
        for node in model.nodes:
            w = model.widget_for_node(node)
            if getattr(w, "isInvalidated", lambda: False)() or getattr(w, "isBlocking", lambda: False)():
                busy = True
            task = getattr(w, "task", None)
            if task is not None:
                busy = True
        if busy:
            calm_since = None
        elif calm_since is None:
            calm_since = time.time()
        elif time.time() - calm_since > quiet:
            return True
    return False


def by_title(model, title):
    for n in model.nodes:
        if n.title == title:
            return n, model.widget_for_node(n)
    raise KeyError(title)


def test_score_table(w):
    m = w.score_table.model
    hdr = [m.headerData(c, 1) for c in range(m.columnCount())]   # 1 == Qt.Horizontal
    rows = [[m.data(m.index(r, c)) for c in range(m.columnCount())] for r in range(m.rowCount())]
    return hdr, rows


def messages(model):
    out = []
    for node in model.nodes:
        for msg in node.state_messages():
            if msg.contents:
                out.append((node.title, msg.severity, msg.contents))
    return out
