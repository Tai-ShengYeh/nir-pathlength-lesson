"""Write exercises/orange/ch07_design.csv: the chapter-7 measurement-order table (deterministic, no random seed).

    python exercises/orange/tools/make_ch07_design.py

Two designs of 24 runs (8 concentration levels 0..7 %, each three times), stacked into 48 rows:
  - "由低到高"  : levels 0,0,0,1,1,1,...,7,7,7 in run order  (floor((run-1)/3))
  - "分3塊打亂" : the lesson's example order, CH7.order.example in site/nir-pathlength/lib/ch07data.js
                 (3 time blocks of 8 runs, every level once per block)
signal = level_pct + drift, drift = 0.5 * (run - 1) / 23  (illustrative instrument drift, the same
linear drift as tools/nir-pathlength/make_data_ch07.py; NOT measured data).
Header uses Orange's flag#name format: c = class (signal), m = meta, C = numeric, D = discrete."""
import os, json, csv
here = os.path.dirname(os.path.abspath(__file__))
orange_dir = os.path.normpath(os.path.join(here, ".."))
lesson = os.path.join(here, "..", "..", "..", "site", "nir-pathlength", "lib", "ch07data.js")
js = open(lesson, encoding="utf-8").read()
CH7 = json.loads(js[js.index("{"): js.rindex("}") + 1])
example = CH7["order"]["example"]
assert len(example) == 24 and sorted(example) == sorted(list(range(8)) * 3)

designs = [("由低到高", [(r - 1) // 3 for r in range(1, 25)]), ("分3塊打亂", example)]
rows = []
for name, levels in designs:
    for run in range(1, 25):
        lv = levels[run - 1]
        drift = 0.5 * (run - 1) / 23
        rows.append([run, lv, name, "%.6f" % (lv + drift), (run - 1) // 8 + 1, "%.6f" % drift])

out = os.path.join(orange_dir, "ch07_design.csv")
with open(out, "w", encoding="utf-8-sig", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["mC#run", "C#level_pct", "mD#order_type", "cC#signal", "mC#block", "mC#drift"])
    w.writerows(rows)
print("wrote", os.path.basename(out), len(rows), "rows")
