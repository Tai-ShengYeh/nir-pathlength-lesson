# 第 3 章「波長挑選的陷阱」用的數字：寫成 site/nir-pathlength/lib/ch03data.js
# 全部是確定性計算（沒有亂數），R / Python / 網頁三邊應該完全一致。
import os, json, itertools
import numpy as np, pandas as pd
from scipy.signal import savgol_filter

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.join(here, '..', '..')
d = pd.read_csv(os.path.join(root, 'data', 'spectra_long.csv'))
d = d[(d.path_mm == 1.0) & (d.replicate == 1)]
g = d.pivot_table(index='acid_pct', columns='wavelength_nm', values='absorbance').sort_index()
wl = g.columns.values; X = g.values; y = g.index.values            # y = 0, 3.32, 4.87, 7.45
D = savgol_filter(X, 9, 2, deriv=2, axis=1)                         # 9 點、二次、點距為 1
usable = (X < 1.5).all(axis=0)                                      # 4 條光譜都 A < 1.5

def r2(a, b): return float(np.corrcoef(a, b)[0, 1] ** 2)
def r2_all(rows, yy): return np.array([r2(D[rows, j], yy) for j in range(len(wl))])
def inv_pred(tr, k, j):                                             # 用 3 點擬合 d2 = a + b·c，反推第 k 點的濃度
    b, a = np.polyfit(y[tr], D[tr, j], 1); return float((D[k, j] - a) / b)

R = r2_all(slice(None), y)
j_all = int(R.argmax())                                             # 142 個波長全挑
Ru = np.where(usable, R, -1.0); j_use = int(Ru.argmax())            # 只在可用波段挑
j2258 = int(np.argmin(abs(wl - 2258)))

perms = sorted(set(itertools.permutations(range(4))))
shuffle = []
for p in perms:
    yp = y[list(p)]
    rp = np.where(usable, r2_all(slice(None), yp), -1.0)
    shuffle.append({'perm': list(p), 'true': list(p) == [0, 1, 2, 3],
                    'best': round(float(rp.max()), 4), 'best_wl': round(float(wl[rp.argmax()]), 1),
                    'fixed': round(r2(D[:, j2258], yp), 4)})

folds_in, e_out, e_fix = [], [], []
for k in range(4):
    tr = [i for i in range(4) if i != k]
    rr = np.where(usable, r2_all(tr, y[tr]), -1.0); o = np.argsort(-rr); j = int(o[0])
    pk = inv_pred(tr, k, j)
    folds_in.append({'left_out': float(y[k]), 'pick_wl': round(float(wl[j]), 1), 'r2': round(float(rr[j]), 5),
                     'r2_2nd': round(float(rr[o[1]]), 5), 'pred': round(pk, 3), 'err': round(pk - y[k], 3)})
    e_out.append(inv_pred(tr, k, j_use) - y[k]); e_fix.append(inv_pred(tr, k, j2258) - y[k])
rmse = lambda e: round(float(np.sqrt(np.mean(np.square(e)))), 3)

out = {
    'wl': [round(float(w), 1) for w in wl], 'usable': [bool(u) for u in usable], 'conc': [float(v) for v in y],
    'r2': [round(float(v), 4) for v in R],
    'pick_all': {'wl': round(float(wl[j_all]), 1), 'r2': round(float(R[j_all]), 4), 'A': [round(float(v), 3) for v in X[:, j_all]]},
    'pick_use': {'wl': round(float(wl[j_use]), 1), 'r2': round(float(R[j_use]), 4)},
    'fixed': {'wl': round(float(wl[j2258]), 1), 'r2': round(float(R[j2258]), 4)},
    'n_usable': int(usable.sum()), 'n_r2_over_095': int((R > 0.95).sum()),
    'shuffle': shuffle,
    'wrong_best_over_09': sum(s['best'] > 0.9 for s in shuffle if not s['true']),
    'wrong_fixed_over_09': sum(s['fixed'] > 0.9 for s in shuffle if not s['true']),
    'folds_inside': folds_in,
    'rmsecv': {'outside': rmse(e_out), 'inside': rmse([f['err'] for f in folds_in]), 'fixed': rmse(e_fix)},
    'err_outside': [round(e, 3) for e in e_out], 'err_fixed': [round(e, 3) for e in e_fix],
}
js = '// 由 tools/nir-pathlength/make_data_ch03.py 產生（1 mm、二階導數、波長挑選）\nconst CH3 = ' + json.dumps(out, ensure_ascii=False) + ';\n'
open(os.path.join(root, 'site', 'nir-pathlength', 'lib', 'ch03data.js'), 'w', encoding='utf-8').write(js)
for k in ['pick_all', 'pick_use', 'fixed', 'n_usable', 'n_r2_over_095', 'wrong_best_over_09', 'wrong_fixed_over_09', 'rmsecv', 'err_outside', 'err_fixed']:
    print(k, out[k])
for f in folds_in: print(f)
