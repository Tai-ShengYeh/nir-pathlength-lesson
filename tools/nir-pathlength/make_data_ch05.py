# 第 5 章「預測未知樣品與不確定度」用的數字：寫成 site/nir-pathlength/lib/ch05data.js（確定性計算，無亂數）
# 單一波長 2258 nm、1 mm、二階導數（×10⁻³，同第 1 章），反向校正：訊號 y = a + b·c，預測 c0 = (y0 − a)/b。
# 區間用 Fieller（精確）：滿足 |y0 − a − b·c| ≤ t·s·√(1/m + 1/n + (c − c̄)²/Sxx) 的 c；g = t²s²/(b²Sxx) ≥ 1 時區間無界。
import os, json
import numpy as np, pandas as pd
from scipy.signal import savgol_filter
from scipy import stats

here = os.path.dirname(os.path.abspath(__file__)); root = os.path.join(here, '..', '..')
d = pd.read_csv(os.path.join(root, 'data', 'spectra_long.csv')); d = d[(d.replicate == 1) & (d.path_mm == 1.0)]
g = d.pivot_table(index='acid_pct', columns='wavelength_nm', values='absorbance').sort_index()
wl = g.columns.values; atago = g.index.values
D = savgol_filter(g.values, 9, 2, deriv=2, axis=1); j = int(np.argmin(abs(wl - 2258)))
sig = D[:, j] * 1e3
titr = np.array([0] + [v * 1.07875 * 0.24 for v in (8.8, 17, 44)])
Y0 = 4.20                                                   # 假想未知樣品的訊號（×10⁻³）

def calib(xs, ys):
    n = len(xs); b, a = np.polyfit(xs, ys, 1); r = ys - (a + b * xs)
    s = float(np.sqrt((r ** 2).sum() / (n - 2))); sxx = float(((xs - xs.mean()) ** 2).sum())
    t = float(stats.t.ppf(.975, n - 2)); seb = s / np.sqrt(sxx)
    return dict(n=n, a=float(a), b=float(b), s=s, sxx=sxx, xbar=float(xs.mean()), ybar=float(ys.mean()), t=t,
                g=float(t ** 2 * s ** 2 / (b ** 2 * sxx)), tb=float(abs(b) / seb))

def predict(c, y0, m=1):
    x0 = (y0 - c['a']) / c['b']
    simple = c['t'] * c['s'] / abs(c['b']) * np.sqrt(1 / m + 1 / c['n'] + (y0 - c['ybar']) ** 2 / (c['b'] ** 2 * c['sxx']))
    # Fieller：解 (y0 − a − b·x)² = t²s²(1/m + 1/n + (x − x̄)²/Sxx) 的二次方程
    k = c['t'] ** 2 * c['s'] ** 2
    A = c['b'] ** 2 - k / c['sxx']; B = -2 * c['b'] * (y0 - c['a']) + 2 * k * c['xbar'] / c['sxx']
    C = (y0 - c['a']) ** 2 - k * (1 / m + 1 / c['n']) - k * c['xbar'] ** 2 / c['sxx']
    if A <= 0: lo = hi = None
    else:
        disc = B * B - 4 * A * C; r = np.sqrt(disc); lo, hi = sorted([(-B - r) / (2 * A), (-B + r) / (2 * A)])
    return dict(x0=float(x0), simple=float(simple), lo=None if lo is None else float(lo), hi=None if hi is None else float(hi))

out = {'signal': [round(float(v), 4) for v in sig], 'atago': [float(v) for v in atago], 'titr': [round(float(v), 3) for v in titr], 'y0': Y0}
# 1. 藏起第三瓶（Atago 4.87、滴定 4.40），用其他 3 點校正
tr = [0, 1, 3]
c3 = calib(atago[tr], sig[tr]); p3 = predict(c3, sig[2])
out['hide3'] = {'calib': c3, 'pred': p3}
# 2. 4 點校正（Atago），假想未知樣品 y0
c4 = calib(atago, sig); p4 = predict(c4, Y0); p4m3 = predict(c4, Y0, 3)
out['four'] = {'calib': c4, 'pred': p4, 'pred_m3': p4m3}
# 3. 換成滴定參考值（4 點）
c4t = calib(titr, sig); p4t = predict(c4t, Y0)
out['four_titr'] = {'calib': c4t, 'pred': p4t}
# 4. 殘差換算的偵測／定量下限（與第 1 章的儀器雜訊 LOD 不同）
out['lod_res'] = {'atago': [3 * c4['s'] / abs(c4['b']), 10 * c4['s'] / abs(c4['b'])], 'titr': [3 * c4t['s'] / abs(c4t['b']), 10 * c4t['s'] / abs(c4t['b'])]}
# 5. 設計示意：第 2 章 0–7 %、8 濃度 × 3 次；假設 s、b 不變（示意，偏樂觀）
def design(c, m):
    xs = np.repeat(np.arange(8.0), 3); n = len(xs); sxx = float(((xs - xs.mean()) ** 2).sum()); t = float(stats.t.ppf(.975, n - 2))
    cc = dict(c, n=n, sxx=sxx, xbar=float(xs.mean()), ybar=float(c['a'] + c['b'] * xs.mean()), t=t, g=float(t ** 2 * c['s'] ** 2 / (c['b'] ** 2 * sxx)))
    return dict(g=cc['g'], **predict(cc, Y0, m))
out['design'] = {'atago_m1': design(c4, 1), 'atago_m3': design(c4, 3), 'titr_m3': design(c4t, 3)}

def rnd(o):
    if isinstance(o, dict): return {k: rnd(v) for k, v in o.items()}
    if isinstance(o, float): return round(o, 4)
    return o
out = rnd(out)
open(os.path.join(root, 'site', 'nir-pathlength', 'lib', 'ch05data.js'), 'w', encoding='utf-8').write(
    '// 由 tools/nir-pathlength/make_data_ch05.py 產生（預測與不確定度）\nconst CH5 = ' + json.dumps(out, ensure_ascii=False) + ';\n')
print(json.dumps(out, ensure_ascii=False, indent=1))
