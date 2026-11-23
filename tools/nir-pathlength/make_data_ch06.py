# 第 6 章「資料品質與離群值」用的數字：寫成 site/nir-pathlength/lib/ch06data.js（確定性計算，無亂數）
# 判斷規則都在看結果前用物理定義：
#   水帶高度 = 1440–1460 nm 平均吸光度 − 1580–1650 nm（平坦區，同第 1 章估雜訊的區域）平均吸光度。
#   醋裡的水最多少約 7 %，所以同一光徑的四瓶，水帶高度應該在該光徑中位數的 ±10 % 以內；超出就是「可疑」。
#   1450 nm 水帶在 A ≥ 1.5 時已飽和（第 1 章），該光徑不能用這個檢查。
import os, json
import numpy as np, pandas as pd
from scipy.signal import savgol_filter
from sklearn.decomposition import PCA

here = os.path.dirname(os.path.abspath(__file__)); root = os.path.join(here, '..', '..')
d = pd.read_csv(os.path.join(root, 'data', 'spectra_long.csv'))
r1 = d[d.replicate == 1]
w = r1.pivot_table(index=['path_mm', 'acid_pct'], columns='wavelength_nm', values='absorbance').sort_index()
wl = w.columns.values; X = w.values; idx = w.index.to_frame().values; P, C = idx[:, 0], idx[:, 1]
CONC = [0.0, 3.32, 4.87, 7.45]; PATHS = [0.25, 0.5, 1.0, 2.0, 5.0]
WB = (wl >= 1440) & (wl <= 1460); FLAT = (wl >= 1580) & (wl <= 1650)
row = lambda p, c: int(np.where((P == p) & (C == c))[0][0])
out = {'wl': [round(float(v), 1) for v in wl]}

# ---- 1. 量尺：唯一的一對重複（0.25 mm 水量兩次）
rep = d[(d.path_mm == 0.25) & (d['sample'] == 'water')].pivot_table(index='replicate', columns='wavelength_nm', values='absorbance').values
diff = rep[1] - rep[0]
out['rep'] = {'a': [round(float(v), 4) for v in rep[0]], 'b': [round(float(v), 4) for v in rep[1]],
              'rms': float(np.sqrt(np.mean(diff ** 2))), 'max': float(np.abs(diff).max())}

# ---- 2. 水帶高度檢查
wb = {}
for p in PATHS:
    hs = [float(X[row(p, c), WB].mean() - X[row(p, c), FLAT].mean()) for c in CONC]
    peak = float(max(X[row(p, c), WB].max() for c in CONC)); med = float(np.median(hs))
    ok = peak < 1.5
    wb[str(p)] = {'h': hs, 'ratio': [h / med for h in hs], 'peak': peak, 'usable': ok,
                  'flag': [bool(ok and abs(h / med - 1) > 0.10) for h in hs]}
out['wb'] = wb
out['spec'] = {str(p): [[round(float(v), 4) for v in X[row(p, c)]] for c in CONC] for p in (0.25, 0.5, 1.0)}

# ---- 3. 自動離群偵測的陷阱：20 條不同光徑的原始光譜一起做 PCA（2 個成分），依 Q 殘差排名
pca = PCA(2).fit(X); S = pca.transform(X); Q = ((X - pca.mean_ - S @ pca.components_) ** 2).sum(1)
order = np.argsort(-Q)
out['pcaQ'] = {'q': [float(v) for v in Q], 'path': [float(v) for v in P], 'conc': [float(v) for v in C],
               'rank': [int(i) for i in order],
               'rank_of_flagged': {'0.25-4.87': int(np.where(order == row(0.25, 4.87))[0][0]) + 1,
                                   '0.5-3.32': int(np.where(order == row(0.5, 3.32))[0][0]) + 1},
               'q_in_path_025': [float(Q[row(0.25, c)]) for c in CONC]}

# ---- 4. 1 mm、2258 nm：大殘差不等於壞光譜
D = savgol_filter(X, 9, 2, deriv=2, axis=1); j = int(np.argmin(abs(wl - 2258)))
m1 = [row(1.0, c) for c in CONC]; y = D[m1, j] * 1e3; atago = np.array(CONC)
titr = np.array([0] + [v * 1.07875 * 0.24 for v in (8.8, 17, 44)])

def loo(ys, cs):
    e = []
    for i in range(len(cs)):
        k = np.arange(len(cs)) != i; b, a = np.polyfit(cs[k], ys[k], 1); e.append(float((ys[i] - a) / b - cs[i]))
    return float(np.sqrt(np.mean(np.square(e)))), e

def r2(ys, cs): return float(np.corrcoef(cs, ys)[0, 1] ** 2)

ra, ea = loo(y, atago); rt, et = loo(y, titr)
out['one'] = {'signal': [float(v) for v in y], 'atago': CONC, 'titr': [float(v) for v in titr],
              'loo_atago': ea, 'rmse_atago': ra, 'loo_titr': et, 'rmse_titr': rt,
              'r2_atago': r2(y, atago), 'r2_titr': r2(y, titr),
              'drop_r2': [r2(np.delete(y, i), np.delete(atago, i)) for i in range(4)],
              'wb_ratio': wb['1.0']['ratio']}

# ---- 5. 0.25 mm：把可疑光譜「修正」或刪掉，分數會怎樣（規則事先定好：訊號除以它的水帶比值）
m025 = [row(0.25, c) for c in CONC]; y25 = D[m025, j] * 1e3
fix = y25.copy(); fix[2] /= wb['0.25']['ratio'][2]
out['q25'] = {'signal': [float(v) for v in y25], 'fixed': [float(v) for v in fix],
              'rmse': loo(y25, atago)[0], 'rmse_fixed': loo(fix, atago)[0],
              'rmse_drop': loo(np.delete(y25, 2), np.delete(atago, 2))[0],
              'r2': r2(y25, atago), 'r2_fixed': r2(fix, atago), 'r2_drop': r2(np.delete(y25, 2), np.delete(atago, 2))}

def rnd(o):
    if isinstance(o, float): return float('%.4g' % o)
    if isinstance(o, dict): return {k: rnd(v) for k, v in o.items()}
    if isinstance(o, list): return [rnd(v) for v in o]
    return o
out = rnd(out)
open(os.path.join(root, 'site', 'nir-pathlength', 'lib', 'ch06data.js'), 'w', encoding='utf-8').write(
    '// 由 tools/nir-pathlength/make_data_ch06.py 產生（資料品質與離群值）\nconst CH6 = ' + json.dumps(out, ensure_ascii=False) + ';\n')
show = {k: v for k, v in out.items() if k not in ('wl', 'spec', 'rep')}
show['rep'] = {k: out['rep'][k] for k in ('rms', 'max')}
print(json.dumps(show, ensure_ascii=False, indent=1))
