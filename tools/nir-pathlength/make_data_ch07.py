# 第 7 章「設計下一次實驗」用的數字：寫成 site/nir-pathlength/lib/ch07data.js
# 除了「隨機順序」的模擬（固定亂數種子 7）之外，都是確定性計算。
# 校正模型同第 5 章：1 mm、2258 nm、SG 9 點二次二階導數 ×10⁻³，反向校正，Fieller 區間。
import os, json
import numpy as np, pandas as pd
from scipy.signal import savgol_filter
from scipy import stats

here = os.path.dirname(os.path.abspath(__file__)); root = os.path.join(here, '..', '..')
d = pd.read_csv(os.path.join(root, 'data', 'spectra_long.csv'))
out = {}

# ---- 0. 現有的 1 mm 校正（同第 5 章）
g = d[(d.replicate == 1) & (d.path_mm == 1.0)].pivot_table(index='acid_pct', columns='wavelength_nm', values='absorbance').sort_index()
wl = g.columns.values; atago = g.index.values
D = savgol_filter(g.values, 9, 2, deriv=2, axis=1); j = int(np.argmin(abs(wl - 2258))); sig = D[:, j] * 1e3
titr = np.array([0] + [v * 1.07875 * 0.24 for v in (8.8, 17, 44)])

def fit(xs, ys):
    n = len(xs); b, a = np.polyfit(xs, ys, 1); s = float(np.sqrt(((ys - a - b * xs) ** 2).sum() / (n - 2)))
    return dict(a=float(a), b=float(b), s=s, s_conc=s / abs(b), sxx=float(((xs - xs.mean()) ** 2).sum()))

fa, ft = fit(atago, sig), fit(titr, sig)
out['now'] = {'atago': fa, 'titr': ft}

# ---- 1. 掃描雜訊 vs 校正殘差（只能估數量級）
# 唯一的一對重複是 0.25 mm 的水：二階導數在 2200–2300 nm 的差 /√2 = 單次掃描雜訊（0.25 mm 實測）。
# 1 mm 在 2258 nm 的吸光度較高；假設透光率雜訊固定，吸光度雜訊約放大 10^(ΔA) 倍 → 1 mm 的粗估上限。
rep = d[(d.path_mm == 0.25) & (d['sample'] == 'water')].pivot_table(index='replicate', columns='wavelength_nm', values='absorbance')
R = savgol_filter(rep.values, 9, 2, deriv=2, axis=1) * 1e3; win = (wl >= 2200) & (wl <= 2300)
scan_sd = float((R[1] - R[0])[win].std() / np.sqrt(2))
A1 = float(g.values[:, j].mean()); A025 = float(rep.values[0, j])
scale = float(10 ** (A1 - A025))
out['scan'] = {'sd025': scan_sd, 'conc025': scan_sd / abs(fa['b']), 'A1': A1, 'A025': A025, 'scale': scale,
               'conc1_est': scan_sd * scale / abs(fa['b']), 'resid_atago': fa['s_conc'], 'resid_titr': ft['s_conc']}

# ---- 2. 預試驗：重新填樣 k 次得到標準差 s，真正的 σ 落在 lo·s 到 hi·s 之間（95 %，卡方，自由度 k − 1）
#      現在的校正殘差自由度只有 2，等於 k = 3 那一列
out['pilot'] = [{'k': k, 'lo': float(np.sqrt((k - 1) / stats.chi2.ppf(.975, k - 1))),
                 'hi': float(np.sqrt((k - 1) / stats.chi2.ppf(.025, k - 1)))} for k in (2, 3, 5, 10, 20)]

# ---- 3. 區間寬度取決於還不知道的誤差（8 濃度 × 3 次填樣，0–7 %，斜率用滴定校正，未知樣品在 4 %）
# 寬度只和 b、s、n、Sxx、m 有關，令截距 a = 0。
def half_width(s_conc, m, xs, b=ft['b'], x0=4.0):
    n = len(xs); s = s_conc * abs(b); sxx = ((xs - xs.mean()) ** 2).sum(); t = stats.t.ppf(.975, n - 2)
    k = t * t * s * s; xb = xs.mean(); y0 = b * x0
    A = b * b - k / sxx; B = -2 * b * y0 + 2 * k * xb / sxx; C = y0 * y0 - k * (1 / m + 1 / n) - k * xb * xb / sxx
    r = np.sqrt(B * B - 4 * A * C); lo, hi = sorted([(-B - r) / (2 * A), (-B + r) / (2 * A)])
    return float((hi - lo) / 2)
X83 = np.repeat(np.arange(8.0), 3); X8 = np.arange(8.0)
SS = [0.1, 0.2, 0.4, 0.6, 1.0]
# (a) 24 點當成彼此獨立：只在配製誤差為 0 時成立
out['sens'] = {'s': SS, 'm1': [half_width(s, 1, X83) for s in SS], 'm3': [half_width(s, 3, X83) for s in SS]}
# (b) 同一個稀釋濃度的 3 次填樣共用配製／參考值誤差 s_ref：改用 8 個濃度平均（n = 8，自由度 6），
#     每點的變異 = s_ref² + s_fill²/3；未知樣品填 3 次 → s_fill²/3 + s_ref²（它的參考值誤差不進預測，但另一瓶的配製誤差會進斜率）
SF = 0.4; SR = [0.0, 0.2, 0.4]
def hw_levels(s_fill, s_ref, m=3):
    s_pt = np.sqrt(s_ref ** 2 + s_fill ** 2 / 3)             # 濃度平均點的標準差
    # 未知樣品 m 次填樣的平均：變異 s_fill²/m；用 s_pt 當共同 s 時 1/m 項要換算
    n = 8; b = ft['b']; xs = X8; s = s_pt * abs(b); sxx = ((xs - xs.mean()) ** 2).sum(); t = stats.t.ppf(.975, n - 2)
    k = t * t * s * s; xb = xs.mean(); y0 = b * 4.0; w_unk = (s_fill ** 2 / m) / s_pt ** 2
    A = b * b - k / sxx; B = -2 * b * y0 + 2 * k * xb / sxx; C = y0 * y0 - k * (w_unk + 1 / n) - k * xb * xb / sxx
    r = np.sqrt(B * B - 4 * A * C); lo, hi = sorted([(-B - r) / (2 * A), (-B + r) / (2 * A)])
    return float((hi - lo) / 2)
out['sens_ref'] = {'s_fill': SF, 's_ref': SR, 'hw': [hw_levels(SF, r) for r in SR],
                   'hw_fill6': [hw_levels(SF, r) if r == 0 else None for r in SR]}
# 填樣次數加倍（3 → 6）能不能壓過配製誤差？s_ref = 0.4 時
out['sens_ref']['more_fills'] = {'fills': [1, 3, 6, 12], 'pt_sd': [float(np.sqrt(0.4 ** 2 + SF ** 2 / f)) for f in (1, 3, 6, 12)]}

# ---- 4. 濃度怎麼排：同樣 24 次測量，斜率標準誤（相對於現在 4 點、各 1 次）
designs = {'now': atago, 'ends': np.repeat([0.0, 7.0], 12), 'three': np.repeat([0.0, 3.5, 7.0], 8), 'eight': X83}
sxx = {k: float(((v - v.mean()) ** 2).sum()) for k, v in designs.items()}
# 線性檢查（缺適性檢定）需要：濃度數 − 2 > 0，而且純誤差自由度（測量數 − 濃度數）> 0
out['spacing'] = {k: {'n': int(len(designs[k])), 'sxx': sxx[k], 'se_rel': float(np.sqrt(sxx['now'] / sxx[k])),
                      'levels': int(len(set(designs[k]))), 'lof_df': int(len(set(designs[k])) - 2),
                      'pure_df': int(len(designs[k]) - len(set(designs[k])))} for k in designs}
# 彎曲的真實反應（示意）：只量兩端時，中間的偏差看不到
cx = np.linspace(0, 7, 29); curve = cx - 0.04 * cx * (7 - cx)          # 中間最多偏 0.49 %
out['spacing']['curve'] = {'x': [float(v) for v in cx], 'y': [float(v) for v in curve], 'mid_bias': float(0.04 * 3.5 * 3.5)}

# ---- 5. 測量順序（示意：一次 24 個測量，儀器線性漂移，整段相當於 0.5 % 濃度）
# 斜率偏差 = Σ(x − x̄)·漂移 / Sxx。完全隨機排列時平均為 0，標準差 = √(S_dd / ((n − 1)·S_xx))（公式，不用亂數）。
# 分 3 個區塊、每塊 8 個濃度各一次、塊內隨機：塊和塊之間的漂移被抵銷，只剩塊內的部分。
drift_total = 0.5; n = 24; dr = drift_total * np.arange(n) / (n - 1)
Sxx = float(((X83 - X83.mean()) ** 2).sum()); Sdd = float(((dr - dr.mean()) ** 2).sum())
asc = np.sort(X83); asc_bias = float(((asc - asc.mean()) * dr).sum() / Sxx * 100)
sd_rand = float(np.sqrt(Sdd / ((n - 1) * Sxx)) * 100)
blk = [dr[i * 8:(i + 1) * 8] for i in range(3)]; sxx_b = float(((X8 - X8.mean()) ** 2).sum())
sd_block = float(np.sqrt(sum(sxx_b * ((b_ - b_.mean()) ** 2).sum() / 7 for b_ in blk)) / Sxx * 100)
# 範例順序：每塊 8 個濃度的一種排法（拉丁方式輪替，任何人都能照抄，不靠亂數）
base = [3, 6, 0, 5, 2, 7, 1, 4]
example = [base[(i + 3 * blk_) % 8] for blk_ in range(3) for i in range(8)]
out['order'] = {'drift_total': drift_total, 'asc_bias': asc_bias, 'rand_sd': sd_rand, 'rand_95': 1.96 * sd_rand,
                'block_sd': sd_block, 'block_95': 1.96 * sd_block, 'example': example,
                'example_bias': float(((np.array(example) - 3.5) * dr).sum() / Sxx * 100)}

# ---- 6. 課堂分組：8 組，每組做 4 個濃度（低中高都有）＋ 同一瓶共同對照
groups = [sorted({(gi + o) % 8 for o in (0, 1, 3, 6)}) for gi in range(8)]   # 示意：每個濃度剛好 4 組做，每組有高有低
out['class'] = {'groups': groups, 'per_level': [sum(lv in gr for gr in groups) for lv in range(8)]}

def rnd(o):
    if isinstance(o, float): return float('%.4g' % o)
    if isinstance(o, dict): return {k: rnd(v) for k, v in o.items()}
    if isinstance(o, list): return [rnd(v) for v in o]
    return o
out = rnd(out)
open(os.path.join(root, 'site', 'nir-pathlength', 'lib', 'ch07data.js'), 'w', encoding='utf-8').write(
    '// 由 tools/nir-pathlength/make_data_ch07.py 產生（設計下一次實驗）\nconst CH7 = ' + json.dumps(out, ensure_ascii=False) + ';\n')
show = {k: v for k, v in out.items()}; show['spacing'] = {k: v for k, v in out['spacing'].items() if k != 'curve'}
print(json.dumps(show, ensure_ascii=False, indent=1))
