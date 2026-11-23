# 第 8 章「參考方法的誤差」用的數字：寫成 site/nir-pathlength/lib/ch08data.js（確定性計算，無亂數）
# Atago 讀值：ATAGO PAL-Easy ACID 181（導電度式酸度計）；PAL-30S 是折射計，直接顯示酸度 %（11/17 量，使用者提供：2.8、4.8、10.4）。
# 滴定資料（另一組學生，11/03）：取 10 mL 醋稀釋到 100 mL，再取 25 mL，以 NaOH 滴定；NaOH 由稱重配製（0.4315 g / 40 / 0.1 L → f = 1.07875，未標定）。
# 醋酸當量 0.006 g/mL（0.1 N × 60.05 g/mol ÷ 1000）。酸度 (w/v %) = V × f × 0.006 × (100/25) × (100/10)。
# 器材允差採 A 級玻璃器皿的典型值（示意），標準不確定度 = 允差 / √3（矩形分布）。
import os, json
import numpy as np, pandas as pd
from scipy.signal import savgol_filter

here = os.path.dirname(os.path.abspath(__file__)); root = os.path.join(here, '..', '..')
out = {}

# ---- 1. 三種參考值
V = np.array([8.8, 17.0, 44.0]); F = 0.4315 / 40 / 0.1 * 10           # 1.07875
K = F * 0.006 * (100 / 25) * (100 / 10)                                # 每 mL NaOH 的酸度 %
titr = V * K
label = np.array([3.0, 5.0, 10.0]); pal = np.array([2.8, 4.8, 10.4]); at1 = np.array([3.32, 4.87, 7.45]); at2 = np.array([3.01, 4.12, 6.86])
out['refs'] = {'label': label.tolist(), 'pal30s': pal.tolist(), 'atago_1102': at1.tolist(), 'atago_1103': at2.tolist(), 'titr': titr.tolist(),
               'V': V.tolist(), 'f': F, 'K': K, 'ratio_titr_atago1103': (titr / at2).tolist(),
               'atago_day_diff_pct': ((at1 - at2) / at1 * 100).tolist()}

# ---- 2. 計算式拆解
out['calc'] = {'eq_g_per_mL': 0.006, 'eq_exact': 0.1 * 60.05 / 1000,   # 講義用 0.006（差 0.08 %）；結果是 w/v %
                'aliquot_factor': 100 / 25, 'dilution_factor': 100 / 10, 'K': K}

# ---- 3. 誤差預算（單次滴定，相對標準不確定度 %）
u_read = 0.02                       # 每次讀滴定管 ±0.02 mL（估到 0.01–0.02）
tol = {'buret': 0.05, 'endpoint': 0.05, 'pip10': 0.02, 'flask100': 0.08, 'pip25': 0.03}
u = {k: v / np.sqrt(3) for k, v in tol.items()}
uV = float(np.sqrt(2 * u_read ** 2 + u['buret'] ** 2 + u['endpoint'] ** 2))   # 滴定體積（兩次讀數 + 允差 + 終點一滴）
budget = []
for v, c in zip(V, titr):
    rel = {'titre': uV / v * 100, 'pip10': u['pip10'] / 10 * 100, 'flask100': u['flask100'] / 100 * 100, 'pip25': u['pip25'] / 25 * 100}
    tot = float(np.sqrt(sum(r ** 2 for r in rel.values())))
    budget.append({'V': float(v), 'c': float(c), 'rel': {k: float(r) for k, r in rel.items()}, 'total_rel': tot, 'total_abs': float(c * tot / 100)})
out['budget'] = {'uV': uV, 'u_read': u_read, 'tol': tol, 'rows': budget,
                 'khp_f_rel': 0.2,                  # 以 KHP 標定、三重複後力價的相對不確定度（典型值，示意）
                 'naoh_bias_range': [2, 5]}         # 未標定 NaOH 吸濕、吸 CO2 → 真正力價比算的低 2–5 %（示意）→ 算出的酸度同比例偏高

# ---- 4. 三種參考值彼此差多少？換成滴定體積（mL）看，和讀數不確定度 uV 比
gaps = {}
for name, ref in (('label', label), ('atago_1103', at2), ('atago_1102', at1), ('pal30s', pal)):
    need = ref / K
    gaps[name] = {'need_mL': need.tolist(), 'gap_mL': (need - V).tolist(), 'gap_in_uV': ((need - V) / uV).tolist()}
out['agree'] = gaps

# NIR 的旁證（很弱）：1 mm、2258 nm，用水＋5 %＋10 % 瓶（滴定值）建 3 點線，預測 3 % 瓶，附 Fieller 95 % 區間（自由度 1）
from scipy import stats
d = pd.read_csv(os.path.join(root, 'data', 'spectra_long.csv')); d = d[(d.replicate == 1) & (d.path_mm == 1.0)]
g = d.pivot_table(index='acid_pct', columns='wavelength_nm', values='absorbance').sort_index()
wl = g.columns.values; D = savgol_filter(g.values, 9, 2, deriv=2, axis=1); j = int(np.argmin(abs(wl - 2258))); sig = D[:, j] * 1e3
ref4 = np.array([0.0, *titr]); keep = [0, 2, 3]
xs, ys = ref4[keep], sig[keep]; n = 3; b, a = np.polyfit(xs, ys, 1); s = float(np.sqrt(((ys - a - b * xs) ** 2).sum() / (n - 2)))
sxx = ((xs - xs.mean()) ** 2).sum(); tq = stats.t.ppf(.975, n - 2); kk = tq * tq * s * s; y0 = sig[1]
gq = float(kk / (b * b * sxx))
A = b * b - kk / sxx; B = -2 * b * (y0 - a) + 2 * kk * xs.mean() / sxx; C = (y0 - a) ** 2 - kk * (1 + 1 / n) - kk * xs.mean() ** 2 / sxx
lo = hi = None
if A > 0:
    r = np.sqrt(B * B - 4 * A * C); lo, hi = sorted([float((-B - r) / (2 * A)), float((-B + r) / (2 * A))])
fit = {}
for name, ref in (('label', label), ('atago_1102', at1), ('atago_1103', at2), ('titr', titr), ('pal30s', pal)):
    x = np.array([0.0, *ref]); e = []
    for i in range(4):
        mm = np.arange(4) != i; bb, aa = np.polyfit(x[mm], sig[mm], 1); e.append(float((sig[i] - aa) / bb - x[i]))
    fit[name] = {'r2': float(np.corrcoef(x, sig)[0, 1] ** 2), 'rmsecv': float(np.sqrt(np.mean(np.square(e))))}
out['nirfit'] = fit
out['nir3'] = {'pred': float((y0 - a) / b), 'g': gq, 'lo': lo, 'hi': hi, 'df': 1}

# ---- 5. 力價的比例誤差：參考值全部乘 (1 + e)，R² 不變，斜率與預測同比例偏
e = -0.03
def r2(x, y): return float(np.corrcoef(x, y)[0, 1] ** 2)
b0, a0 = np.polyfit(ref4, sig, 1); b1, a1 = np.polyfit(ref4 * (1 + e), sig, 1)
y0 = 4.20
out['factor'] = {'e': e, 'r2': r2(ref4, sig), 'r2_scaled': r2(ref4 * (1 + e), sig), 'slope': float(b0), 'slope_scaled': float(b1),
                 'pred': float((y0 - a0) / b0), 'pred_scaled': float((y0 - a1) / b1), 'y0': y0,
                 'khp_g_for_20mL': 0.1 * 0.020 * 204.22, 'khp_M': 204.22}

# ---- 6. 滴定體積的選擇：改變稀釋（D mL 醋 → 100 mL），仍取 25 mL 滴定；目標 15–40 mL（50 mL 滴定管）
def titre(c, D_):
    return c / (F * 0.006 * (100 / 25) * (100 / D_))
concs = [3, 5, 7, 11.39]; dil = [5, 10, 20]
out['aliquot'] = {'conc': concs, 'D': dil, 'V': [[float(titre(c, D_)) for D_ in dil] for c in concs],
                  'rel_titre': [[float(uV / titre(c, D_) * 100) for D_ in dil] for c in concs], 'target': [15, 40]}

# ---- 7. 參考值誤差和 NIR 殘差比一比（濃度 %）——注意每一項的依據不同
b0_, a0_ = np.polyfit(ref4, sig, 1); res = sig - np.polyval([b0_, a0_], ref4)
out['tonir'] = {'nir_resid_titr': float(np.sqrt((res ** 2).sum() / 2) / abs(b0_)),          # 實測（4 點、自由度 2）
                'titr_floor_abs': [r['total_abs'] for r in budget],                           # 器材允差推估的下限（假設）
                'atago_day_sd_abs': float(np.sqrt(np.mean((at1 - at2) ** 2) / 2)),          # 實測（兩天各一次）
                'titr_3pct_gap_abs': float(3.0 - titr[0]),                                   # 3 % 瓶滴定與標示的差
                'triplicate_factor': float(1 / np.sqrt(3))}

def rnd(o):
    if isinstance(o, float): return float('%.4g' % o)
    if isinstance(o, dict): return {k: rnd(v) for k, v in o.items()}
    if isinstance(o, list): return [rnd(v) for v in o]
    return o
out = rnd(out)
open(os.path.join(root, 'site', 'nir-pathlength', 'lib', 'ch08data.js'), 'w', encoding='utf-8').write(
    '// 由 tools/nir-pathlength/make_data_ch08.py 產生（參考方法的誤差）\nconst CH8 = ' + json.dumps(out, ensure_ascii=False) + ';\n')
print(json.dumps(out, ensure_ascii=False, indent=1))
